#ifndef ALEPHV0NEW_H
#define ALEPHV0NEW_H

/*
  WP2: standalone improved V0 reconstruction (SV/V0 revisit, handoff §3).
  Keeps the reference VertexFinderLCFIPlus::get_V0s untouched for comparisons;
  same event-level output type (VertexingUtils::FCCAnalysesV0) so downstream
  plumbing (assign_V0s_to_jets, property getters) works unchanged.

  Design changes vs the reference (each motivated by a measured failure mode):
   - SINGLE consistent fit: one VertexFitter_Tk call with rescale_cm_mm=false.
     The ALEPH cm-native trackstates run through the primary fit as a
     self-consistent homothety (all lengths numerically cm); vertex position,
     chi2, phases and momentum DIRECTIONS are exact, and momentum MAGNITUDES
     from that same fit need only the exact factor 10 (cm-as-mm homothety;
     the Bz/2 factor is already applied inside the fitter). This removes the
     reference's second fit with mixed units (par mm / cov cm) - the dominant
     FP-instability amplifier - and its mixed-fit outputs.
   - SINGLE best-hypothesis assignment per pair (no Ks+Lambda+Lambdabar
     multi-booking): among windows passed, smallest normalised mass distance.
   - QUALITY-RANKED global claiming: all pairs are fitted first, candidates
     sorted by vertex chi2, tracks claimed in that order (removes the
     track-ordering dependence of the reference's greedy first-come loop).
   - Displacement window [dis_lo, dis_hi] in cm (upper bound new: junk lives
     at large dxyz), cosPointing cut, optional per-track fit-chi2 cut.
   - reco_ind is FILLED (the alltracks overload) so candidate->track-pair
     association needs no replica bookkeeping.

  Sign conventions verified 2026-07-22 (root-causer): input = flipD0_copy'ed
  trackstates, params AND covariance consistently ALEPH->LCIO transformed.
  Distance units: all fit-chain positions are numerically cm (skeptic-verified).

  Cut values are TUNABLE arguments; defaults below are starting points from
  the WP2 ROC scans - adoption of tuned values requires sign-off.
*/

#include <algorithm>
#include <cmath>
#include <numeric>

#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "edm4hep/TrackState.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

namespace FCCAnalyses {
namespace AlephV0New {

using ROOT::VecOps::RVec;

const double m_pi_ = 0.13957039;
const double m_p_  = 0.93827208;
const double MKS   = 0.497611;
const double MLAM  = 1.115683;

// ---------------------------------------------------------------------------
// Cut package: single named source. TIGHT = the adopted package (2026-07-22);
// the findV0s parameter defaults reference these, and candTight re-evaluates
// them offline. LOOSE = the ML-training superset tier (signed off 2026-07-25):
// flat pointing, widened AP bands, relaxed Lambda qT veto. Mass windows, chi2
// and displacement window are COMMON to both tiers.
// ---------------------------------------------------------------------------
constexpr double KS_M_LO = 0.40, KS_M_HI = 0.60;
constexpr double LAM_M_LO = 1.08, LAM_M_HI = 1.20;
constexpr double DIS_LO = 0.1, DIS_HI = 150.;
constexpr double CHI2_CUT = 10.;
constexpr double TIGHT_COS_KS_LOWP = 0.999, TIGHT_COS_KS_MIDP = 0.9995,
                 TIGHT_COS_KS_HIGHP = 0.9999;
constexpr double TIGHT_COS_LAM = 0.99995, TIGHT_QT_MIN_LAM = 0.04;
constexpr double AP_BAND_KS = 0.05, AP_LAM_LO = 0.10, AP_LAM_HI = 0.20;
constexpr double TIGHT_NSIG_KS_LOWP = 3., TIGHT_NSIG_KS_HIGHP = 4.;
constexpr double LOOSE_COS_POINT = 0.999;
constexpr double LOOSE_QT_MIN_LAM = 0.02;
constexpr double LOOSE_NSIG_KS = 6.;
constexpr double LOOSE_LAM_BAND_LO = 0.20, LOOSE_LAM_BAND_HI = 0.40;
constexpr double LAM_P_LO = 8., LAM_P_HI = 20.;

// Shared per-hypothesis acceptance helpers, used by findV0s (both tiers) and
// candTight (offline re-evaluation of the booked hypothesis).
inline double ksPointThr(double pmag, double lowp, double midp, double highp) {
  return (pmag < 2.) ? lowp : (pmag < 4.) ? midp : highp;
}
inline double ksBandEll(double alpha, double qt, double pmag) {
  const double PSTAR_K = 0.20582, ESTAR_K = 0.248806;
  double beta = pmag / std::sqrt(pmag * pmag + MKS * MKS);
  double amax = PSTAR_K / (beta * ESTAR_K);
  return std::sqrt(std::pow(alpha / amax, 2) + std::pow(qt / PSTAR_K, 2));
}
inline double ksBandThr(double pmag, double floor_, double nsig_lo, double nsig_hi) {
  // resolution-scaled width (signed off 2026-07-23): sigma_ell ~ 0.007+0.0015p;
  // floor_ acts as the low-p floor (bit-identical below ~9.5 GeV at nsig=3)
  double nsig = (pmag < 15.) ? nsig_lo : nsig_hi;
  return std::max(floor_, nsig * (0.007 + 0.0015 * pmag));
}
inline double lamBandEll(double alpha, double qt, double pmag) {
  const double PSTAR_L = 0.1005, ALPHA0_L = 0.69157;
  double beta = pmag / std::sqrt(pmag * pmag + MLAM * MLAM);
  double amp = 2. * PSTAR_L / (beta * MLAM);
  return std::sqrt(std::pow((std::abs(alpha) - ALPHA0_L) / amp, 2) +
                   std::pow(qt / PSTAR_L, 2));
}
inline double lamBandThr(double pmag, double lo, double hi) {
  return (pmag < LAM_P_LO) ? lo : (pmag < LAM_P_HI) ? lo + (hi - lo) * (pmag - LAM_P_LO) / (LAM_P_HI - LAM_P_LO) : hi;
}

// momenta of the two tracks at the fitted vertex (already rescaled to the true
// GeV scale inside findV0s, so every downstream consumer sees consistent values)
inline void pairMomenta(const VertexingUtils::FCCAnalysesVertex& v,
                        TVector3& p1, TVector3& p2) {
  p1 = v.updated_track_momentum_at_vertex[0];
  p2 = v.updated_track_momentum_at_vertex[1];
}

inline double invMass(const TVector3& p1, double m1, const TVector3& p2, double m2) {
  double e1 = std::sqrt(p1.Mag2() + m1 * m1);
  double e2 = std::sqrt(p2.Mag2() + m2 * m2);
  TVector3 p = p1 + p2;
  double e = e1 + e2;
  return std::sqrt(std::max(0., e * e - p.Mag2()));
}

// ---------------------------------------------------------------------------
// The finder.
//   np_tracks : flipD0_copy'ed non-primary trackstates (LCIO-consistent)
//   PV        : fitted primary vertex (positions numerically cm)
// TWO-TIER selection (2026-07-25): every pair is evaluated against the TIGHT
// (adopted) package first; only tight-failing pairs enter the LOOSE training
// tier. Tight candidates claim tracks first, so filtering the output to
// tight candidates (candTight==1) reproduces the historical tight-only module
// output EXACTLY. Returns candidates in claim order (tight block first, chi2
// ascending within each tier); pdgAbs holds the single best hypothesis (310 or
// 3122); invM the mass under that hypothesis.
// ---------------------------------------------------------------------------
inline VertexingUtils::FCCAnalysesV0 findV0s(
    const RVec<edm4hep::TrackState>& np_tracks,
    const VertexingUtils::FCCAnalysesVertex& PV,
    double solenoidBz = 1.5,
    double cos_point_ks = TIGHT_COS_KS_HIGHP,          // Ks pointing, p>=4 GeV tier; first cut arg so the stage1 --v0nKsPointing override passes ONLY this value (NB candTight assumes the default)
    double ks_m_lo = KS_M_LO, double ks_m_hi = KS_M_HI,        // Ks mass window [GeV] (sidebands kept; common to both tiers)
    double lam_m_lo = LAM_M_LO, double lam_m_hi = LAM_M_HI,    // Lambda mass window [GeV] (common to both tiers)
    double dis_lo = DIS_LO, double dis_hi = DIS_HI,            // displacement window [cm] (common)
    double cos_point_lam = TIGHT_COS_LAM,              // tight pointing cut, Lambda hyp (adopted 2026-07-22)
    double qt_min_lam = TIGHT_QT_MIN_LAM,              // tight Armenteros qT conversion veto, Lambda hyp [GeV]
    double cos_ks_lowp = TIGHT_COS_KS_LOWP,            // tight Ks pointing, p<2 GeV tier (p-dependent: flat cut crushed low-p Ks)
    double cos_ks_midp = TIGHT_COS_KS_MIDP,            // tight Ks pointing, 2<=p<4 GeV tier
    double ap_band_ks = AP_BAND_KS,                    // Ks exact-locus AP band |ell-1| floor (<=0 disables the band in BOTH tiers)
    double ap_lam_lo = AP_LAM_LO, double ap_lam_hi = AP_LAM_HI, // tight Lambda band thr(p): lo below 8 GeV, linear to hi at 20 GeV (<=0 off)
    double chi2_cut = CHI2_CUT,                        // vertex chi2 (ndf=1, common)
    double trk_chi2_cut = -1.,                         // per-track chi2 (<=0 off, common)
    bool lam_point_ks_tiers = false) {                 // sizing variant (user-approved 2026-07-26): tight Lambda pointing uses the Ks p-tiers instead of cos_point_lam; candTight still encodes the ADOPTED package

  VertexingUtils::FCCAnalysesV0 result;
  const int nTr = np_tracks.size();
  if (nTr < 2) return result;

  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);

  struct Cand {
    VertexingUtils::FCCAnalysesVertex vtx;
    int i, j;
    int pdg;       // best hypothesis
    double m;      // mass under best hypothesis
    double chi2;
    bool tight;    // passed the tight (adopted) package
  };
  std::vector<Cand> cands;

  RVec<edm4hep::TrackState> tr_pair(2);
  for (int i = 0; i < nTr - 1; ++i) {
    tr_pair[0] = np_tracks[i];
    for (int j = i + 1; j < nTr; ++j) {
      if (np_tracks[i].omega * np_tracks[j].omega > 0) continue; // same charge
      tr_pair[1] = np_tracks[j];

      // single consistent fit; alltracks passed so reco_ind is filled
      auto v = VertexFitterSimple::VertexFitter_Tk(
          0, tr_pair, np_tracks, false, 0., 0., 0., 0., 0., 0., solenoidBz, false);
      if (v.updated_track_momentum_at_vertex.size() != 2) continue;
      // cm-as-mm homothety: single-fit momentum magnitudes are 10x too small
      // (Bz/2 already applied inside the fitter); rescale ONCE at the source so
      // masses, p, qt and every ntuple branch downstream are consistent
      for (auto& tp : v.updated_track_momentum_at_vertex) tp *= 10.;
      double chi2 = v.vertex.chi2; // normalised, ndf=1
      if (chi2 >= chi2_cut || !(chi2 == chi2)) continue;
      if (trk_chi2_cut > 0 && v.reco_chi2.size() == 2 &&
          (v.reco_chi2[0] > trk_chi2_cut || v.reco_chi2[1] > trk_chi2_cut))
        continue;

      // displacement window (cm) + pointing
      TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
      TVector3 d = x - pv;
      double dis = d.Mag();
      if (dis < dis_lo || dis > dis_hi) continue;

      TVector3 p1, p2;
      pairMomenta(v, p1, p2);
      TVector3 p = p1 + p2;
      double pmag = p.Mag();
      if (pmag <= 0) continue;
      double cp = d.Dot(p) / (dis * pmag);
      double qt = p1.Cross(p.Unit()).Mag(); // Armenteros qT (same for either daughter)
      // Armenteros alpha: physical charge = -sign(omega) (flipped collection)
      double la = p1.Dot(p) / pmag, lb = p2.Dot(p) / pmag;
      double q1 = (np_tracks[i].omega < 0) ? 1. : -1.;
      double lplus = (q1 > 0) ? la : lb, lminus = (q1 > 0) ? lb : la;
      double alpha = (lplus + lminus != 0.) ? (lplus - lminus) / (lplus + lminus) : 0.;

      // hypothesis masses: Ks(pipi), Lambda(p pi) with proton = higher-|p| track
      // (in a Lambda decay the baryon carries most of the momentum)
      double mks = invMass(p1, m_pi_, p2, m_pi_);
      double mlam = (p1.Mag() > p2.Mag()) ? invMass(p1, m_p_, p2, m_pi_)
                                          : invMass(p1, m_pi_, p2, m_p_);

      // TIGHT (adopted 2026-07-22) package first. Arbitration among the
      // tight-passing hypotheses only, so the tight subset is EXACTLY the
      // historical tight-only module output.
      // Ks: p-dependent pointing tiers + resolution-scaled exact-locus AP band
      // (sigma-scaling rationale: fixed width shrank to 0.85 sigma at 30-40 GeV
      // and caused the high-p efficiency deficit; p>=15 GeV runs at 4 sigma).
      // Lambda: pointing + qT conversion veto + p-dependent exact-locus band
      // (protects the 10-20 GeV s-tagging bins).
      bool inWinKs = (mks > ks_m_lo && mks < ks_m_hi);
      bool inWinLam = (mlam > lam_m_lo && mlam < lam_m_hi);
      bool okKs = inWinKs && cp > ksPointThr(pmag, cos_ks_lowp, cos_ks_midp, cos_point_ks);
      if (okKs && ap_band_ks > 0)
        okKs = std::abs(ksBandEll(alpha, qt, pmag) - 1.) <
               ksBandThr(pmag, ap_band_ks, TIGHT_NSIG_KS_LOWP, TIGHT_NSIG_KS_HIGHP);
      double cos_lam_thr = lam_point_ks_tiers
          ? ksPointThr(pmag, cos_ks_lowp, cos_ks_midp, cos_point_ks)
          : cos_point_lam;
      bool okLam = inWinLam && cp > cos_lam_thr && qt > qt_min_lam;
      if (okLam && ap_lam_lo > 0)
        okLam = std::abs(lamBandEll(alpha, qt, pmag) - 1.) <
                lamBandThr(pmag, ap_lam_lo, ap_lam_hi);
      bool tight = okKs || okLam;
      if (!tight) {
        // LOOSE training tier (signed off 2026-07-25): flat pointing, widened
        // AP bands, relaxed Lambda qT veto; windows/chi2/displacement common.
        okKs = inWinKs && cp > LOOSE_COS_POINT;
        if (okKs && ap_band_ks > 0)
          okKs = std::abs(ksBandEll(alpha, qt, pmag) - 1.) <
                 ksBandThr(pmag, ap_band_ks, LOOSE_NSIG_KS, LOOSE_NSIG_KS);
        okLam = inWinLam && cp > LOOSE_COS_POINT && qt > LOOSE_QT_MIN_LAM;
        if (okLam && ap_lam_lo > 0)
          okLam = std::abs(lamBandEll(alpha, qt, pmag) - 1.) <
                  lamBandThr(pmag, LOOSE_LAM_BAND_LO, LOOSE_LAM_BAND_HI);
        if (!okKs && !okLam) continue;
      }
      double dks = std::abs(mks - MKS) / (0.5 * (ks_m_hi - ks_m_lo));
      double dlam = std::abs(mlam - MLAM) / (0.5 * (lam_m_hi - lam_m_lo));
      int pdg; double m;
      if (okKs && (!okLam || dks <= dlam)) { pdg = 310; m = mks; }
      else                                  { pdg = 3122; m = mlam; }

      cands.push_back({v, i, j, pdg, m, chi2, tight});
    }
  }

  // quality-ranked global claiming: tight candidates claim first (preserving
  // the historical tight-only output), then loose; best chi2 first within a tier
  std::vector<size_t> order(cands.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    if (cands[a].tight != cands[b].tight) return cands[a].tight;
    return cands[a].chi2 < cands[b].chi2;
  });

  std::vector<bool> used(nTr, false);
  for (size_t k : order) {
    const Cand& c = cands[k];
    if (used[c.i] || used[c.j]) continue;
    used[c.i] = true;
    used[c.j] = true;
    result.vtx.push_back(c.vtx);
    result.pdgAbs.push_back(c.pdg);
    result.invM.push_back(c.m);
  }
  return result;
}

// ---------------------------------------------------------------------------
// Sizing-variant entry point (user-approved 2026-07-26): identical to the
// adopted findV0s defaults except the tight-tier Lambda pointing is aligned to
// the Ks p-tiers (0.999 / 0.9995 / 0.9999 for p<2 / 2-4 / >=4 GeV) — the
// Study-B "aligned" scenario. Defaults are spelled here in C++ next to their
// constants (single source, no Python-side hand-sync). NOT for adopted
// productions. NB (skeptic-verified 2026-07-27): the v0n_tight branch in its
// output STILL encodes the ADOPTED package (candTight re-derivation) — NOT
// variant-tier membership. Variant-tight must be re-derived offline from
// kinematics; trusting v0n_tight would silently drop the tier migrants.
// ---------------------------------------------------------------------------
inline VertexingUtils::FCCAnalysesV0 findV0sLamKsPointing(
    const RVec<edm4hep::TrackState>& np_tracks,
    const VertexingUtils::FCCAnalysesVertex& PV,
    double solenoidBz = 1.5) {
  return findV0s(np_tracks, PV, solenoidBz, TIGHT_COS_KS_HIGHP, KS_M_LO, KS_M_HI,
                 LAM_M_LO, LAM_M_HI, DIS_LO, DIS_HI, TIGHT_COS_LAM,
                 TIGHT_QT_MIN_LAM, TIGHT_COS_KS_LOWP, TIGHT_COS_KS_MIDP,
                 AP_BAND_KS, AP_LAM_LO, AP_LAM_HI, CHI2_CUT, -1.,
                 /*lam_point_ks_tiers=*/true);
}

// ---------------------------------------------------------------------------
// Truth-free per-candidate diagnostics (work on data; reco_ind is filled by
// this module, momenta are already at the true GeV scale).
// ---------------------------------------------------------------------------

inline RVec<float> candAlpha(const VertexingUtils::FCCAnalysesV0& v0s,
                             const RVec<edm4hep::TrackState>& secondaries) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) {
    if (v.reco_ind.size() < 2 || v.updated_track_momentum_at_vertex.size() < 2 ||
        v.reco_ind[0] < 0 || v.reco_ind[0] >= (int)secondaries.size()) {
      out.push_back(-99.);
      continue;
    }
    TVector3 pa = v.updated_track_momentum_at_vertex[0];
    TVector3 pb = v.updated_track_momentum_at_vertex[1];
    TVector3 p = pa + pb;
    double la = pa.Dot(p) / p.Mag(), lb = pb.Dot(p) / p.Mag();
    double q1 = (secondaries[v.reco_ind[0]].omega < 0) ? 1. : -1.;
    double lp = (q1 > 0) ? la : lb, lm = (q1 > 0) ? lb : la;
    out.push_back((lp + lm != 0.) ? (lp - lm) / (lp + lm) : -99.);
  }
  return out;
}

// Offline tight-package flag: 1 if the candidate's BOOKED hypothesis passes the
// adopted (2026-07-22) tight package, 0 if it entered via the loose training
// tier. Uses the same shared helpers/constants as findV0s (single source), so
// selecting candTight==1 reproduces the historical tight-only module output
// exactly (tight candidates claim tracks first). Assumes the adopted package
// (i.e. no --v0nKsPointing override in the production). Works on data.
inline RVec<int> candTight(const VertexingUtils::FCCAnalysesV0& v0s,
                           const VertexingUtils::FCCAnalysesVertex& PV,
                           const RVec<edm4hep::TrackState>& secondaries) {
  RVec<int> out;
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  for (size_t c = 0; c < v0s.vtx.size(); ++c) {
    const auto& v = v0s.vtx[c];
    if (v.reco_ind.size() < 2 || v.updated_track_momentum_at_vertex.size() < 2 ||
        v.reco_ind[0] < 0 || v.reco_ind[0] >= (int)secondaries.size()) {
      out.push_back(0);
      continue;
    }
    TVector3 p1 = v.updated_track_momentum_at_vertex[0];
    TVector3 p2 = v.updated_track_momentum_at_vertex[1];
    TVector3 p = p1 + p2;
    double pmag = p.Mag();
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    TVector3 d = x - pv;
    double dis = d.Mag();
    if (pmag <= 0 || dis <= 0) { out.push_back(0); continue; }
    double cp = d.Dot(p) / (dis * pmag);
    double qt = p1.Cross(p.Unit()).Mag();
    double la = p1.Dot(p) / pmag, lb = p2.Dot(p) / pmag;
    double q1 = (secondaries[v.reco_ind[0]].omega < 0) ? 1. : -1.;
    double lplus = (q1 > 0) ? la : lb, lminus = (q1 > 0) ? lb : la;
    double alpha = (lplus + lminus != 0.) ? (lplus - lminus) / (lplus + lminus) : 0.;
    double m = v0s.invM[c];
    bool ok;
    if (v0s.pdgAbs[c] == 310) {
      ok = (m > KS_M_LO && m < KS_M_HI) &&
           cp > ksPointThr(pmag, TIGHT_COS_KS_LOWP, TIGHT_COS_KS_MIDP, TIGHT_COS_KS_HIGHP);
      if (ok)
        ok = std::abs(ksBandEll(alpha, qt, pmag) - 1.) <
             ksBandThr(pmag, AP_BAND_KS, TIGHT_NSIG_KS_LOWP, TIGHT_NSIG_KS_HIGHP);
    } else {
      ok = (m > LAM_M_LO && m < LAM_M_HI) && cp > TIGHT_COS_LAM && qt > TIGHT_QT_MIN_LAM;
      if (ok)
        ok = std::abs(lamBandEll(alpha, qt, pmag) - 1.) <
             lamBandThr(pmag, AP_LAM_LO, AP_LAM_HI);
    }
    out.push_back(ok ? 1 : 0);
  }
  return out;
}

// Pointing significance: chi2-like significance of the displacement component
// PERPENDICULAR to the candidate momentum, using candidate + PV position
// covariance (all in the consistent cm space). A well-pointing candidate has
// sig ~ O(1) regardless of how precisely it is measured — the basis for
// replacing the fixed cosPointing threshold with a measurement-aware cut.
// Returns -1 if the transverse covariance is singular/non-positive.
inline RVec<float> candPointSig(const VertexingUtils::FCCAnalysesV0& v0s,
                                const VertexingUtils::FCCAnalysesVertex& PV) {
  RVec<float> out;
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  const auto& cP = PV.vertex.covMatrix; // packed lower-tri: xx,xy,yy,xz,yz,zz
  for (const auto& v : v0s.vtx) {
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    TVector3 p(0., 0., 0.);
    for (const auto& tp : v.updated_track_momentum_at_vertex) p += tp;
    TVector3 d = x - pv;
    if (p.Mag() <= 0 || d.Mag() <= 0) { out.push_back(-1.); continue; }
    TVector3 ph = p.Unit();
    TVector3 u1 = ph.Orthogonal().Unit();
    TVector3 u2 = ph.Cross(u1);
    const auto& cV = v.vertex.covMatrix;
    double C[3][3] = {
      {double(cV[0]) + cP[0], double(cV[1]) + cP[1], double(cV[3]) + cP[3]},
      {double(cV[1]) + cP[1], double(cV[2]) + cP[2], double(cV[4]) + cP[4]},
      {double(cV[3]) + cP[3], double(cV[4]) + cP[4], double(cV[5]) + cP[5]}};
    auto quad = [&](const TVector3& a, const TVector3& b) {
      double s = 0.;
      double av[3] = {a.X(), a.Y(), a.Z()}, bv[3] = {b.X(), b.Y(), b.Z()};
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) s += av[i] * C[i][j] * bv[j];
      return s;
    };
    double c11 = quad(u1, u1), c22 = quad(u2, u2), c12 = quad(u1, u2);
    double det = c11 * c22 - c12 * c12;
    if (det <= 0. || c11 <= 0. || c22 <= 0.) { out.push_back(-1.); continue; }
    double d1 = d.Dot(u1), d2 = d.Dot(u2);
    double sig2 = (d1 * (c22 * d1 - c12 * d2) + d2 * (c11 * d2 - c12 * d1)) / det;
    out.push_back(sig2 > 0. ? std::sqrt(sig2) : 0.);
  }
  return out;
}

inline RVec<float> candQt(const VertexingUtils::FCCAnalysesV0& v0s) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) {
    TVector3 pa = v.updated_track_momentum_at_vertex[0];
    TVector3 p = pa + v.updated_track_momentum_at_vertex[1];
    out.push_back(pa.Cross(p.Unit()).Mag());
  }
  return out;
}

} // namespace AlephV0New
} // namespace FCCAnalyses

#endif
