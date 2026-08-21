#ifndef ALEPHPHIKK_H
#define ALEPHPHIKK_H

/*
  phi(1020) -> K+ K- reconstruction from the FULL baseline-selected track list:
  no primary/secondary masking, no displacement or pointing requirement (phi
  from b/c decays are displaced), no exclusive track claiming; same-charge
  pairs are reconstructed too and flagged as the combinatorial control.
  Output = a kaon sample defined on kinematics alone, so NO dE/dx quantity
  enters any selection; the daughters' dE/dx is stored, never cut on.

  Input = flipD0_copy'ed trackstates (physical charge = +sign(omega)); lengths
  cm, momenta GeV; the single VertexFitter_Tk call is made with
  rescale_cm_mm=false and its momenta rescaled once by 10, as in AlephV0New.
*/

#include <algorithm>
#include <cmath>

#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "edm4hep/TrackState.h"
#include "edm4hep/MCParticleData.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

#include "aleph_units.h"

namespace FCCAnalyses {
namespace AlephPhiKK {

using ROOT::VecOps::RVec;

// PDG 2024 central values [GeV]
constexpr double M_K   = 0.493677;   // charged kaon
constexpr double M_PHI = 1.019461;   // phi(1020); Gamma = 4.25 MeV
constexpr int    PDG_PHI = 333;
constexpr int    PDG_K   = 321;

// daughter momentum and energy in the phi rest frame
inline double pstarKK() {
  return std::sqrt(0.25 * M_PHI * M_PHI - M_K * M_K);  // 0.1269181 GeV
}
constexpr double estarKK() { return 0.5 * M_PHI; }

// Storage defaults: deliberately loose, the sample is meant to be re-cut
// offline. Near-collinear KK pairs leave the vertex unconstrained along the
// flight direction, hence the wide DPV_FID.
constexpr double M_LO = 0.98, M_HI = 1.10;   // stored KK mass window [GeV]
constexpr double PRE_MARGIN = 0.02;          // pre-fit window margin [GeV]
constexpr double CHI2_CUT = 25.;             // loose sanity only (ndf = 1)
constexpr double DPV_FID = 50.;              // |vtx-PV| storage fiducial [cm]
// applied to the PERIGEE momentum, while the stored trk p is the at-vertex one
constexpr double P_MIN_DEF = 0.3;            // per-track |p| floor [GeV]

// Working points stored as flags, evaluated on the stored quantities (post-fit
// mass, at-vertex daughter momenta). Charge-blind by design: signal =
// wp && !same_sign, control = wp && same_sign under identical cuts.
constexpr double WP_DM = 0.012, TIGHT_DM = 0.005;  // |m - m_phi| [GeV]
constexpr double WP_PDAU = 1.0;                    // daughter |p| [GeV]
constexpr double WP_DPV = 1.0;                     // |vtx-PV| [cm]
constexpr double TIGHT_SIGD0 = 0.01;               // daughter sigma(d0) [cm]

// Armenteros-Podolanski band variable for the equal-mass locus
// (alpha/alpha_max)^2 + (qT/p*)^2, centred on alpha = 0; 1 on the exact
// phi -> K+K- ellipse. For equal masses it reparametrises the mass window.
inline double phiBandEll(double alpha, double qt, double pmag) {
  const double ps = pstarKK();
  double beta = pmag / std::sqrt(pmag * pmag + M_PHI * M_PHI);
  if (beta <= 0.) return -1.;
  double amax = ps / (beta * estarKK());
  return std::sqrt(alpha * alpha / (amax * amax) + qt * qt / (ps * ps));
}

inline double invMassKK(const TVector3& p1, const TVector3& p2) {
  double e1 = std::sqrt(p1.Mag2() + M_K * M_K);
  double e2 = std::sqrt(p2.Mag2() + M_K * M_K);
  TVector3 p = p1 + p2;
  double e = e1 + e2;
  return std::sqrt(std::max(0., e * e - p.Mag2()));
}

// Momentum at the perigee (pT = kPtPerTeslaCm*Bz/|omega|); used only for the
// pre-fit mass window, every stored quantity uses the fitted momenta.
inline TVector3 perigeeMomentum(const edm4hep::TrackState& t, double Bz) {
  double om = std::abs(t.omega);
  if (om <= 0.) return TVector3(0., 0., 0.);
  double pt = AlephUnits::kPtPerTeslaCm * Bz / om;
  return TVector3(pt * std::cos(t.phi), pt * std::sin(t.phi), pt * t.tanLambda);
}

// 3D compatibility significance of two vertex positions: sqrt of the chi2 of
// d = x1 - x2 under the summed position covariances (packed lower triangle
// xx, yx, yy, zx, zy, zz). -1 when the summed covariance is not invertible.
template <typename CovA, typename CovB>
inline float vertexDistSig(const TVector3& d, const CovA& ca, const CovB& cb) {
  double C[3][3] = {
    {double(ca[0]) + cb[0], double(ca[1]) + cb[1], double(ca[3]) + cb[3]},
    {double(ca[1]) + cb[1], double(ca[2]) + cb[2], double(ca[4]) + cb[4]},
    {double(ca[3]) + cb[3], double(ca[4]) + cb[4], double(ca[5]) + cb[5]}};
  double det = C[0][0] * (C[1][1] * C[2][2] - C[1][2] * C[2][1])
             - C[0][1] * (C[1][0] * C[2][2] - C[1][2] * C[2][0])
             + C[0][2] * (C[1][0] * C[2][1] - C[1][1] * C[2][0]);
  if (!(std::abs(det) > 0.) || !std::isfinite(det)) return -1.;
  double inv[3][3];
  inv[0][0] = (C[1][1] * C[2][2] - C[1][2] * C[2][1]) / det;
  inv[0][1] = (C[0][2] * C[2][1] - C[0][1] * C[2][2]) / det;
  inv[0][2] = (C[0][1] * C[1][2] - C[0][2] * C[1][1]) / det;
  inv[1][0] = (C[1][2] * C[2][0] - C[1][0] * C[2][2]) / det;
  inv[1][1] = (C[0][0] * C[2][2] - C[0][2] * C[2][0]) / det;
  inv[1][2] = (C[0][2] * C[1][0] - C[0][0] * C[1][2]) / det;
  inv[2][0] = (C[1][0] * C[2][1] - C[1][1] * C[2][0]) / det;
  inv[2][1] = (C[0][1] * C[2][0] - C[0][0] * C[2][1]) / det;
  inv[2][2] = (C[0][0] * C[1][1] - C[0][1] * C[1][0]) / det;
  double dv[3] = {d.X(), d.Y(), d.Z()};
  double s2 = 0.;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) s2 += dv[i] * inv[i][j] * dv[j];
  return (s2 > 0. && std::isfinite(s2)) ? std::sqrt(s2) : -1.;
}

// ---------------------------------------------------------------------------
// Candidate container. One entry per accepted track pair; nothing is claimed
// exclusively, so a track may appear in several candidates.
// ---------------------------------------------------------------------------
struct PhiKKCands {
  RVec<float> invM;       // KK invariant mass at the fitted vertex [GeV]
  RVec<float> p, px, py, pz;  // pair momentum at the fitted vertex [GeV]
  RVec<float> alpha, qt;  // Armenteros-Podolanski (qt in GeV)
  RVec<float> bandEll;    // equal-mass AP band variable (1 = exact locus)
  RVec<float> chi2;       // vertex fit chi2, normalised (ndf = 1)
  RVec<float> vx, vy, vz; // fitted vertex position [cm]
  RVec<float> dpv;        // |vertex - PV| [cm]
  RVec<float> dpvSig;     // 3D significance of the same, -1 if undefined
  RVec<int>   same_sign;  // 1 = same-charge pair (combinatorial control)
  // per-daughter block; index 1 is the higher-|p| track of the pair
  RVec<int>   trk1_origIdx, trk2_origIdx;  // index into the original Tracks
  RVec<int>   trk1_q, trk2_q;              // physical charge, +sign(omega)
  RVec<float> trk1_p, trk2_p;              // |p| at the fitted vertex [GeV]
  RVec<float> trk1_costheta, trk2_costheta;// pz/|p| at the fitted vertex
  RVec<float> trk1_d0, trk2_d0;            // perigee d0 [cm]
  RVec<float> trk1_z0, trk2_z0;            // perigee z0 [cm]
  RVec<float> trk1_sigd0, trk2_sigd0;      // sqrt(cov[0]) [cm]
  RVec<int>   trk1_nvdet, trk2_nvdet;      // VDET hits
  RVec<int>   trk1_nitc,  trk2_nitc;       // ITC hits
  RVec<float> trk1_chi2ndf, trk2_chi2ndf;  // track fit chi2/ndf
  RVec<int>   trk1_isprim, trk2_isprim;    // 1 = in the fitted primary set
  RVec<int>   wp, tight;                   // working-point flags, charge-blind
                                           // (see the constants above)
};

// ---------------------------------------------------------------------------
// The finder.
//   tracks     : flipD0_copy'ed baseline-selected trackstates
//   orig_idx   : original-Tracks index of each entry of `tracks`
//   nvdet/nitc : VDET/ITC hit counts of each entry (-1 = unknown)
//   chi2ndf    : track fit chi2/ndf of each entry (-1 = unknown)
//   isprim     : 1 if the entry is in the fitted primary set (stored, never cut on)
//   PV         : fitted primary vertex (positions numerically cm)
//   veto_orig  : original-Tracks indices to exclude from pairing (opt-in);
//                empty by default, i.e. nothing restricts the pairing
// A cut argument <= 0 disables that cut (except the mass window, which is the
// storage definition).
// ---------------------------------------------------------------------------
inline PhiKKCands findPhiKK(
    const RVec<edm4hep::TrackState>& tracks,
    const RVec<int>& orig_idx,
    const RVec<int>& nvdet,
    const RVec<int>& nitc,
    const RVec<float>& chi2ndf,
    const RVec<int>& isprim,
    const VertexingUtils::FCCAnalysesVertex& PV,
    double solenoidBz = 1.5,
    double m_lo = M_LO, double m_hi = M_HI,   // stored KK mass window [GeV]
    double chi2_cut = CHI2_CUT,               // vertex chi2 (ndf=1); loose sanity
    double ap_band = -1.,                     // |bandEll-1| cut; <=0 = off
    double dpv_max = -1.,                     // |vtx-PV| cut [cm]; <=0 = off
    double dpv_fid = DPV_FID,                 // |vtx-PV| storage fiducial [cm]; <=0 = off
    double dpv_sig_max = -1.,                 // significance cut; <=0 = off
    double sigd0_max = -1.,                   // per-track sigma(d0) cap [cm]; <=0 = off
    int    min_vdet_itc = 0,                  // per-track min (nVDET + nITC)
    double trk_chi2ndf_max = -1.,             // per-track chi2/ndf cap; <=0 = off
    double p_min = P_MIN_DEF,                 // per-track |p| floor [GeV]; <=0 = off
    bool   do_same_sign = true,               // also reconstruct same-charge pairs
    double pre_margin = PRE_MARGIN,           // pre-fit mass window margin [GeV]
    const RVec<int>& veto_orig = {}) {        // original-Tracks indices to exclude (opt-in)

  PhiKKCands out;
  const int nTr = tracks.size();
  if (nTr < 2) return out;

  // track-quality prefilter: one pass, then pairs are formed among survivors
  std::vector<int> good;
  good.reserve(nTr);
  std::vector<TVector3> pper(nTr);
  for (int i = 0; i < nTr; ++i) {
    const auto& t = tracks[i];
    if (!veto_orig.empty()) {
      int o = (i < (int)orig_idx.size()) ? orig_idx[i] : -1;
      if (o >= 0 && std::find(veto_orig.begin(), veto_orig.end(), o) != veto_orig.end())
        continue;
    }
    if (sigd0_max > 0.) {
      double c0 = t.covMatrix[0];
      if (!(c0 > 0.) || std::sqrt(c0) > sigd0_max) continue;
    }
    if (min_vdet_itc > 0) {
      int nv = (i < (int)nvdet.size()) ? nvdet[i] : -1;
      int ni = (i < (int)nitc.size()) ? nitc[i] : -1;
      if (nv < 0 || ni < 0 || nv + ni < min_vdet_itc) continue;
    }
    if (trk_chi2ndf_max > 0.) {
      float c = (i < (int)chi2ndf.size()) ? chi2ndf[i] : -1.f;
      if (!(c >= 0.) || c > trk_chi2ndf_max) continue;
    }
    pper[i] = perigeeMomentum(t, solenoidBz);
    if (p_min > 0. && pper[i].Mag() < p_min) continue;
    good.push_back(i);
  }
  if (good.size() < 2) return out;

  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  const double pre_lo = m_lo - pre_margin, pre_hi = m_hi + pre_margin;

  RVec<edm4hep::TrackState> tr_pair(2);
  for (size_t a = 0; a + 1 < good.size(); ++a) {
    const int i = good[a];
    for (size_t b = a + 1; b < good.size(); ++b) {
      const int j = good[b];
      const bool ss = (tracks[i].omega * tracks[j].omega > 0);
      if (ss && !do_same_sign) continue;

      // pre-fit KK mass from the perigee momenta: removes the bulk of the pair
      // combinatorics before any fit, at the cost of a small tail at the low
      // mass edge (the fit can still move pT by more than the margin).
      double m_pre = invMassKK(pper[i], pper[j]);
      if (m_pre < pre_lo || m_pre > pre_hi) continue;

      tr_pair[0] = tracks[i];
      tr_pair[1] = tracks[j];
      auto v = VertexFitterSimple::VertexFitter_Tk(
          0, tr_pair, tracks, false, 0., 0., 0., 0., 0., 0., solenoidBz, false);
      if (v.updated_track_momentum_at_vertex.size() != 2) continue;
      // cm-as-mm homothety: rescale ONCE at the source (see AlephV0New)
      for (auto& tp : v.updated_track_momentum_at_vertex) tp *= 10.;

      double chi2 = v.vertex.chi2;
      if (!(chi2 == chi2)) continue;
      if (chi2_cut > 0. && chi2 >= chi2_cut) continue;

      TVector3 pa = v.updated_track_momentum_at_vertex[0];
      TVector3 pb = v.updated_track_momentum_at_vertex[1];
      double m = invMassKK(pa, pb);
      if (m < m_lo || m > m_hi) continue;

      TVector3 p = pa + pb;
      double pmag = p.Mag();
      if (!(pmag > 0.)) continue;
      double qt = pa.Cross(p.Unit()).Mag();
      double la = pa.Dot(p) / pmag, lb = pb.Dot(p) / pmag;
      // physical charge = +sign(omega) for the flipD0_copy'ed collection
      int qi = (tracks[i].omega > 0) ? 1 : -1;
      // alpha is (pL+ - pL-)/(pL+ + pL-); for a same-sign pair the labels are
      // conventional, so the ordering by charge is replaced by track order
      double lplus = ss ? la : ((qi > 0) ? la : lb);
      double lminus = ss ? lb : ((qi > 0) ? lb : la);
      double alpha = (lplus + lminus != 0.) ? (lplus - lminus) / (lplus + lminus) : 0.;
      double ell = phiBandEll(alpha, qt, pmag);
      if (ap_band > 0. && !(std::abs(ell - 1.) < ap_band)) continue;

      TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
      TVector3 d = x - pv;
      double dis = d.Mag();
      if (dpv_fid > 0. && dis > dpv_fid) continue;
      if (dpv_max > 0. && dis > dpv_max) continue;
      float dsig = vertexDistSig(d, v.vertex.covMatrix, PV.vertex.covMatrix);
      if (dpv_sig_max > 0. && !(dsig >= 0. && dsig < dpv_sig_max)) continue;

      // daughter 1 = higher-|p| track, so the per-daughter branches have a
      // definite meaning without needing the charge to disambiguate
      int i1 = i, i2 = j;
      TVector3 q1 = pa, q2 = pb;
      if (pb.Mag() > pa.Mag()) { i1 = j; i2 = i; q1 = pb; q2 = pa; }
      auto oidx = [&](int k) {
        return (k >= 0 && k < (int)orig_idx.size()) ? orig_idx[k] : -1;
      };
      auto nv = [&](const RVec<int>& src, int k) {
        return (k >= 0 && k < (int)src.size()) ? src[k] : -1;
      };
      auto c2 = [&](int k) {
        return (k >= 0 && k < (int)chi2ndf.size()) ? chi2ndf[k] : -1.f;
      };

      out.invM.push_back(m);
      out.p.push_back(pmag);
      out.px.push_back(p.X()); out.py.push_back(p.Y()); out.pz.push_back(p.Z());
      out.alpha.push_back(alpha);
      out.qt.push_back(qt);
      out.bandEll.push_back(ell);
      out.chi2.push_back(chi2);
      out.vx.push_back(x.X()); out.vy.push_back(x.Y()); out.vz.push_back(x.Z());
      out.dpv.push_back(dis);
      out.dpvSig.push_back(dsig);
      out.same_sign.push_back(ss ? 1 : 0);
      out.trk1_origIdx.push_back(oidx(i1));
      out.trk2_origIdx.push_back(oidx(i2));
      out.trk1_q.push_back((tracks[i1].omega > 0) ? 1 : -1);
      out.trk2_q.push_back((tracks[i2].omega > 0) ? 1 : -1);
      out.trk1_p.push_back(q1.Mag());
      out.trk2_p.push_back(q2.Mag());
      out.trk1_costheta.push_back(q1.Mag() > 0. ? q1.Z() / q1.Mag() : -99.);
      out.trk2_costheta.push_back(q2.Mag() > 0. ? q2.Z() / q2.Mag() : -99.);
      out.trk1_d0.push_back(tracks[i1].D0);
      out.trk2_d0.push_back(tracks[i2].D0);
      out.trk1_z0.push_back(tracks[i1].Z0);
      out.trk2_z0.push_back(tracks[i2].Z0);
      out.trk1_sigd0.push_back(tracks[i1].covMatrix[0] > 0. ? std::sqrt(tracks[i1].covMatrix[0]) : -1.);
      out.trk2_sigd0.push_back(tracks[i2].covMatrix[0] > 0. ? std::sqrt(tracks[i2].covMatrix[0]) : -1.);
      out.trk1_nvdet.push_back(nv(nvdet, i1));
      out.trk2_nvdet.push_back(nv(nvdet, i2));
      out.trk1_nitc.push_back(nv(nitc, i1));
      out.trk2_nitc.push_back(nv(nitc, i2));
      out.trk1_chi2ndf.push_back(c2(i1));
      out.trk2_chi2ndf.push_back(c2(i2));
      out.trk1_isprim.push_back(nv(isprim, i1));
      out.trk2_isprim.push_back(nv(isprim, i2));

      const double dm = std::abs(m - M_PHI);
      const bool pdau = (q1.Mag() > WP_PDAU && q2.Mag() > WP_PDAU);
      const bool prompt = (dis < WP_DPV);
      const double sd1 = out.trk1_sigd0.back(), sd2 = out.trk2_sigd0.back();
      out.wp.push_back((dm < WP_DM && pdau && prompt) ? 1 : 0);
      out.tight.push_back((dm < TIGHT_DM && pdau && prompt &&
                           sd1 < TIGHT_SIGD0 && sd2 < TIGHT_SIGD0) ? 1 : 0);
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// Per-track auxiliary quantities, aligned with a selected trackstate
// collection through its original-Tracks index map.
// ---------------------------------------------------------------------------

// Component `which` (0 = VDET, 1 = ITC, 2 = TPC) of the per-track
// subdetectorHitNumbers block. -1 when the block is missing or too short.
inline RVec<int> subdetHits(const RVec<int>& orig_idx,
                            const RVec<unsigned int>& begin,
                            const RVec<unsigned int>& end,
                            const RVec<int>& values, int which) {
  RVec<int> out;
  for (int o : orig_idx) {
    int v = -1;
    if (o >= 0 && o < (int)begin.size() && o < (int)end.size()) {
      unsigned int b = begin[o], e = end[o];
      if (b + which < e && b + which < values.size()) v = values[b + which];
    }
    out.push_back(v);
  }
  return out;
}

// 1 for each entry whose original-Tracks index appears in `set_orig`; backs
// the "daughter was in the fitted primary set" flag (stored, never cut on).
inline RVec<int> flagInSet(const RVec<int>& orig_idx, const RVec<int>& set_orig) {
  RVec<int> out;
  for (int o : orig_idx)
    out.push_back((o >= 0 && std::find(set_orig.begin(), set_orig.end(), o) !=
                   set_orig.end()) ? 1 : 0);
  return out;
}

// Original-Tracks indices of the daughters of the selected V0 candidates
// (select with the tight flag, or pass an all-ones flag for every candidate).
// Only consumed by the OPT-IN pairing veto.
inline RVec<int> claimedOrigIdx(const RVec<int>& d1, const RVec<int>& d2,
                                const RVec<int>& keep) {
  RVec<int> out;
  for (size_t i = 0; i < d1.size() && i < d2.size(); ++i) {
    if (i < keep.size() && !keep[i]) continue;
    if (d1[i] >= 0) out.push_back(d1[i]);
    if (d2[i] >= 0) out.push_back(d2[i]);
  }
  return out;
}

inline RVec<float> trackChi2Ndf(const RVec<int>& orig_idx,
                                const RVec<float>& chi2,
                                const RVec<int>& ndf) {
  RVec<float> out;
  for (int o : orig_idx) {
    float v = -1.f;
    if (o >= 0 && o < (int)chi2.size() && o < (int)ndf.size() && ndf[o] != 0)
      v = chi2[o] / float(ndf[o]);
    out.push_back(v);
  }
  return out;
}

// ---------------------------------------------------------------------------
// MC truth. The MCParticles daughter relations are EMPTY here; provenance
// comes from generatorStatus = 10000*KS + line, with `line` the 1-BASED LUND
// line of the MOTHER (MCParticles entry i is LUND line i+1). Geant
// secondaries carry generatorStatus 0 and have no LUND mother.
// ---------------------------------------------------------------------------
inline int lundMotherIdx(int genStatus) {
  return (genStatus > 0) ? (genStatus % 10000) - 1 : -1;
}

// heaviest quark in a hadron PDG code: 5 (b), 4 (c), else 0
inline int heavyFlavour(int pdg) {
  int a = std::abs(pdg);
  if (a >= 1000000000) return 0;  // nucleus
  int nq3 = (a / 10) % 10, nq2 = (a / 100) % 10, nq1 = (a / 1000) % 10;
  int mx = std::max({nq1, nq2, nq3});
  if (mx >= 5) return (mx == 5) ? 5 : 0;  // top-flavoured codes not expected
  return (mx == 4) ? 4 : 0;
}

struct TruePhis {
  RVec<int>   idx;        // MCParticles index of the phi
  RVec<int>   dauPlus;    // MC index of the K+
  RVec<int>   dauMinus;   // MC index of the K-
  RVec<int>   mothPdg;    // PDG of the immediate LUND mother (0 = none)
  RVec<int>   origin;     // 5 = from b hadron, 4 = from c hadron, 0 = other
  RVec<float> p, pt, costheta, px, py, pz;
  RVec<float> vx, vy, vz; // production point [cm]
  RVec<float> dauPlus_p, dauMinus_p;
  RVec<int>   nmatched;   // daughters with >= 1 linked track (0/1/2)
};

inline TruePhis findTruePhis(const RVec<edm4hep::MCParticleData>& mc,
                             const RVec<RVec<int>>& mcToTracks) {
  TruePhis out;
  const int n = mc.size();
  if (n == 0) return out;
  // one pass: LUND mother of every entry
  std::vector<int> moth(n);
  for (int i = 0; i < n; ++i) {
    int m = lundMotherIdx(mc[i].generatorStatus);
    moth[i] = (m >= 0 && m < n) ? m : -1;
  }
  for (int i = 0; i < n; ++i) {
    if (mc[i].PDG != PDG_PHI) continue;
    int kp = -1, km = -1;
    bool extra = false;
    for (int k = 0; k < n; ++k) {
      if (moth[k] != i) continue;
      if (mc[k].PDG == PDG_K) { if (kp < 0) kp = k; else extra = true; }
      else if (mc[k].PDG == -PDG_K) { if (km < 0) km = k; else extra = true; }
      else extra = true;
    }
    if (kp < 0 || km < 0 || extra) continue;  // not phi -> K+K-

    int mp = (moth[i] >= 0) ? mc[moth[i]].PDG : 0;
    int orig = 0;
    for (int a = moth[i], guard = 0; a >= 0 && guard < 200; a = moth[a], ++guard) {
      int hf = heavyFlavour(mc[a].PDG);
      if (hf == 5) { orig = 5; break; }
      if (hf == 4 && orig == 0) orig = 4;  // keep looking: a b ancestor wins
    }
    TVector3 pv3(mc[i].momentum.x, mc[i].momentum.y, mc[i].momentum.z);
    int nm = 0;
    if (kp < (int)mcToTracks.size() && !mcToTracks[kp].empty()) ++nm;
    if (km < (int)mcToTracks.size() && !mcToTracks[km].empty()) ++nm;

    out.idx.push_back(i);
    out.dauPlus.push_back(kp);
    out.dauMinus.push_back(km);
    out.mothPdg.push_back(mp);
    out.origin.push_back(orig);
    out.p.push_back(pv3.Mag());
    out.pt.push_back(pv3.Perp());
    out.costheta.push_back(pv3.Mag() > 0. ? pv3.Z() / pv3.Mag() : 0.);
    out.px.push_back(pv3.X()); out.py.push_back(pv3.Y()); out.pz.push_back(pv3.Z());
    out.vx.push_back(mc[i].vertex.x);
    out.vy.push_back(mc[i].vertex.y);
    out.vz.push_back(mc[i].vertex.z);
    out.dauPlus_p.push_back(TVector3(mc[kp].momentum.x, mc[kp].momentum.y, mc[kp].momentum.z).Mag());
    out.dauMinus_p.push_back(TVector3(mc[km].momentum.x, mc[km].momentum.y, mc[km].momentum.z).Mag());
    out.nmatched.push_back(nm);
  }
  return out;
}

// Per-candidate truth. `cls` is 1 when BOTH daughter tracks link to kaons
// sharing the same true phi -> K+K- mother, else 0; the daughter MC PDG and
// its LUND mother PDG make the kaon-sample composition measurable offline.
struct PhiKKTruth {
  RVec<int> cls;          // 1 = true phi -> K+K-
  RVec<int> truephi_idx;  // index into the TruePhis lists, -1 if none
  RVec<int> trk1_mcpdg, trk2_mcpdg;      // PDG of the linked MC particle (0 = unlinked)
  RVec<int> trk1_mothpdg, trk2_mothpdg;  // PDG of its LUND mother (0 = none)
};

inline PhiKKTruth classifyPhiKK(const PhiKKCands& c,
                                const RVec<RVec<int>>& trackToMCs,
                                const RVec<edm4hep::MCParticleData>& mc,
                                const TruePhis& tp) {
  PhiKKTruth out;
  const int n = mc.size();
  auto firstMC = [&](int trk) {
    if (trk < 0 || trk >= (int)trackToMCs.size() || trackToMCs[trk].empty()) return -1;
    int m = trackToMCs[trk][0];
    return (m >= 0 && m < n) ? m : -1;
  };
  for (size_t k = 0; k < c.invM.size(); ++k) {
    int m1 = firstMC(c.trk1_origIdx[k]);
    int m2 = firstMC(c.trk2_origIdx[k]);
    int pdg1 = (m1 >= 0) ? mc[m1].PDG : 0;
    int pdg2 = (m2 >= 0) ? mc[m2].PDG : 0;
    int mo1 = (m1 >= 0) ? lundMotherIdx(mc[m1].generatorStatus) : -1;
    int mo2 = (m2 >= 0) ? lundMotherIdx(mc[m2].generatorStatus) : -1;
    out.trk1_mcpdg.push_back(pdg1);
    out.trk2_mcpdg.push_back(pdg2);
    out.trk1_mothpdg.push_back((mo1 >= 0 && mo1 < n) ? mc[mo1].PDG : 0);
    out.trk2_mothpdg.push_back((mo2 >= 0 && mo2 < n) ? mc[mo2].PDG : 0);
    int ti = -1;
    if (std::abs(pdg1) == PDG_K && std::abs(pdg2) == PDG_K && pdg1 == -pdg2 &&
        mo1 >= 0 && mo1 == mo2)
      for (size_t t = 0; t < tp.idx.size(); ++t)
        if (tp.idx[t] == mo1) { ti = (int)t; break; }
    out.truephi_idx.push_back(ti);
    out.cls.push_back(ti >= 0 ? 1 : 0);
  }
  return out;
}

// 1 if the true phi at position t was reconstructed by at least one candidate.
inline RVec<int> truePhiFound(const TruePhis& tp, const PhiKKTruth& info) {
  RVec<int> out(tp.idx.size(), 0);
  for (int t : info.truephi_idx)
    if (t >= 0 && t < (int)out.size()) out[t] = 1;
  return out;
}

} // namespace AlephPhiKK
} // namespace FCCAnalyses

#endif
