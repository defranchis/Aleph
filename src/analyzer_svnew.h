#ifndef ALEPH_SVNEW_ANALYZERS_H
#define ALEPH_SVNEW_ANALYZERS_H

// Standalone secondary-vertex finder on the secondary-track collection, run
// after the V0 module has claimed its tracks (V0-first pipeline).

#include "ROOT/RVec.hxx"
#include "TVector3.h"
#include "edm4hep/TrackState.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

#include <algorithm>
#include <numeric>
#include <vector>

namespace FCCAnalyses {
namespace AlephSVNew {

using ROOT::VecOps::RVec;

constexpr double SVN_MPI = 0.13957039;

// Adopted secondary-vertex selection: the single source for these values.
constexpr double SVN_CHI2 = 10.;        // normalised vertex chi2 (seed and growth)
constexpr double SVN_DIS_LO = 0.03;     // PV displacement window, low edge [cm]
constexpr double SVN_DIS_HI = 3.;       // PV displacement window, high edge [cm]
constexpr double SVN_SIGL_MAX = 0.10;   // longitudinal vertex sigma guard [cm]
constexpr int SVN_MAX_TRK = 8;          // maximum tracks per candidate (growth cap)
constexpr double SVN_TRK_CHI2 = 5.;     // per-track chi2 contribution cap (<=0 off)
constexpr double SVN_COS_POINT = 0.7;   // minimum cosPointing
constexpr int SVN_CLAIM_MODE = 0;       // seed ordering: best-chi2 seed first
constexpr double SVN_GROW_SHIFT = 0.;   // max fitted-vertex shift per growth step [cm]; 0 = off

// V0-track masking modes for findSVs
constexpr int SVN_MASK_NONE = 0;        // mask nothing (unmasked control twin)
constexpr int SVN_MASK_MODE = 1;        // adopted: mask the tight-claimed V0 tracks

// sigma along direction u from the edm4hep lower-triangular covMatrix
// (xx, yx, yy, zx, zy, zz); fit runs in the cm-as-mm homothety, positions cm.
template <typename Cov>
inline double sigmaAlong(const Cov& c, const TVector3& u) {
  double x = u.x(), y = u.y(), z = u.z();
  double var = c[0] * x * x + c[2] * y * y + c[5] * z * z +
               2. * (c[1] * x * y + c[3] * x * z + c[4] * y * z);
  return (var > 0.) ? std::sqrt(var) : 0.;
}

// findSVs: v0s/v0_tight mask V0-claimed tracks; mask_mode 0 = none, 1 = tight
// only, 2 = all V0-claimed. Returns FCCAnalysesV0 (pdgAbs = 0, invM = N-pion
// mass); reco_ind indexes the SECONDARY collection, NOT the original Tracks
// index space of v0n_trk1/trk2 (map through sec2origIdx to compare).
inline VertexingUtils::FCCAnalysesV0 findSVs(
    const RVec<edm4hep::TrackState>& np_tracks,
    const VertexingUtils::FCCAnalysesVertex& PV,
    const VertexingUtils::FCCAnalysesV0& v0s,
    const RVec<int>& v0_tight,
    int mask_mode,
    double solenoidBz,
    double chi2_cut = SVN_CHI2,
    double dis_lo = SVN_DIS_LO, double dis_hi = SVN_DIS_HI,
    double sigl_max = SVN_SIGL_MAX,
    int max_tracks = SVN_MAX_TRK,
    double trk_chi2 = SVN_TRK_CHI2,
    double cos_point_min = SVN_COS_POINT,
    int claim_mode = SVN_CLAIM_MODE,
    double grow_shift = SVN_GROW_SHIFT) {

  VertexingUtils::FCCAnalysesV0 result;
  const int nTr = np_tracks.size();
  if (nTr < 2) return result;

  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);

  std::vector<bool> masked(nTr, false);
  if (mask_mode > 0) {
    for (size_t iv = 0; iv < v0s.vtx.size(); ++iv) {
      if (mask_mode == 1 && (iv >= v0_tight.size() || v0_tight[iv] != 1)) continue;
      for (int ti : v0s.vtx[iv].reco_ind)
        if (ti >= 0 && ti < nTr) masked[ti] = true;
    }
  }

  struct Cand {
    VertexingUtils::FCCAnalysesVertex vtx;
    std::vector<int> trk;
    double chi2;   // normalised
    double mass;
  };

  auto fitGroup = [&](const std::vector<int>& idx,
                      VertexingUtils::FCCAnalysesVertex& out) -> bool {
    RVec<edm4hep::TrackState> group;
    for (int k : idx) group.push_back(np_tracks[k]);
    auto v = VertexFitterSimple::VertexFitter_Tk(
        0, group, np_tracks, false, 0., 0., 0., 0., 0., 0., solenoidBz, false);
    if ((int)v.updated_track_momentum_at_vertex.size() != (int)idx.size()) return false;
    // cm-as-mm homothety: momentum magnitudes 10x too small — rescale once.
    for (auto& tp : v.updated_track_momentum_at_vertex) tp *= 10.;
    double chi2 = v.vertex.chi2;  // normalised
    if (!(chi2 == chi2) || chi2 >= chi2_cut) return false;
    // per-track compatibility: every track must individually fit the vertex
    if (trk_chi2 > 0 && v.reco_chi2.size() == idx.size())
      for (float rc : v.reco_chi2)
        if (rc > trk_chi2) return false;
    out = v;
    return true;
  };

  auto passWindows = [&](const VertexingUtils::FCCAnalysesVertex& v,
                         double& mass_out) -> bool {
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    TVector3 d = x - pv;
    double dis = d.Mag();
    if (dis < dis_lo || dis > dis_hi) return false;
    TVector3 psum(0., 0., 0.);
    double esum = 0.;
    for (const auto& tp : v.updated_track_momentum_at_vertex) {
      psum += tp;
      esum += std::sqrt(tp.Mag2() + SVN_MPI * SVN_MPI);
    }
    if (psum.Mag() <= 0.) return false;
    // pointing: the SV momentum must not be anti-aligned with the flight line
    if (d.Dot(psum) / (dis * psum.Mag()) < cos_point_min) return false;
    // collinear-degeneracy guard: fit constrained along the track bundle?
    if (sigmaAlong(v.vertex.covMatrix, psum.Unit()) > sigl_max) return false;
    double m2 = esum * esum - psum.Mag2();
    mass_out = (m2 > 0.) ? std::sqrt(m2) : 0.;
    return true;
  };

  // ---- seed pass: all unmasked pairs, any charge combination -------------
  // pairok caches which pairs fit together; growth only tries pair-linked tracks.
  std::vector<Cand> seeds;
  std::vector<char> pairok((size_t)nTr * nTr, 0);
  for (int i = 0; i < nTr - 1; ++i) {
    if (masked[i]) continue;
    for (int j = i + 1; j < nTr; ++j) {
      if (masked[j]) continue;
      VertexingUtils::FCCAnalysesVertex v;
      if (!fitGroup({i, j}, v)) continue;
      pairok[(size_t)i * nTr + j] = pairok[(size_t)j * nTr + i] = 1;
      double m;
      if (!passWindows(v, m)) continue;
      seeds.push_back({v, {i, j}, v.vertex.chi2, m});
    }
  }

  // ---- growth of ONE candidate ------------------------------------------
  // repeatedly attach the available (unblocked, pair-linked) track giving the
  // best refit chi2, while windows/guards still pass.
  auto growCand = [&](Cand c, const std::vector<bool>& blocked) {
    bool grew = true;
    while (grew && (int)c.trk.size() < max_tracks) {
      grew = false;
      Cand best = c;
      TVector3 xc(c.vtx.vertex.position[0], c.vtx.vertex.position[1],
                  c.vtx.vertex.position[2]);
      for (int k = 0; k < nTr; ++k) {
        if (blocked[k]) continue;
        if (std::find(c.trk.begin(), c.trk.end(), k) != c.trk.end()) continue;
        bool linked = false;
        for (int m0 : c.trk)
          if (pairok[(size_t)m0 * nTr + k]) { linked = true; break; }
        if (!linked) continue;
        std::vector<int> trial = c.trk;
        trial.push_back(k);
        VertexingUtils::FCCAnalysesVertex v;
        if (!fitGroup(trial, v)) continue;
        double m;
        if (!passWindows(v, m)) continue;
        // position-stability guard: reject growth steps moving the fitted
        // vertex by more than grow_shift (cm). 0 = guard off.
        if (grow_shift > 0.) {
          TVector3 xn(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
          if ((xn - xc).Mag() > grow_shift) continue;
        }
        if (!grew || v.vertex.chi2 < best.chi2) {
          best = {v, trial, v.vertex.chi2, m};
          grew = true;
        }
      }
      if (grew) c = best;
    }
    return c;
  };

  // ---- seed ordering -----------------------------------------------------
  // claim_mode 0 = best normalised chi2 first; 1 = densest seed first (clique
  // size from pairok, chi2 tie-break); 2 = grow every seed with all tracks
  // available, then claim by (ntracks desc, chi2 asc).
  std::vector<size_t> order(seeds.size());
  std::iota(order.begin(), order.end(), 0);
  if (claim_mode == 1) {
    std::vector<int> clq(seeds.size(), 0);
    for (size_t s = 0; s < seeds.size(); ++s) {
      int i = seeds[s].trk[0], j = seeds[s].trk[1], n = 0;
      for (int k = 0; k < nTr; ++k)
        if (!masked[k] && k != i && k != j &&
            pairok[(size_t)i * nTr + k] && pairok[(size_t)j * nTr + k]) ++n;
      clq[s] = n;
    }
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
      if (clq[a] != clq[b]) return clq[a] > clq[b];
      return seeds[a].chi2 < seeds[b].chi2;
    });
  } else {
    std::stable_sort(order.begin(), order.end(),
                     [&](size_t a, size_t b) { return seeds[a].chi2 < seeds[b].chi2; });
  }

  auto emit = [&](const Cand& c) {
    result.vtx.push_back(c.vtx);
    result.pdgAbs.push_back(0);
    result.invM.push_back(c.mass);
  };

  std::vector<bool> used(nTr, false);
  if (claim_mode <= 1) {
    std::vector<bool> blocked = masked;
    for (size_t s : order) {
      Cand c = seeds[s];
      if (used[c.trk[0]] || used[c.trk[1]]) continue;
      c = growCand(c, blocked);
      for (int t : c.trk) { used[t] = true; blocked[t] = true; }
      emit(c);
    }
  } else {
    std::vector<Cand> grown;
    grown.reserve(seeds.size());
    for (size_t s : order) grown.push_back(growCand(seeds[s], masked));
    std::vector<size_t> ord2(grown.size());
    std::iota(ord2.begin(), ord2.end(), 0);
    std::stable_sort(ord2.begin(), ord2.end(), [&](size_t a, size_t b) {
      if (grown[a].trk.size() != grown[b].trk.size())
        return grown[a].trk.size() > grown[b].trk.size();
      return grown[a].chi2 < grown[b].chi2;
    });
    for (size_t s : ord2) {
      bool clash = false;
      for (int t : grown[s].trk) if (used[t]) { clash = true; break; }
      if (clash) continue;
      for (int t : grown[s].trk) used[t] = true;
      emit(grown[s]);
    }
  }
  return result;
}

// svn-specific getters (generic candChi2/candDxyz/candP/candCosPointing/
// candPointSig come from AlephTruth/AlephV0New on the shared struct).
inline RVec<int> candNtracks(const VertexingUtils::FCCAnalysesV0& svs) {
  RVec<int> out;
  for (const auto& v : svs.vtx) out.push_back((int)v.reco_ind.size());
  return out;
}

// SV position components relative to the PV (comp 0/1/2 = x/y/z), cm.
inline RVec<float> candDcomp(const VertexingUtils::FCCAnalysesV0& svs,
                             const VertexingUtils::FCCAnalysesVertex& PV,
                             int comp) {
  RVec<float> out;
  for (const auto& v : svs.vtx)
    out.push_back((float)(v.vertex.position[comp] - PV.vertex.position[comp]));
  return out;
}

inline RVec<float> candSigL(const VertexingUtils::FCCAnalysesV0& svs) {
  RVec<float> out;
  for (const auto& v : svs.vtx) {
    TVector3 psum(0., 0., 0.);
    for (const auto& tp : v.updated_track_momentum_at_vertex) psum += tp;
    out.push_back(psum.Mag() > 0.
                      ? (float)sigmaAlong(v.vertex.covMatrix, psum.Unit())
                      : -1.f);
  }
  return out;
}

// Flat candidate<->track association: candTrkSV[k] = candidate index,
// candTrkIdx[k] = its track index in the SECONDARY collection (NOT the
// v0n_trk1/2 original-Tracks space; map through sec2origIdx to compare).
inline RVec<int> candTrkSV(const VertexingUtils::FCCAnalysesV0& svs) {
  RVec<int> out;
  for (size_t i = 0; i < svs.vtx.size(); ++i)
    for (size_t k = 0; k < svs.vtx[i].reco_ind.size(); ++k) out.push_back((int)i);
  return out;
}

inline RVec<int> candTrkIdx(const VertexingUtils::FCCAnalysesV0& svs) {
  RVec<int> out;
  for (const auto& v : svs.vtx)
    for (int t : v.reco_ind) out.push_back(t);
  return out;
}

}  // namespace AlephSVNew
}  // namespace FCCAnalyses

#endif
