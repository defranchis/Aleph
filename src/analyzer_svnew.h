#ifndef ALEPH_SVNEW_ANALYZERS_H
#define ALEPH_SVNEW_ANALYZERS_H

// ---------------------------------------------------------------------------
// WP3 standalone secondary-vertex module (PROTOTYPE, 2026-07-26).
//
// Design (user-directed): V0-FIRST pipeline — the new V0 module claims its
// tracks first; this finder runs on the remaining secondary tracks. Same
// construction principles as the V0 module (analyzer_v0new.h): one consistent
// VertexFitter_Tk fit per candidate (no mixed-unit second fit — the known
// resolution defect of the reference SV chain), alltracks passed so reco_ind
// is FILLED (track-content truth matching; the old sv_* block leaves it
// empty), greedy quality-ranked track claiming.
//
// Old-finder pathologies explicitly guarded against:
//  - degenerate collinear 2-track fits far from the PV (18% of old SVs, 85%
//    truth-unmatched): sigma_L guard (projected vertex covariance along the
//    candidate momentum) + displacement window;
//  - dR_prefilter=0.8 seed geometry: no dR prefilter here at all.
//
// All cut values are PROTOTYPE engineering values, not adopted physics —
// tabled for user sign-off in the WP3 design round.
// ---------------------------------------------------------------------------

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

// The cut values are caller-supplied with no defaults here: stage1.py's argparse
// is the single source of the production values.
constexpr double SVN_MPI = 0.13957039;

// sigma along direction u from the edm4hep lower-triangular covMatrix
// (xx, yx, yy, zx, zy, zz); fit runs in the cm-as-mm homothety, positions cm.
template <typename Cov>
inline double sigmaAlong(const Cov& c, const TVector3& u) {
  double x = u.x(), y = u.y(), z = u.z();
  double var = c[0] * x * x + c[2] * y * y + c[5] * z * z +
               2. * (c[1] * x * y + c[3] * x * z + c[4] * y * z);
  return (var > 0.) ? std::sqrt(var) : 0.;
}

// ---------------------------------------------------------------------------
// findSVs: standalone SV finder on the secondary-track collection.
//   v0s / v0_tight: output of AlephV0New::findV0s + candTight flags, used to
//   MASK V0-claimed tracks before SV finding.
//   mask_mode: 0 = no masking (control twin), 1 = mask tight-claimed tracks
//   (default), 2 = mask all V0-claimed tracks (incl. loose tier).
// Returns FCCAnalysesV0 (struct reused: pdgAbs = 0, invM = N-pion mass) so the
// existing candChi2/candDxyz/candP/candCosPointing/candPointSig getters apply.
// reco_ind indexes the SECONDARY track collection (np_tracks). NOTE this is
// NOT the index space of v0n_trk1/trk2, which are ORIGINAL Tracks indices
// (classifyV0s maps them through sec2orig, analyzer_truth.h) — walk reco_ind
// through sec2origIdx before comparing against them.
// ---------------------------------------------------------------------------
inline VertexingUtils::FCCAnalysesV0 findSVs(
    const RVec<edm4hep::TrackState>& np_tracks,
    const VertexingUtils::FCCAnalysesVertex& PV,
    const VertexingUtils::FCCAnalysesV0& v0s,
    const RVec<int>& v0_tight,
    int mask_mode,
    double solenoidBz,
    double chi2_cut,
    double dis_lo, double dis_hi,
    double sigl_max,
    int max_tracks,
    double trk_chi2,
    double cos_point_min,
    int claim_mode,
    double grow_shift) {

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
    // cm-as-mm homothety: momentum magnitudes 10x too small — rescale ONCE
    // at the source (identical to the V0 module).
    for (auto& tp : v.updated_track_momentum_at_vertex) tp *= 10.;
    double chi2 = v.vertex.chi2;  // normalised
    if (!(chi2 == chi2) || chi2 >= chi2_cut) return false;
    // per-track compatibility: the global normalised chi2 dilutes as tracks
    // are added — every track must individually fit the vertex
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
  // pairok caches which pairs fit together at all — growth later only tries
  // tracks that pair-fit with a current member (speed + junk suppression).
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
        // position-stability guard: a track that DRAGS the vertex is a track
        // from a different decay (the b->c cascade sits ~0.06 cm downstream of
        // its parent b vertex).  Reject growth steps that move the fitted
        // vertex by more than grow_shift.  0 = guard off.
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
  // claim_mode 0 (default, as prototyped): best normalised chi2 first.
  // claim_mode 1: densest seed first — order by the number of tracks that
  //   pair-fit with BOTH seed members ("clique size", free: read off pairok),
  //   chi2 as tie-break.  Rationale: chi2-ascending actively prefers the
  //   degenerate 2-track junk pairs (their chi2 q10 is 0.006 vs 0.19 for
  //   truth-matched candidates), which then claim tracks before real vertices.
  // claim_mode 2: two-phase — grow EVERY seed with all tracks available, then
  //   claim the grown candidates by (ntracks desc, chi2 asc).  Removes the
  //   dependence on seed order entirely at the price of many more fits.
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

// ---------------------------------------------------------------------------
// svn-specific getters (generic candChi2/candDxyz/candP/candCosPointing/
// candPointSig come from AlephTruth/AlephV0New on the shared struct).
// ---------------------------------------------------------------------------
inline RVec<int> candNtracks(const VertexingUtils::FCCAnalysesV0& svs) {
  RVec<int> out;
  for (const auto& v : svs.vtx) out.push_back((int)v.reco_ind.size());
  return out;
}

// SV position components relative to the PV (comp 0/1/2 = x/y/z), cm.
// Needed offline: candDxyz gives only the magnitude, while truth position
// matching against the true SV needs the full displacement vector.
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

// flat candidate<->track association (variable-length track lists):
// candTrkSV[k] = candidate index of the k-th association, candTrkIdx[k] = its
// track index in the SECONDARY collection (NOT the v0n_trk1/2 space — those
// are ORIGINAL Tracks indices; map through sec2origIdx to compare).
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
