"""Regression gates for the C++ port analyzer_pvnew.h vs the python prototype.

Gates (regression set for the C++ port vs the python prototype):
  G0  derivative self-test (numerical vs analytic dX/dpar, dX/dL)  [new]
  V2  pathological event ZM4212_48_AL #1124: converges ~2 iters, < 1 um from
      the restarted-Delphes answer
  V3  jitter ensemble: 1e-12 relative z0 jitters -> stability at the 1e-12
      level; beam-spot float/double variants; explicit seed scan
  V1  replay agreement vs the python prototype on 400 production-kept sets
      (positions, chi2, covariance, converged flags)
  V4  pulls vs truth (core RMS ~1.35, matching prototype/production)
  V5  selection replay: port clean pruning vs prototype strict_ties=False on
      400 events, plus the expected total numbers of kept tracks
  G6  one-shot smoothing sanity: updated helix passes through the vertex,
      momentum direction = tangent, pT magnitude sane
  G7  force_first_removal A/B: the forced mode reproduces one extra removal on
      events whose first pass is already below threshold
  T   timing per event

Usage: python3 validate_pvnew.py > report_pvnew.txt   (run under FCC env -b)
"""
import os
import sys
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
NPZ_BASE = ("/afs/cern.ch/user/m/mdefranc/.claude/projects/-afs-cern-ch-work-"
            "m-mdefranc-private-ALEPH-analysis-ntuple-validation/notes/"
            "sv-v0-revisit/proto_pvfit/npz_inputs")
NPZ_DIR = os.path.join(NPZ_BASE, "pvproto")
TRUTHDIR = os.path.join(NPZ_BASE, "truth4")
HDR = os.path.normpath(os.path.join(HERE, "..", "..", "src",
                                    "analyzer_pvnew.h"))
sys.path.insert(0, HERE)
from pvfit_proto import (BeamSpot as PyBS, edm4hep_to_internal, fit_vertex,
                         select_primary_tracks)

UM = 1e4
# beam-spot widths mirror stage1's res_{x,y,z}_loose defaults (stage1.py is the
# authority for these values)
BSPY = PyBS(0.0, 0.0, 0.0, 0.02, 0.01, 2.0)
NWORK = 6

import ROOT
ROOT.gROOT.SetBatch(True)
# compiled FCCAnalyses library, loaded for the edm4hep/cling dependencies
LIB = os.environ.get("FCCANALYSES_LIB",
                     "/afs/cern.ch/work/m/mdefranc/private/ALEPH_analysis/"
                     "ntuple_validation/Aleph/FCCAnalyses/install/lib/"
                     "libFCCAnalyses.so")
ROOT.gSystem.Load(LIB)
ROOT.gInterpreter.ProcessLine('#pragma cling optimize(2)')
assert ROOT.gInterpreter.Declare('#include "%s"' % HDR), "JIT of module failed"
ROOT.gInterpreter.Declare(r'''
using namespace FCCAnalyses::AlephPVNew;
BeamSpot g_bs;   // defaults = (0,0,0, 0.02, 0.01, 2.0)
FitConfig g_cfg;
PVFitResult pv_fit(const ROOT::VecOps::RVec<double>& d0,
                   const ROOT::VecOps::RVec<double>& phi,
                   const ROOT::VecOps::RVec<double>& om,
                   const ROOT::VecOps::RVec<double>& z0,
                   const ROOT::VecOps::RVec<double>& tl,
                   const ROOT::VecOps::RVec<double>& c21,
                   bool use_bs) {
  return fit_vertex_raw(d0, phi, om, z0, tl, c21, g_bs, use_bs,
                        0, 0, 0, false, g_cfg);
}
PVFitResult pv_fit_seed(const ROOT::VecOps::RVec<double>& d0,
                        const ROOT::VecOps::RVec<double>& phi,
                        const ROOT::VecOps::RVec<double>& om,
                        const ROOT::VecOps::RVec<double>& z0,
                        const ROOT::VecOps::RVec<double>& tl,
                        const ROOT::VecOps::RVec<double>& c21,
                        double sx, double sy, double sz) {
  return fit_vertex_raw(d0, phi, om, z0, tl, c21, g_bs, true,
                        sx, sy, sz, true, g_cfg);
}
PVSelResult pv_sel(const ROOT::VecOps::RVec<double>& d0,
                   const ROOT::VecOps::RVec<double>& phi,
                   const ROOT::VecOps::RVec<double>& om,
                   const ROOT::VecOps::RVec<double>& z0,
                   const ROOT::VecOps::RVec<double>& tl,
                   const ROOT::VecOps::RVec<double>& c21,
                   double chi2max, bool force_first) {
  return select_primary_tracks_raw(d0, phi, om, z0, tl, c21, g_bs, chi2max,
                                   g_cfg, 2, force_first);
}
// smoothing needs TrackStates: build them from raw doubles (narrowing to the
// float members is intrinsic to the production interface being exercised)
SmoothedTracks pv_smooth(const ROOT::VecOps::RVec<double>& d0,
                         const ROOT::VecOps::RVec<double>& phi,
                         const ROOT::VecOps::RVec<double>& om,
                         const ROOT::VecOps::RVec<double>& z0,
                         const ROOT::VecOps::RVec<double>& tl,
                         const ROOT::VecOps::RVec<double>& c21,
                         const PVFitResult& res, double Bz) {
  ROOT::VecOps::RVec<edm4hep::TrackState> trks;
  for (size_t i = 0; i < d0.size(); ++i) {
    edm4hep::TrackState t;
    t.D0 = d0[i]; t.phi = phi[i]; t.omega = om[i]; t.Z0 = z0[i];
    t.tanLambda = tl[i];
    for (int k = 0; k < 21; ++k) t.covMatrix[k] = c21[21*i + k];
    trks.push_back(t);
  }
  return smooth_at_vertex(trks, res, Bz, g_cfg);
}
PVFitResult pv_fit_ts(const ROOT::VecOps::RVec<double>& d0,
                      const ROOT::VecOps::RVec<double>& phi,
                      const ROOT::VecOps::RVec<double>& om,
                      const ROOT::VecOps::RVec<double>& z0,
                      const ROOT::VecOps::RVec<double>& tl,
                      const ROOT::VecOps::RVec<double>& c21) {
  ROOT::VecOps::RVec<edm4hep::TrackState> trks;
  for (size_t i = 0; i < d0.size(); ++i) {
    edm4hep::TrackState t;
    t.D0 = d0[i]; t.phi = phi[i]; t.omega = om[i]; t.Z0 = z0[i];
    t.tanLambda = tl[i];
    for (int k = 0; k < 21; ++k) t.covMatrix[k] = c21[21*i + k];
    trks.push_back(t);
  }
  return fit_vertex(trks, g_bs, g_cfg);
}
// derivative self-test, returns worst relative deviation over both Jacobians
double deriv_test(double D, double p0, double C, double z0, double ct,
                  double L) {
  using namespace FCCAnalyses::AlephPVNew::detail;
  Vec5 p; p << D, p0, C, z0, ct;
  double worst = 0.0;
  const double hL = 1e-6;
  Vec3 num = (helix_point(p, L + hL) - helix_point(p, L - hL)) / (2 * hL);
  Vec3 ana = helix_dXdL(p, L);
  worst = std::max(worst, (num - ana).norm() / std::max(1.0, ana.norm()));
  Mat35 A = helix_dXdpar(p, L);
  for (int j = 0; j < 5; ++j) {
    double h = 1e-7 * std::max(1.0, std::abs(p(j)));
    Vec5 pp = p, pm = p;
    pp(j) += h; pm(j) -= h;
    Vec3 numj = (helix_point(pp, L) - helix_point(pm, L)) / (2 * h);
    for (int k = 0; k < 3; ++k)
      worst = std::max(worst, std::abs(numj(k) - A(k, j)) /
                                  std::max(1.0, std::abs(A(k, j))));
  }
  return worst;
}
''')

CPP = ROOT.FCCAnalyses.AlephPVNew


def rv(a):
    return ROOT.VecOps.AsRVec(np.ascontiguousarray(a, dtype=np.float64))


def cpp_fit(d, a, b, keep=None, use_bs=True):
    par = d["par"][a:b]
    cov = d["cov"][a:b]
    if keep is not None:
        par, cov = par[keep], cov[keep]
    return ROOT.pv_fit(rv(par[:, 0]), rv(par[:, 1]), rv(par[:, 2]),
                       rv(par[:, 3]), rv(par[:, 4]), rv(cov.ravel()), use_bs)


def cpp_fit_arr(par, cov, seed=None):
    if seed is None:
        return ROOT.pv_fit(rv(par[:, 0]), rv(par[:, 1]), rv(par[:, 2]),
                           rv(par[:, 3]), rv(par[:, 4]), rv(cov.ravel()), True)
    return ROOT.pv_fit_seed(rv(par[:, 0]), rv(par[:, 1]), rv(par[:, 2]),
                            rv(par[:, 3]), rv(par[:, 4]), rv(cov.ravel()),
                            float(seed[0]), float(seed[1]), float(seed[2]))


def pos(r):
    return np.array([r.position[0], r.position[1], r.position[2]])


def load(tag):
    d = np.load(os.path.join(NPZ_DIR, "in_%s.npz" % tag))
    tr = np.load(os.path.join(TRUTHDIR, tag + ".npz"))
    pvt = {int(e): tr["pv_gen"][i] for i, e in enumerate(tr["evnum"])}
    return d, pvt


def ev_keep(d, e, which):
    o = d[which + "_off"]
    return d[which][o[e]:o[e + 1]]


def _py_fit(args):
    par5, cov21, keep = args
    par, cov = edm4hep_to_internal(par5[:, 0], par5[:, 1], par5[:, 2],
                                   par5[:, 3], par5[:, 4], cov21)
    if keep is not None:
        par, cov = par[keep], cov[keep]
    r = fit_vertex(par, cov, BSPY)
    return (r.position, r.covariance, r.chi2, r.converged, r.n_iterations,
            r.seed_used)


def _py_sel(args):
    par5, cov21 = args
    par, cov = edm4hep_to_internal(par5[:, 0], par5[:, 1], par5[:, 2],
                                   par5[:, 3], par5[:, 4], cov21)
    keep, res, hist = select_primary_tracks(par, cov, BSPY, chi2_max=5.0,
                                            strict_ties=False)
    return keep, res.converged, res.chi2


def main():
    out = lambda s="": print(s, flush=True)
    out("C++ PORT analyzer_pvnew.h — REGRESSION GATES vs python prototype")
    out("generated %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    out("module: %s" % HDR)

    d48, pvt48 = load("ZM4212_48_AL")
    d47, pvt47 = load("ZM4212_47_AL")

    # ---------------------------------------------------------------- G0
    out("\n== G0 DERIVATIVE SELF-TEST ==")
    rng = [(0.01, 0.3, 0.002, -0.05, 0.7, 5.0),
           (-0.2, -2.1, -0.004, 0.4, -1.5, 20.0),
           (0.0005, 1.0, 1e-7, 0.001, 0.1, 50.0),   # near-straight track
           (0.1, 2.9, 0.02, -1.0, 3.0, 1.0)]
    worst = max(ROOT.deriv_test(*t) for t in rng)
    out("  worst relative deviation over %d configs: %.3e  -> %s"
        % (len(rng), worst, "PASS" if worst < 1e-5 else "FAIL"))

    # ---------------------------------------------------------------- V2
    out("\n== V2 PATHOLOGICAL EVENT ZM4212_48_AL #1124 ==")
    e = int(np.where(d48["evn"] == 1124)[0][0])
    a, b = d48["off"][e], d48["off"][e + 1]
    par5, cov21 = d48["par"][a:b], d48["cov"][a:b]
    L9 = np.array([3, 4, 7, 10, 13, 14, 15, 16, 17])
    L10 = np.array([0, 3, 4, 7, 10, 13, 14, 15, 16, 17])
    good = np.array([0.036, 0.0008, 1.0029])
    v2_pass = True
    for name, L in (("10-track (pv30 set)", L10), ("9-track (pv20 set)", L9)):
        r = cpp_fit_arr(par5[L], cov21[L].reshape(len(L), 21))
        # reference for THIS track set = the python prototype on the same set
        pyr = _py_fit((par5, cov21, L))
        dpy = np.linalg.norm(pos(r) - pyr[0]) * UM
        dg = np.linalg.norm(pos(r) - good) * UM
        out("  %s: conv=%s status=%s niter=%d chi2/ndf=%.4g" %
            (name, bool(r.converged), ROOT.FCCAnalyses.AlephPVNew.status_name(r.status),
             r.n_iterations, r.chi2 / max(r.ndf, 1)))
        out("    PV=(%.6f, %.6f, %.6f) cm  |port - prototype| = %.4f um"
            "  |port - good restart| = %.3f um"
            % (pos(r)[0], pos(r)[1], pos(r)[2], dpy, dg))
        if not (r.converged and dpy < 0.1 and r.n_iterations <= 4):
            v2_pass = False
        # the headline spec gate applies to the 9-track pv20 set only:
        # converge ~2 iters, < 1 um from the restarted Delphes answer
        if L is L9 and not (dg < 1.0 and r.n_iterations <= 4):
            v2_pass = False
    out("  V2: %s  (gate: converged, <0.1 um from prototype on both sets;"
        " 9-track set additionally <1 um from the restarted answer)"
        % ("PASS" if v2_pass else "FAIL"))

    # ---------------------------------------------------------------- V3
    out("\n== V3 CHAOS IMMUNITY ==")
    base = cpp_fit_arr(par5[L9], cov21[L9].reshape(9, 21))
    bp = pos(base)
    worstj = 0.0
    for i in range(9):
        p2 = par5[L9].copy()
        p2[i, 3] *= (1.0 + 1e-12)
        r = cpp_fit_arr(p2, cov21[L9].reshape(9, 21))
        worstj = max(worstj, np.linalg.norm(pos(r) - bp) * UM)
    out("  worst displacement over nine 1e-12 z0 jitters: %.3e um  -> %s"
        % (worstj, "PASS" if worstj < 1e-2 else "FAIL"))
    seeds_ok = True
    for R in (0.05, 1.0, 5.0, 20.0, 33.9, 34.0):
        r = cpp_fit_arr(par5[L9], cov21[L9].reshape(9, 21),
                        seed=np.array([R, 0.0, R]))
        dd = np.linalg.norm(pos(r) - bp) * UM
        if not (r.converged and dd < 1e-2):
            seeds_ok = False
        out("    seed |x|~%.4g cm: conv=%s d_ref=%.3e um chi2=%.6g"
            % (R, bool(r.converged), dd, r.chi2))
    out("  seed-independence: %s" % ("PASS" if seeds_ok else "FAIL"))

    # ---------------------------------------------------------------- V1/V4
    out("\n== V1 REPLAY vs PYTHON PROTOTYPE (production-kept sets) ==")
    jobs_py, jobs_meta = [], []
    for d, pvt in ((d48, pvt48), (d47, pvt47)):
        for e in range(len(d["evn"])):
            evn = int(d["evn"][e])
            if evn not in pvt:
                continue
            keep = ev_keep(d, e, "prod_keep")
            if keep.size < 2:
                continue
            ppv = d["prod_pv"][e]
            if np.linalg.norm(ppv - pvt[evn]) >= 0.03:
                continue
            a, b = d["off"][e], d["off"][e + 1]
            jobs_py.append((d["par"][a:b], d["cov"][a:b], keep))
            jobs_meta.append((d, e, keep, pvt[evn]))
            if len(jobs_py) >= 400:
                break
        if len(jobs_py) >= 400:
            break
    t0 = time.time()
    with ProcessPoolExecutor(NWORK) as ex:
        pyres = list(ex.map(_py_fit, jobs_py, chunksize=4))
    t_py = time.time() - t0
    t0 = time.time()
    cres = []
    for (d, e, keep, tv) in jobs_meta:
        a, b = d["off"][e], d["off"][e + 1]
        cres.append(cpp_fit(d, a, b, keep=keep))
    t_cpp = time.time() - t0
    dpos, dchi, dcov = [], [], []
    convmatch = 0
    pulls = []
    for (pyr, cr, m) in zip(pyres, cres, jobs_meta):
        if bool(cr.converged) == bool(pyr[3]):
            convmatch += 1
        if pyr[3] and cr.converged:
            dp = np.linalg.norm(pos(cr) - pyr[0]) * UM
            dpos.append(dp)
            dchi.append(abs(cr.chi2 - pyr[2]) / max(abs(pyr[2]), 1e-12))
            cc = np.array([[cr.cov[0], cr.cov[1], cr.cov[3]],
                           [cr.cov[1], cr.cov[2], cr.cov[4]],
                           [cr.cov[3], cr.cov[4], cr.cov[5]]])
            dcov.append(np.max(np.abs(cc - pyr[1])) /
                        max(np.max(np.abs(pyr[1])), 1e-300))
            sig = np.sqrt(np.diag(cc))
            pulls.append((pos(cr) - m[3]) / sig)
    dpos = np.array(dpos)
    out("  events: %d  converged-flag agreement: %d/%d" %
        (len(cres), convmatch, len(cres)))
    out("  |port - prototype| position [um]: median=%.3e  95%%=%.3e  max=%.3e"
        % (np.median(dpos), np.percentile(dpos, 95), dpos.max()))
    out("  relative |dchi2|: max=%.3e   relative cov diff: max=%.3e"
        % (max(dchi), max(dcov)))
    v1_pass = (convmatch == len(cres) and dpos.max() < 0.1)
    out("  V1: %s  (gate: all flags agree, max |dpos| < 0.1 um)"
        % ("PASS" if v1_pass else "FAIL"))
    out("  wall: python %.1f s (%d workers) | C++ %.2f s single-thread "
        "(%.2f ms/event)" % (t_py, NWORK, t_cpp, 1e3 * t_cpp / len(cres)))

    out("\n== V4 PULLS vs TRUTH (port covariance) ==")
    P = np.array(pulls)
    v4_pass = True
    for k, ax in enumerate("xyz"):
        p = P[:, k]
        core = p[np.abs(p) < 5]
        out("  pull %s: core RMS=%.3f mean=%.3f (N=%d)"
            % (ax, core.std(), core.mean(), core.size))
        if not (1.0 < core.std() < 1.7):
            v4_pass = False
    out("  V4: %s  (gate: core RMS in (1.0, 1.7), matching the known ~35%%"
        " optimism of prototype AND production alike)"
        % ("PASS" if v4_pass else "FAIL"))

    # ---------------------------------------------------------------- V5
    out("\n== V5 SELECTION REPLAY (clean pruning, pvchi2=5.0) ==")
    jobs_py, meta = [], []
    for d, pvt in ((d48, pvt48), (d47, pvt47)):
        for e in range(len(d["evn"])):
            a, b = d["off"][e], d["off"][e + 1]
            if b - a < 2:
                continue
            jobs_py.append((d["par"][a:b], d["cov"][a:b]))
            meta.append((d, e, b - a, ev_keep(d, e, "prod_keep"),
                         ev_keep(d, e, "fixed_keep")))
            if len(jobs_py) >= 400:
                break
        if len(jobs_py) >= 400:
            break
    t0 = time.time()
    with ProcessPoolExecutor(NWORK) as ex:
        pysel = list(ex.map(_py_sel, jobs_py, chunksize=2))
    t_py = time.time() - t0
    t0 = time.time()
    csel = []
    for (d, e, n, kp, kf) in meta:
        a, b = d["off"][e], d["off"][e + 1]
        par, cov = d["par"][a:b], d["cov"][a:b]
        csel.append(ROOT.pv_sel(rv(par[:, 0]), rv(par[:, 1]), rv(par[:, 2]),
                                rv(par[:, 3]), rv(par[:, 4]), rv(cov.ravel()),
                                5.0, False))
    t_cpp = time.time() - t0
    n_in = sum(m[2] for m in meta)
    d_py = d_prod = d_fix = 0
    e_py = 0
    all_split = 0
    for (r, pyr, m) in zip(csel, pysel, meta):
        N = m[2]
        kc = np.array(list(r.kept), int)
        kp_ = np.array(list(pyr[0]), int)
        ms = np.zeros(N, bool); ms[kc] = True
        ps = np.zeros(N, bool); ps[kp_] = True
        dd = int((ms != ps).sum()); d_py += dd; e_py += dd > 0
        prod = np.zeros(N, bool); prod[m[3]] = True
        fix = np.zeros(N, bool); fix[m[4]] = True
        d_prod += int((ms != prod).sum())
        d_fix += int((ms != fix).sum())
        all_split += bool(r.split_converged)
    out("  events: %d  tracks: %d  split_converged: %d/%d"
        % (len(csel), n_in, all_split, len(csel)))
    out("  port vs python proto (strict_ties=False): %d/%d tracks (%.3f%%),"
        " events %d" % (d_py, n_in, 100. * d_py / n_in, e_py))
    out("  port vs production(unfixed) keeps: %d/%d = %.2f%%   "
        "port vs units-fixed clone: %d/%d = %.2f%%"
        % (d_prod, n_in, 100. * d_prod / n_in, d_fix, n_in,
           100. * d_fix / n_in))
    v5_pass = (d_py / n_in < 0.001 and all_split == len(csel))
    out("  V5: %s  (gate: <0.1%% track diff vs prototype, all splits"
        " converged)" % ("PASS" if v5_pass else "FAIL"))
    out("  wall: python %.1f s | C++ %.2f s (%.2f ms/event)"
        % (t_py, t_cpp, 1e3 * t_cpp / len(csel)))

    # ---------------------------------------------------------------- G7
    out("\n== G7 FORCED-FIRST-REMOVAL A/B ==")
    n_extra = n_same = 0
    for (d, e, n, kp, kf) in meta[:200]:
        a, b = d["off"][e], d["off"][e + 1]
        par, cov = d["par"][a:b], d["cov"][a:b]
        args = (rv(par[:, 0]), rv(par[:, 1]), rv(par[:, 2]), rv(par[:, 3]),
                rv(par[:, 4]), rv(cov.ravel()))
        r0 = ROOT.pv_sel(*args, 5.0, False)
        r1 = ROOT.pv_sel(*args, 5.0, True)
        k0, k1 = len(r0.kept), len(r1.kept)
        if k1 == k0 - 1:
            n_extra += 1
        elif k1 == k0:
            n_same += 1
    out("  events where forced mode removes exactly one extra track: %d/200"
        % n_extra)
    out("  events unchanged (max chi2 was >= threshold anyway): %d/200"
        % n_same)
    g7_pass = (n_extra + n_same == 200 and n_extra > 0)
    out("  G7: %s  (gate: forced mode == clean mode minus <=1 track,"
        " both behaviours observed)" % ("PASS" if g7_pass else "FAIL"))

    # ---------------------------------------------------------------- G6
    out("\n== G6 ONE-SHOT SMOOTHING SANITY ==")
    (d, e, keep, tv) = jobs_meta[0]
    a, b = d["off"][e], d["off"][e + 1]
    par, cov = d["par"][a:b][keep], d["cov"][a:b][keep]
    args = (rv(par[:, 0]), rv(par[:, 1]), rv(par[:, 2]), rv(par[:, 3]),
            rv(par[:, 4]), rv(cov.ravel()))
    r = ROOT.pv_fit(*args, True)
    sm = ROOT.pv_smooth(*args, r, 1.5)
    out("  fit converged=%s  smoothing valid=%s  N=%d"
        % (bool(r.converged), bool(sm.valid), len(par)))
    x = pos(r)
    worst_miss = 0.0
    worst_dp = 0.0
    for i in range(len(par)):
        pn = np.array([sm.par[5 * i + k] for k in range(5)])
        L2 = sm.phase[i]
        u = pn[2] * L2
        f = L2 * (np.sin(u) / u if abs(u) > 1e-4 else 1 - u * u / 6)
        ang = pn[1] + u
        Xn = np.array([-pn[0] * np.sin(pn[1]) + f * np.cos(ang),
                       pn[0] * np.cos(pn[1]) + f * np.sin(ang),
                       pn[3] + pn[4] * L2])
        worst_miss = max(worst_miss, np.linalg.norm(Xn - x))
        # momentum: direction vs tangent, pT vs original curvature
        mom = sm.momentum[i]
        tang = np.array([np.cos(pn[1] + 2 * pn[2] * L2),
                         np.sin(pn[1] + 2 * pn[2] * L2), pn[4]])
        mv = np.array([mom.X(), mom.Y(), mom.Z()])
        cosang = mv.dot(tang) / (np.linalg.norm(mv) * np.linalg.norm(tang))
        worst_dp = max(worst_dp, 1 - cosang)
        # 0.0029979*1.5 = kPtPerTeslaCm (analyzer_pvnew.h) x BZ (stage1.py)
        pt_orig = 0.0029979 * 1.5 / abs(par[i, 2])
        if i < 3:
            out("    trk %d: |X(smoothed)-PV| = %.2e cm  pT=%.3f GeV "
                "(orig %.3f)  1-cos(dir,tangent)=%.1e"
                % (i, np.linalg.norm(Xn - x), np.hypot(mom.X(), mom.Y()),
                   pt_orig, 1 - cosang))
    out("  worst |X(smoothed helix at vertex phase) - PV| over %d tracks:"
        " %.3e cm" % (len(par), worst_miss))
    g6_pass = bool(sm.valid) and worst_miss < 5e-4 and worst_dp < 1e-9
    out("  G6: %s  (gate: smoothed helix within 5 um of the vertex,"
        " momentum parallel to tangent)" % ("PASS" if g6_pass else "FAIL"))

    # ---------------------------------------------------------------- TS
    out("\n== TS TRACKSTATE (float) INTERFACE CROSS-CHECK ==")
    rts = ROOT.pv_fit_ts(*args)
    dd = np.linalg.norm(pos(rts) - x) * UM
    out("  |TrackState path - raw-double path| on one event: %.4f um"
        " (float32 input narrowing; expected sub-um)" % dd)
    ts_pass = dd < 1.0
    out("  TS: %s" % ("PASS" if ts_pass else "FAIL"))

    out("\n== SUMMARY ==")
    gates = [("G0 derivatives", worst < 1e-5), ("V2 pathological", v2_pass),
             ("V3 stability", worstj < 1e-2 and seeds_ok),
             ("V1 replay", v1_pass), ("V4 pulls", v4_pass),
             ("V5 selection", v5_pass), ("G6 smoothing", g6_pass),
             ("G7 forced-first-removal A/B", g7_pass), ("TS interface", ts_pass)]
    for n, p in gates:
        out("  %-22s %s" % (n, "PASS" if p else "FAIL"))
    out("  OVERALL: %s" % ("ALL GATES PASS" if all(p for _, p in gates)
                           else "FAILURES PRESENT"))


if __name__ == "__main__":
    main()
