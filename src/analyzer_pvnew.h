#ifndef ALEPHPVNEW_H
#define ALEPHPVNEW_H

/*
  Standalone primary-vertex (PV) fitter: damped Gauss-Newton with an optional
  Gaussian beam-spot constraint, a deterministic seed ladder, and iterative
  pvchi2 track pruning through the same fitter entry point.

  Units: cm / rad / 1/cm everywhere. Inputs are ALEPH-flipped, cm-native
  edm4hep track states. Momentum: pT [GeV] = kPtPerTeslaCm * Bz[T] / |omega|.

  Every result carries an explicit `converged` flag and a failure status; a
  non-converged fit is never returned silently as a vertex.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "edm4hep/TrackState.h"
#include "edm4hep/EDM4hepVersion.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "aleph_units.h"

namespace FCCAnalyses {
namespace AlephPVNew {

using ROOT::VecOps::RVec;

// pT [GeV] = kPtPerTeslaCm * Bz [T] / |omega [1/cm]| (single source: aleph_units.h)
constexpr double kPtPerTeslaCm = AlephUnits::kPtPerTeslaCm;
using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
using Vec5 = Eigen::Matrix<double, 5, 1>;
using Mat5 = Eigen::Matrix<double, 5, 5>;
using Mat35 = Eigen::Matrix<double, 3, 5>;

// ---------------------------------------------------------------------------
// configuration and result types
// ---------------------------------------------------------------------------

struct BeamSpot {
  // Gaussian luminous-region constraint, PHYSICAL cm. Widths are normally
  // supplied by the caller; these defaults serve the offline harness.
  double x = 0.0, y = 0.0, z = 0.0;
  double sigma_x = 0.02;   // 200 um
  double sigma_y = 0.01;   // 100 um
  double sigma_z = 2.0;    // 2 cm
  Vec3 center() const { return Vec3(x, y, z); }
  Mat3 cov_inv() const {
    Mat3 m = Mat3::Zero();
    m(0, 0) = 1.0 / (sigma_x * sigma_x);
    m(1, 1) = 1.0 / (sigma_y * sigma_y);
    m(2, 2) = 1.0 / (sigma_z * sigma_z);
    return m;
  }
};

struct FitConfig {
  int max_iter = 60;            // hard cap on damped Gauss-Newton iterations
  double tol_step = 1e-8;       // convergence: dx^T H dx (dimensionless)
  double tol_dchi2 = 1e-7;      // convergence: |chi2_new - chi2_old|, absolute
  double tol_dchi2_rel = 1e-7;  // ... plus this times the current chi2; must sit
                                // above the chi2-evaluation noise floor
  double lm_lambda0 = 1e-6;     // initial Levenberg-Marquardt damping
  double lm_up = 10.0;          // damping increase on a rejected step
  double lm_down = 0.1;         // damping decrease on an accepted step
  double lm_lambda_max = 1e12;  // give up on the seed above this
  int max_reject = 30;          // cumulative rejected steps before giving up
  int phase_iter = 6;           // inner Newton steps for the per-track phase
  double phase_tol = 1e-12;     // inner phase convergence, cm
  double rcond = 1e-12;         // eigenvalue floor (relative) in reg_inv
  double max_radius = 100.0;    // cm; a step leaving this ball HARD-fails the
                                // seed run (status diverged_radius)
};

enum PVStatus : int {
  kOk = 0,
  kMaxIter = 1,
  kLmStall = 2,
  kDivergedRadius = 3,
  kSingularHessian = 4,
  kLinalgFail = 5,
  kTooFewTracks = 6
};

inline const char* status_name(int s) {
  switch (s) {
    case kOk: return "ok";
    case kMaxIter: return "max_iter";
    case kLmStall: return "lm_stall";
    case kDivergedRadius: return "diverged_radius";
    case kSingularHessian: return "singular_hessian";
    case kLinalgFail: return "linalg_fail";
    case kTooFewTracks: return "too_few_tracks";
  }
  return "unknown";
}

enum PVSeed : int {
  kSeedNone = -1,
  kSeedUser = 0,
  kSeedLinear = 1,
  kSeedBeamspot = 2,
  kSeedPerigeeMedian = 3
};

inline const char* seed_name(int s) {
  switch (s) {
    case kSeedUser: return "user";
    case kSeedLinear: return "linear";
    case kSeedBeamspot: return "beamspot";
    case kSeedPerigeeMedian: return "perigee_median";
  }
  return "none";
}

struct PVFitResult {
  // ALWAYS fully populated. converged == false means the numbers are a best
  // effort and MUST NOT be used as a vertex.
  std::array<double, 3> position{{std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::quiet_NaN()}};
  std::array<double, 6> cov{{0, 0, 0, 0, 0, 0}};  // lower triangle:
                                                  // xx yx yy zx zy zz
  double chi2 = std::numeric_limits<double>::quiet_NaN();  // incl. BS term
  double chi2_beamspot = 0.0;
  int ndf = 0;      // 2N - 3 (production convention)
  int ndf_bs = 0;   // 2N (beam-spot dof counted, when constrained)
  int n_iterations = 0;
  bool converged = false;
  int status = kTooFewTracks;
  RVec<double> track_chi2;
  RVec<double> track_phase;   // transverse arc length L at the vertex, cm
  int n_tracks = 0;
  bool beamspot_used = false;
  double cond_hessian = std::numeric_limits<double>::infinity();
  double cond_track_max = std::numeric_limits<double>::infinity();
  int n_rejected_steps = 0;
  int seed_used = kSeedNone;
  std::string seeds_tried;    // e.g. "linear:ok" / "linear:lm_stall,beamspot:ok"

  double chi2_per_ndf() const {
    return ndf > 0 ? chi2 / ndf : std::numeric_limits<double>::quiet_NaN();
  }
  TVector3 position_v3() const {
    return TVector3(position[0], position[1], position[2]);
  }
  double cov_xx() const { return cov[0]; }
  double cov_yx() const { return cov[1]; }
  double cov_yy() const { return cov[2]; }
  double cov_zx() const { return cov[3]; }
  double cov_zy() const { return cov[4]; }
  double cov_zz() const { return cov[5]; }
};

struct PVSelResult {
  // Result of the iterative pvchi2 pruning.
  RVec<int> kept;           // indices (into the input track list) kept primary
  PVFitResult fit;          // the fit of the final kept set
  bool split_converged = false;  // EVERY pruning pass converged (fit.converged
                                 // covers only the final fit)
  bool trivial = false;  // fewer than min_tracks tracks entered the fit, so the
                         // vertex carries no event information even when both
                         // converged flags are true; check before using it.
  int n_passes = 0;
};

struct SmoothedTracks {
  // One-shot post-convergence smoothing output.
  bool valid = false;                 // false when the fit did not converge
  RVec<double> par;                   // 5 per track: (D, phi0, C, z0, ct),
                                      // updated to pass through the vertex
  RVec<TVector3> momentum;            // momentum vector at the vertex, GeV
  RVec<double> phase;                 // arc length of the vertex on the
                                      // UPDATED helix, cm
};

// ---------------------------------------------------------------------------
// input conversion: ALEPH-flipped cm-native edm4hep -> internal (D,phi0,C,z0,ct)
// ---------------------------------------------------------------------------

namespace detail {

// packed lower-triangular index of the edm4hep 21-element covMatrix,
// ordered (d0, phi, omega, z0, tanLambda)
constexpr int kTri[5][5] = {{0, 1, 3, 6, 10},
                            {1, 2, 4, 7, 11},
                            {3, 4, 5, 8, 12},
                            {6, 7, 8, 9, 13},
                            {10, 11, 12, 13, 14}};
// (d0, phi, omega, z0, tanLambda) -> (D, phi0, C, z0, ct)
constexpr double kScale[5] = {1.0, 1.0, -0.5, 1.0, 1.0};

struct TrackSet {
  std::vector<Vec5> par;
  std::vector<Mat5> cov;
  size_t size() const { return par.size(); }
  TrackSet subset(const std::vector<int>& idx) const {
    TrackSet s;
    s.par.reserve(idx.size());
    s.cov.reserve(idx.size());
    for (int i : idx) {
      s.par.push_back(par[i]);
      s.cov.push_back(cov[i]);
    }
    return s;
  }
};

inline void add_track(TrackSet& ts, double d0, double phi, double omega,
                      double z0, double tanl, const double* cov21) {
  Vec5 p;
  p << d0 * kScale[0], phi * kScale[1], omega * kScale[2], z0 * kScale[3],
      tanl * kScale[4];
  Mat5 c;
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      c(i, j) = cov21[kTri[i][j]] * kScale[i] * kScale[j];
  ts.par.push_back(p);
  ts.cov.push_back(c);
}

inline TrackSet convert(const RVec<edm4hep::TrackState>& tracks) {
  TrackSet ts;
  ts.par.reserve(tracks.size());
  ts.cov.reserve(tracks.size());
  double c21[21];
  for (const auto& t : tracks) {
    for (int k = 0; k < 21; ++k) c21[k] = static_cast<double>(t.covMatrix[k]);
    add_track(ts, static_cast<double>(t.D0), static_cast<double>(t.phi),
              static_cast<double>(t.omega), static_cast<double>(t.Z0),
              static_cast<double>(t.tanLambda), c21);
  }
  return ts;
}

// -------------------------------------------------------------------------
// helix model
// -------------------------------------------------------------------------

inline double sinc_(double u) {
  if (std::abs(u) < 1e-4) return 1.0 - u * u / 6.0 + u * u * u * u / 120.0;
  return std::sin(u) / u;
}

inline double dsinc_(double u) {
  if (std::abs(u) < 1e-4) return -u / 3.0 + u * u * u / 30.0;
  return (u * std::cos(u) - std::sin(u)) / (u * u);
}

inline Vec3 helix_point(const Vec5& p, double L) {
  const double D = p(0), p0 = p(1), C = p(2), z0 = p(3), ct = p(4);
  const double u = C * L;
  const double f = L * sinc_(u);  // sin(CL)/C
  const double ang = p0 + u;
  return Vec3(-D * std::sin(p0) + f * std::cos(ang),
              D * std::cos(p0) + f * std::sin(ang), z0 + ct * L);
}

inline Vec3 helix_dXdL(const Vec5& p, double L) {
  const double ang = p(1) + 2.0 * p(2) * L;
  return Vec3(std::cos(ang), std::sin(ang), p(4));
}

inline Mat35 helix_dXdpar(const Vec5& p, double L) {
  const double D = p(0), p0 = p(1), C = p(2), ct = p(4);
  (void)ct;
  const double u = C * L;
  const double s_u = sinc_(u);
  const double f = L * s_u;             // sin(CL)/C
  const double dfdC = L * L * dsinc_(u);  // d/dC [sin(CL)/C]
  const double ang = p0 + u;
  const double ca = std::cos(ang), sa = std::sin(ang);
  const double cp = std::cos(p0), sp = std::sin(p0);
  Mat35 A = Mat35::Zero();
  A(0, 0) = -sp;
  A(1, 0) = cp;
  A(0, 1) = -D * cp - f * sa;
  A(1, 1) = -D * sp + f * ca;
  // d/dC at fixed L: z does not depend on C
  A(0, 2) = dfdC * ca - f * L * sa;
  A(1, 2) = dfdC * sa + f * L * ca;
  A(2, 3) = 1.0;
  A(2, 4) = L;
  return A;
}

// -------------------------------------------------------------------------
// deterministic linear algebra
// -------------------------------------------------------------------------

// Symmetric 3x3 inverse via eigen-decomposition with a relative eigenvalue
// floor. Deterministic: no pivoting, no recursion, no branch on exact zeros.
inline Mat3 reg_inv(const Mat3& M_in, double rcond, double& cond, bool& ok) {
  Mat3 M = 0.5 * (M_in + M_in.transpose());
  // Iterative QR solver, NOT computeDirect(): the analytic 3x3 path loses
  // precision on the ill-conditioned per-track weight matrices.
  Eigen::SelfAdjointEigenSolver<Mat3> es(M);
  if (es.info() != Eigen::Success) {
    cond = std::numeric_limits<double>::infinity();
    ok = false;
    return Mat3::Constant(std::numeric_limits<double>::quiet_NaN());
  }
  const Vec3 w = es.eigenvalues();
  const double wmax = w.cwiseAbs().maxCoeff();
  if (!(wmax > 0) || !std::isfinite(wmax)) {
    cond = std::numeric_limits<double>::infinity();
    ok = false;
    return Mat3::Constant(std::numeric_limits<double>::quiet_NaN());
  }
  const double floor_ = rcond * wmax;
  ok = (w.array() > floor_).all();
  Vec3 wc;
  double wmin = std::numeric_limits<double>::infinity();
  for (int i = 0; i < 3; ++i) {
    wc(i) = w(i) > floor_ ? w(i) : floor_;
    wmin = std::min(wmin, wc(i));
  }
  cond = wmax / wmin;
  const Mat3& V = es.eigenvectors();
  return V * wc.cwiseInverse().asDiagonal() * V.transpose();
}

// Solve H dx = -g. Cholesky first (H is positive definite by construction);
// on failure fall back to the regularised eigen-inverse.
inline Vec3 cholesky_solve(const Mat3& H, const Vec3& g, double rcond) {
  Eigen::LLT<Mat3> llt(H);
  if (llt.info() == Eigen::Success) return llt.solve(-g);
  double cond;
  bool ok;
  Mat3 Hi = reg_inv(H, rcond, cond, ok);
  return -(Hi * g);
}

// -------------------------------------------------------------------------
// chi2 machinery
// -------------------------------------------------------------------------

struct TrackTerms {
  std::vector<Mat3> D;     // projected weights, symmetrised
  std::vector<Vec3> X;     // helix points at the final phases
  std::vector<double> chi2;
  std::vector<double> L;   // refined phases
  double cond_max = 0.0;
};

inline void track_terms(const TrackSet& ts, const Vec3& x,
                        const std::vector<double>& L_in, const FitConfig& cfg,
                        TrackTerms& out) {
  const size_t N = ts.size();
  std::vector<double> L = L_in;
  // inner Newton on the per-track phase at fixed x
  for (int it = 0; it < cfg.phase_iter; ++it) {
    bool all_conv = true;
    for (size_t i = 0; i < N; ++i) {
      const Vec3 X = helix_point(ts.par[i], L[i]);
      const Mat35 A = helix_dXdpar(ts.par[i], L[i]);
      const Vec3 a = helix_dXdL(ts.par[i], L[i]);
      const Mat3 Winv = A * ts.cov[i] * A.transpose();
      double cond;
      bool ok;
      const Mat3 W = reg_inv(Winv, cfg.rcond, cond, ok);
      const Vec3 aw = W * a;
      const double denom = a.dot(aw);
      if (!(denom > 0) || !std::isfinite(denom)) continue;  // conv, dL = 0
      const double dL = (x - X).dot(aw) / denom;
      L[i] += dL;
      if (!(std::abs(dL) < cfg.phase_tol)) all_conv = false;
    }
    if (all_conv) break;
  }
  // final evaluation at the refined phases
  out.D.resize(N);
  out.X.resize(N);
  out.chi2.resize(N);
  out.L = L;
  out.cond_max = 0.0;
  for (size_t i = 0; i < N; ++i) {
    const Vec3 X = helix_point(ts.par[i], L[i]);
    const Mat35 A = helix_dXdpar(ts.par[i], L[i]);
    const Vec3 a = helix_dXdL(ts.par[i], L[i]);
    const Mat3 Winv = A * ts.cov[i] * A.transpose();
    double cond;
    bool ok;
    const Mat3 W = reg_inv(Winv, cfg.rcond, cond, ok);
    out.cond_max = std::max(out.cond_max, cond);
    const Vec3 aw = W * a;
    const double denom = a.dot(aw);
    Mat3 Di;
    if (denom > 0 && std::isfinite(denom))
      Di = W - (aw * aw.transpose()) / denom;
    else
      Di = W;
    out.D[i] = 0.5 * (Di + Di.transpose());
    out.X[i] = X;
    const Vec3 r = x - X;
    out.chi2[i] = r.dot(out.D[i] * r);
  }
}

inline double chi2_total(const TrackSet& ts, const Vec3& x,
                         const std::vector<double>& L, const FitConfig& cfg,
                         const BeamSpot* bs, TrackTerms& terms,
                         double& chi2_bs) {
  track_terms(ts, x, L, cfg, terms);
  double tot = 0.0;
  for (double c : terms.chi2) tot += c;
  chi2_bs = 0.0;
  if (bs) {
    const Vec3 d = x - bs->center();
    chi2_bs = d.dot(bs->cov_inv() * d);
    tot += chi2_bs;
  }
  return tot;
}

// fast linear seed: every track projected at its perigee (L = 0), one solve
inline Vec3 seed_linear(const TrackSet& ts, const BeamSpot* bs,
                        const FitConfig& cfg) {
  const size_t N = ts.size();
  Mat3 H = Mat3::Zero();
  Vec3 b = Vec3::Zero();
  for (size_t i = 0; i < N; ++i) {
    const Vec3 X = helix_point(ts.par[i], 0.0);
    const Mat35 A = helix_dXdpar(ts.par[i], 0.0);
    const Vec3 a = helix_dXdL(ts.par[i], 0.0);
    const Mat3 Winv = A * ts.cov[i] * A.transpose();
    double cond;
    bool ok;
    const Mat3 W = reg_inv(Winv, cfg.rcond, cond, ok);
    const Vec3 aw = W * a;
    const double den = a.dot(aw);
    const Mat3 Di =
        (den > 0) ? Mat3(W - (aw * aw.transpose()) / den) : W;
    H += Di;
    b += Di * X;
  }
  if (bs) {
    H += bs->cov_inv();
    b += bs->cov_inv() * bs->center();
  }
  double cond;
  bool ok;
  const Mat3 Hi = reg_inv(H, cfg.rcond, cond, ok);
  if (ok) return Hi * b;
  return bs ? bs->center() : Vec3::Zero();
}

struct RunResult {
  Vec3 x;
  Mat3 cov;
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  double chi2_bs = 0.0;
  std::vector<double> trk_chi2;
  std::vector<double> L;
  int n_iter = 0;
  int status = kMaxIter;
  bool converged = false;
  double condH = std::numeric_limits<double>::infinity();
  double cond_tr = std::numeric_limits<double>::infinity();
  int nrej = 0;
};

// damped Gauss-Newton (Levenberg-Marquardt) from one seed
inline RunResult run_from_seed(const TrackSet& ts, const BeamSpot* bs,
                               const Vec3& x0, const FitConfig& cfg) {
  const size_t N = ts.size();
  RunResult rr;
  Vec3 x = x0;
  std::vector<double> L(N, 0.0);
  TrackTerms tt;
  double chi2_bs;
  double chi2 = chi2_total(ts, x, L, cfg, bs, tt, chi2_bs);
  L = tt.L;
  double lam = cfg.lm_lambda0;
  int nrej = 0;
  int status = kMaxIter;
  int it = 0;
  for (it = 1; it <= cfg.max_iter; ++it) {
    Mat3 H = Mat3::Zero();
    Vec3 g = Vec3::Zero();
    for (size_t i = 0; i < N; ++i) {
      H += tt.D[i];
      g += tt.D[i] * (x - tt.X[i]);
    }
    if (bs) {
      H += bs->cov_inv();
      g += bs->cov_inv() * (x - bs->center());
    }
    bool accepted = false;
    while (true) {
      Mat3 Hd = H;
      for (int k = 0; k < 3; ++k) Hd(k, k) += lam * H(k, k);
      const Vec3 dx = cholesky_solve(Hd, g, cfg.rcond);
      if (!dx.allFinite()) {
        status = kLinalgFail;
        break;
      }
      const Vec3 xn = x + dx;
      if (xn.norm() > cfg.max_radius) {
        // HARD reject: a PV outside the containment ball is not a candidate
        // solution — fail this seed run instead of re-damping.
        status = kDivergedRadius;
        break;
      }
      TrackTerms ttn;
      double chi2_bs_n;
      const double c2n = chi2_total(ts, xn, L, cfg, bs, ttn, chi2_bs_n);
      const double dstep = dx.dot(H * dx);  // undamped H
      const double dchi2 = std::abs(chi2 - c2n);
      const double tolc = cfg.tol_dchi2 + cfg.tol_dchi2_rel * std::abs(chi2);
      if (std::isfinite(c2n) && c2n <= chi2 + 1e-12) {
        x = xn;
        chi2 = c2n;
        chi2_bs = chi2_bs_n;
        tt = ttn;
        L = ttn.L;
        lam = std::max(lam * cfg.lm_down, 1e-14);
        accepted = true;
        if (dstep < cfg.tol_step && dchi2 < tolc) status = kOk;
        break;
      }
      // A rejected step negligible in both the parameter metric and the chi2
      // means only round-off is left: convergence, not a stall.
      if (std::isfinite(c2n) && dstep < cfg.tol_step && dchi2 < tolc) {
        status = kOk;
        break;
      }
      lam *= cfg.lm_up;
      ++nrej;
      if (lam > cfg.lm_lambda_max || nrej > cfg.max_reject) {
        status = kLmStall;
        break;
      }
    }
    if (status == kOk || status == kLmStall || status == kLinalgFail ||
        status == kDivergedRadius)
      break;
    if (!accepted) break;
  }

  Mat3 H = Mat3::Zero();
  for (size_t i = 0; i < N; ++i) H += tt.D[i];
  if (bs) H += bs->cov_inv();
  double condH;
  bool okH;
  const Mat3 covx = reg_inv(H, cfg.rcond, condH, okH);
  bool converged = (status == kOk) && okH && x.allFinite();
  if (status == kOk && !okH) status = kSingularHessian;

  rr.x = x;
  rr.cov = covx;
  rr.chi2 = chi2;
  rr.chi2_bs = chi2_bs;
  rr.trk_chi2 = tt.chi2;
  rr.L = tt.L;
  rr.n_iter = it;
  rr.status = status;
  rr.converged = converged;
  rr.condH = condH;
  rr.cond_tr = tt.cond_max;
  rr.nrej = nrej;
  return rr;
}

inline PVFitResult fit_core(const TrackSet& ts, const BeamSpot* bs,
                            const Vec3* x_start, const FitConfig& cfg) {
  const size_t N = ts.size();
  PVFitResult out;
  out.n_tracks = static_cast<int>(N);
  out.beamspot_used = (bs != nullptr);
  out.track_chi2.resize(N, 0.0);
  out.track_phase.resize(N, 0.0);

  if (N < 2 && bs == nullptr) {
    out.status = kTooFewTracks;
    out.seed_used = kSeedNone;
    return out;
  }

  // deterministic seed ladder
  std::vector<std::pair<int, Vec3>> seeds;
  if (x_start) seeds.emplace_back(kSeedUser, *x_start);
  seeds.emplace_back(kSeedLinear, seed_linear(ts, bs, cfg));
  seeds.emplace_back(kSeedBeamspot, bs ? bs->center() : Vec3::Zero());
  {
    std::vector<double> cx, cy, cz;
    cx.reserve(N);
    cy.reserve(N);
    cz.reserve(N);
    for (size_t i = 0; i < N; ++i) {
      const Vec3 p = helix_point(ts.par[i], 0.0);
      cx.push_back(p(0));
      cy.push_back(p(1));
      cz.push_back(p(2));
    }
    auto median = [](std::vector<double>& v) {
      const size_t n = v.size();
      std::sort(v.begin(), v.end());
      if (n == 0) return std::numeric_limits<double>::quiet_NaN();
      return (n % 2 == 1) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
    };
    seeds.emplace_back(kSeedPerigeeMedian, Vec3(median(cx), median(cy),
                                                median(cz)));
  }

  std::string tried;
  bool have_best = false;
  RunResult best;
  int best_seed = kSeedNone;
  for (const auto& s : seeds) {
    if (!tried.empty()) tried += ",";
    if (!s.second.allFinite()) {
      tried += std::string(seed_name(s.first)) + ":badseed";
      continue;
    }
    RunResult res = run_from_seed(ts, bs, s.second, cfg);
    tried += std::string(seed_name(s.first)) + ":" + status_name(res.status);
    if (!have_best) {
      best = res;
      best_seed = s.first;
      have_best = true;
    }
    if (res.converged) {
      if (!best.converged || res.chi2 < best.chi2 - 1e-9) {
        best = res;
        best_seed = s.first;
      }
      break;  // first converged seed wins (deterministic ladder order)
    }
  }

  out.position = {best.x(0), best.x(1), best.x(2)};
  out.cov = {best.cov(0, 0), best.cov(1, 0), best.cov(1, 1),
             best.cov(2, 0), best.cov(2, 1), best.cov(2, 2)};
  out.chi2 = best.chi2;
  out.chi2_beamspot = best.chi2_bs;
  out.ndf = 2 * static_cast<int>(N) - 3;
  out.ndf_bs = bs ? 2 * static_cast<int>(N) : out.ndf;
  out.n_iterations = best.n_iter;
  out.converged = best.converged;
  out.status = best.status;
  out.track_chi2.assign(best.trk_chi2.begin(), best.trk_chi2.end());
  out.track_phase.assign(best.L.begin(), best.L.end());
  out.cond_hessian = best.condH;
  out.cond_track_max = best.cond_tr;
  out.n_rejected_steps = best.nrej;
  out.seed_used = best_seed;
  out.seeds_tried = tried;
  return out;
}

inline PVSelResult select_core(const TrackSet& ts, const BeamSpot* bs,
                               double chi2_max, const FitConfig& cfg,
                               int min_tracks, bool force_first_removal) {
  const int N = static_cast<int>(ts.size());
  PVSelResult out;
  out.kept.reserve(N);
  std::vector<int> keep(N);
  for (int i = 0; i < N; ++i) keep[i] = i;

  if (N < min_tracks) {
    out.fit = fit_core(ts, bs, nullptr, cfg);
    out.kept.assign(keep.begin(), keep.end());
    out.split_converged = out.fit.converged;
    out.trivial = true;
    out.n_passes = 1;
    return out;
  }

  bool first_pass = true;
  while (true) {
    const TrackSet sub = ts.subset(keep);
    PVFitResult res = fit_core(sub, bs, nullptr, cfg);
    ++out.n_passes;
    if (!res.converged) {
      // Refuse to prune on a non-converged fit: return the split so far,
      // flagged via split_converged = false.
      out.kept.assign(keep.begin(), keep.end());
      out.fit = res;
      out.split_converged = false;
      return out;
    }
    // full-precision double argmax; strict '>' scan => the LOWEST index wins
    // an exact tie; EXACTLY ONE track removed per pass (clean loop)
    double cmax = -std::numeric_limits<double>::infinity();
    int imax = -1;
    for (size_t i = 0; i < res.track_chi2.size(); ++i) {
      if (res.track_chi2[i] > cmax) {
        cmax = res.track_chi2[i];
        imax = static_cast<int>(i);
      }
    }
    // force_first_removal reproduces the reference loop's unconditional first
    // removal; the default tests the threshold first.
    const bool forced = first_pass && force_first_removal;
    first_pass = false;
    if (imax < 0 || (!forced && !(cmax >= chi2_max)) ||
        static_cast<int>(keep.size()) - 1 < min_tracks) {
      out.kept.assign(keep.begin(), keep.end());
      out.fit = res;
      out.split_converged = true;
      return out;
    }
    keep.erase(keep.begin() + imax);
  }
}

}  // namespace detail

// ---------------------------------------------------------------------------
// public interface — edm4hep track states (ALEPH-flipped, cm-native)
// ---------------------------------------------------------------------------

inline PVFitResult fit_vertex(const RVec<edm4hep::TrackState>& tracks,
                              const BeamSpot& bs,
                              const FitConfig& cfg = FitConfig()) {
  const detail::TrackSet ts = detail::convert(tracks);
  return detail::fit_core(ts, &bs, nullptr, cfg);
}

inline PVFitResult fit_vertex_nobs(const RVec<edm4hep::TrackState>& tracks,
                                   const FitConfig& cfg = FitConfig()) {
  const detail::TrackSet ts = detail::convert(tracks);
  return detail::fit_core(ts, nullptr, nullptr, cfg);
}

inline PVFitResult fit_vertex_seeded(const RVec<edm4hep::TrackState>& tracks,
                                     const BeamSpot& bs, double sx, double sy,
                                     double sz,
                                     const FitConfig& cfg = FitConfig()) {
  const detail::TrackSet ts = detail::convert(tracks);
  const Vec3 seed(sx, sy, sz);
  return detail::fit_core(ts, &bs, &seed, cfg);
}

// Iterative pvchi2 pruning with the SAME fitter and the SAME beam-spot
// constraint as the final fit (identity by construction).
inline PVSelResult select_primary_tracks(
    const RVec<edm4hep::TrackState>& tracks, const BeamSpot& bs,
    double chi2_max, const FitConfig& cfg = FitConfig(),
    int min_tracks = 2, bool force_first_removal = false) {
  const detail::TrackSet ts = detail::convert(tracks);
  return detail::select_core(ts, &bs, chi2_max, cfg, min_tracks,
                             force_first_removal);
}

// Per-track chi2 against a FIXED point (no fit): the beamspot-as-fixed-PV
// fallback classification for events whose pruning did not converge.
inline RVec<double> track_chi2_vs_point(const RVec<edm4hep::TrackState>& tracks,
                                        double px, double py, double pz,
                                        const FitConfig& cfg = FitConfig()) {
  const detail::TrackSet ts = detail::convert(tracks);
  const Vec3 x(px, py, pz);
  std::vector<double> L(ts.size(), 0.0);
  detail::TrackTerms tt;
  detail::track_terms(ts, x, L, cfg, tt);
  return RVec<double>(tt.chi2.begin(), tt.chi2.end());
}

// One-shot post-convergence smoothing: updated track parameters (D, phi0, C,
// z0, ct) and momenta [GeV] at the fitted vertex. Refuses on a non-converged
// fit (returns valid = false).
inline SmoothedTracks smooth_at_vertex(const RVec<edm4hep::TrackState>& tracks,
                                       const PVFitResult& res,
                                       double solenoidBz = 1.5,
                                       const FitConfig& cfg = FitConfig()) {
  SmoothedTracks out;
  if (!res.converged ||
      static_cast<int>(tracks.size()) != res.n_tracks)
    return out;  // valid = false
  const detail::TrackSet ts = detail::convert(tracks);
  const Vec3 x(res.position[0], res.position[1], res.position[2]);
  const size_t N = ts.size();
  out.par.reserve(5 * N);
  out.momentum.reserve(N);
  out.phase.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    // Constrained least squares: min (p - par)^T C^-1 (p - par) subject to the
    // helix at p passing through the fitted vertex x, by iterated linearised
    // projection. The vertex itself is never touched.
    Vec5 pnew = ts.par[i];
    double L2 = res.track_phase[i];
    for (int it = 0; it < 3; ++it) {
      // re-solve the phase of the vertex on the current helix
      for (int k = 0; k < 3; ++k) {
        const Vec3 X2 = detail::helix_point(pnew, L2);
        const Vec3 a2 = detail::helix_dXdL(pnew, L2);
        const double dn = a2.dot(a2);
        if (!(dn > 0)) break;
        L2 += (x - X2).dot(a2) / dn;
      }
      const Vec3 X = detail::helix_point(pnew, L2);
      const Mat35 A = detail::helix_dXdpar(pnew, L2);
      const Vec3 a = detail::helix_dXdL(pnew, L2);
      const Mat3 Winv = A * ts.cov[i] * A.transpose();
      double cond;
      bool ok;
      const Mat3 W = detail::reg_inv(Winv, cfg.rcond, cond, ok);
      const Vec3 aw = W * a;
      const double denom = a.dot(aw);
      const Mat3 D = (denom > 0 && std::isfinite(denom))
                         ? Mat3(W - (aw * aw.transpose()) / denom)
                         : W;
      // residual of the constraint linearised about pnew, referred back to
      // the ORIGINAL parameters: p_next = par + C A^T D (x - X + A*(pnew-par))
      const Vec3 rtil = (x - X) + A * (pnew - ts.par[i]);
      pnew = ts.par[i] + ts.cov[i] * A.transpose() * (D * rtil);
    }
    {
      // final phase refresh on the final helix
      for (int k = 0; k < 3; ++k) {
        const Vec3 X2 = detail::helix_point(pnew, L2);
        const Vec3 a2 = detail::helix_dXdL(pnew, L2);
        const double dn = a2.dot(a2);
        if (!(dn > 0)) break;
        L2 += (x - X2).dot(a2) / dn;
      }
    }
    const double ang = pnew(1) + 2.0 * pnew(2) * L2;
    const double absC = std::abs(pnew(2));
    // omega = -2C
    const double pt = absC > 0 ? kPtPerTeslaCm * solenoidBz / (2.0 * absC)
                               : std::numeric_limits<double>::quiet_NaN();
    for (int k = 0; k < 5; ++k) out.par.push_back(pnew(k));
    out.momentum.emplace_back(pt * std::cos(ang), pt * std::sin(ang),
                              pt * pnew(4));
    out.phase.push_back(L2);
  }
  out.valid = true;
  return out;
}

// ---------------------------------------------------------------------------
// stage1 wiring glue: vertex object and primary-track split fallback.
// ---------------------------------------------------------------------------

// PVSelResult -> FCCAnalysesVertex. Position is ALWAYS the fit's answer; the
// covariance is zeroed when the fit did not converge. chi2 is stored as
// chi2/ndf (the production VertexData convention).
inline VertexingUtils::FCCAnalysesVertex toFCCVertex(const PVSelResult& sel) {
  VertexingUtils::FCCAnalysesVertex out;
  edm4hep::VertexData vd;
  vd.position = edm4hep::Vector3f(float(sel.fit.position[0]),
                                  float(sel.fit.position[1]),
                                  float(sel.fit.position[2]));
  std::array<float, 6> cm{};  // zeros
  if (sel.fit.converged)
    for (int k = 0; k < 6; ++k) cm[k] = float(sel.fit.cov[k]);
  vd.covMatrix = cm;
  vd.chi2 = sel.fit.ndf > 0 ? float(sel.fit.chi2 / sel.fit.ndf) : -1.f;
  vd.algorithmType = 2;  // distinguishes the new fitter from the Delphes one
#if EDM4HEP_BUILD_VERSION <= EDM4HEP_VERSION(0, 10, 5)
  vd.primary = 1;
#else
  vd.type = 1;
#endif
  out.vertex = vd;
  out.ntracks = sel.fit.n_tracks;
  out.mc_ind = -1;
  return out;
}

// The primary-track split: the kept set on a fully converged pruning, else
// every track classified against the beam spot as a FIXED point with the same
// chi2 threshold. Fewer than 2 input tracks returns empty.
inline RVec<edm4hep::TrackState> primaryTracksFromSel(
    const RVec<edm4hep::TrackState>& tracks, const PVSelResult& sel,
    double bx, double by, double bz, double chi2_max,
    const FitConfig& cfg = FitConfig()) {
  RVec<edm4hep::TrackState> out;
  if (tracks.size() < 2) return out;
  if (sel.split_converged) {
    for (int i : sel.kept) out.push_back(tracks[i]);
    return out;
  }
  const detail::TrackSet ts = detail::convert(tracks);
  const Vec3 x(bx, by, bz);
  std::vector<double> L(ts.size(), 0.0);
  detail::TrackTerms tt;
  detail::track_terms(ts, x, L, cfg, tt);
  for (size_t i = 0; i < tracks.size(); ++i)
    if (tt.chi2[i] < chi2_max) out.push_back(tracks[i]);
  return out;
}

// ---------------------------------------------------------------------------
// raw-double interface: bypasses the float edm4hep track-state members so the
// offline gates run at full double fidelity.
// ---------------------------------------------------------------------------

inline PVFitResult fit_vertex_raw(const RVec<double>& d0,
                                  const RVec<double>& phi,
                                  const RVec<double>& omega,
                                  const RVec<double>& z0,
                                  const RVec<double>& tanl,
                                  const RVec<double>& cov21_flat,
                                  const BeamSpot& bs, bool use_bs,
                                  double seed_x, double seed_y, double seed_z,
                                  bool use_seed,
                                  const FitConfig& cfg = FitConfig()) {
  detail::TrackSet ts;
  const size_t N = d0.size();
  for (size_t i = 0; i < N; ++i)
    detail::add_track(ts, d0[i], phi[i], omega[i], z0[i], tanl[i],
                      cov21_flat.data() + 21 * i);
  const Vec3 seed(seed_x, seed_y, seed_z);
  return detail::fit_core(ts, use_bs ? &bs : nullptr,
                          use_seed ? &seed : nullptr, cfg);
}

inline PVSelResult select_primary_tracks_raw(
    const RVec<double>& d0, const RVec<double>& phi, const RVec<double>& omega,
    const RVec<double>& z0, const RVec<double>& tanl,
    const RVec<double>& cov21_flat, const BeamSpot& bs,
    double chi2_max,
    const FitConfig& cfg = FitConfig(), int min_tracks = 2,
    bool force_first_removal = false) {
  detail::TrackSet ts;
  const size_t N = d0.size();
  for (size_t i = 0; i < N; ++i)
    detail::add_track(ts, d0[i], phi[i], omega[i], z0[i], tanl[i],
                      cov21_flat.data() + 21 * i);
  return detail::select_core(ts, &bs, chi2_max, cfg, min_tracks,
                             force_first_removal);
}

}  // namespace AlephPVNew
}  // namespace FCCAnalyses

#endif  // ALEPHPVNEW_H
