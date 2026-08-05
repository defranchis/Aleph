"""Standalone primary-vertex (PV) fitter prototype -- numpy, cm-native.

OFFLINE PROTOTYPE.  Not integrated anywhere.

Unit convention (SINGLE convention, end to end, no switches):
    lengths            cm
    curvature omega    1/cm
    angles             rad
    beam-spot sigmas   cm
There is no `Units_mm` flag anywhere in this file, by design: the units defect
that the inherited chain suffered from (selection fit run with a 1000x looser
beam-spot constraint than the final fit) is impossible here because there is
exactly one fitter entry point and it takes physical cm.

Track model (internal parametrisation, "C-parametrisation", identical to the
Delphes TrackCovariance one up to the phase variable):
    a = (D, phi0, C, z0, ct)
    D    signed transverse impact parameter        [cm]
    phi0 azimuth of the momentum at the perigee    [rad]
    C    = -omega/2  (half signed curvature)       [1/cm]
    z0   longitudinal impact parameter             [cm]
    ct   = tanLambda = cot(theta)                  [-]
Trajectory, parametrised by the TRANSVERSE ARC LENGTH L [cm] (this prototype's
phase variable; Delphes uses s = 2*C*L):
    X(L) = ( -D sin(phi0) + f(L) cos(phi0 + C L),
              D cos(phi0) + f(L) sin(phi0 + C L),
              z0 + ct L )                       with f(L) = sin(C L)/C
L is used instead of s because dX/dL is a unit vector in the transverse plane
for every track, independent of the curvature: the Jacobians stay O(1) instead
of carrying 1/(2C) ~ radius-of-curvature factors.  The fitted vertex, its
covariance and the chi2 are invariant under this reparametrisation: any
A -> A + a b^T leaves the phase-projected chi2 unchanged.

chi2 definition (identical to the Delphes VertexFit one):
    per track   chi2_i(x) = min_t (r_i - t a_i)^T W_i (r_i - t a_i)
                          = r_i^T D_i r_i
    with        r_i = x - X_i(L_i),  a_i = dX_i/dL,
                W_i = (A_i C_i A_i^T)^-1,  A_i = dX_i/da  (3x5),
                D_i = W_i - W_i a_i a_i^T W_i / (a_i^T W_i a_i)
    beam spot   chi2_bs(x) = (x - b)^T Cb^-1 (x - b),  Cb = diag(sigma^2)
    total       chi2 = sum_i chi2_i + [chi2_bs if constrained]
i.e. each track contributes 2 degrees of freedom (the along-track direction is
projected out) and the beam spot 3.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field

__all__ = [
    "BeamSpot", "FitConfig", "PVFitResult",
    "edm4hep_to_internal", "flip_aleph_edm4hep",
    "fit_vertex", "select_primary_tracks",
]

# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class BeamSpot:
    """Beam-spot constraint in PHYSICAL units (cm).  Explicit, no conversions.

    sigma_x/y/z are the Gaussian widths of the luminous region; the constraint
    enters the fit as an extra measurement of the vertex position with
    covariance diag(sigma_x^2, sigma_y^2, sigma_z^2).
    """
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    sigma_x: float = 0.02    # 200 um
    sigma_y: float = 0.01    # 100 um
    sigma_z: float = 2.0     # 2 cm

    @property
    def center(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)

    @property
    def cov_inv(self) -> np.ndarray:
        s = np.array([self.sigma_x, self.sigma_y, self.sigma_z], dtype=float)
        if np.any(s <= 0):
            raise ValueError("BeamSpot sigmas must be > 0 (cm)")
        return np.diag(1.0 / s ** 2)


@dataclass(frozen=True)
class FitConfig:
    max_iter: int = 60           # hard cap on damped Gauss-Newton iterations
    tol_step: float = 1e-8       # convergence: dx^T H dx  (dimensionless)
    tol_dchi2: float = 1e-7      # convergence: |chi2_new - chi2_old|, absolute
    tol_dchi2_rel: float = 1e-7  # ... plus this times the current chi2.
    #   The relative part is REQUIRED: evaluating the chi2 involves an inner
    #   Newton loop on the phases and 3x3 eigen-decompositions, so it carries a
    #   noise floor -- ordinarily ~1e-10 * chi2, but a track whose projected
    #   weight matrix is pinned at the rcond eigenvalue floor inflates it to
    #   ~2e-8 * chi2.  A tolerance below that floor can never be met on such
    #   fits, and the fit reports a stall where it has in fact converged.
    lm_lambda0: float = 1e-6     # initial Levenberg-Marquardt damping
    lm_up: float = 10.0          # damping increase on a rejected step
    lm_down: float = 0.1         # damping decrease on an accepted step
    lm_lambda_max: float = 1e12  # give up on the seed above this
    max_reject: int = 30         # cumulative rejected steps before giving up
    phase_iter: int = 6          # inner Newton steps for the per-track phase
    phase_tol: float = 1e-12     # inner phase convergence, cm
    rcond: float = 1e-12         # eigenvalue floor (relative) in reg_inv
    cond_warn: float = 1e12      # condition number above which we flag
    max_radius: float = 100.0    # cm; a seed run leaving this ball is abandoned


@dataclass
class PVFitResult:
    """Everything the fit knows.  ALWAYS fully populated.

    converged == False means the numbers below are a best effort and must not
    be used as a vertex.  There is no silent last-iterate return.
    """
    position: np.ndarray                 # (3,) cm
    covariance: np.ndarray               # (3,3) cm^2
    chi2: float                          # total, incl. beam-spot term
    ndf: int                             # 2*N - 3   (Delphes/production convention)
    ndf_bs: int                          # 2*N       (beam spot counted, if used)
    n_iterations: int
    converged: bool
    status: str                          # 'ok' | 'max_iter' | 'lm_stall' | ...
    track_chi2: np.ndarray               # (N,)
    track_phase: np.ndarray              # (N,) transverse arc length L, cm
    n_tracks: int
    beamspot_used: bool
    cond_hessian: float                  # condition number of the final H
    cond_track_max: float                # worst per-track W-inversion cond.
    n_rejected_steps: int
    seed_used: str
    seeds_tried: tuple = ()
    chi2_beamspot: float = 0.0

    @property
    def chi2_per_ndf(self) -> float:
        return self.chi2 / self.ndf if self.ndf > 0 else np.nan

    def __str__(self) -> str:
        p = self.position
        return ("PVFitResult(pos=(%.6f, %.6f, %.6f) cm  chi2=%.6g ndf=%d "
                "chi2/ndf=%.4g  niter=%d  converged=%s  status=%s  N=%d  "
                "seed=%s  condH=%.3g)" %
                (p[0], p[1], p[2], self.chi2, self.ndf, self.chi2_per_ndf,
                 self.n_iterations, self.converged, self.status,
                 self.n_tracks, self.seed_used, self.cond_hessian))


# --------------------------------------------------------------------------
# input conversion
# --------------------------------------------------------------------------

# packed lower-triangular index of the edm4hep 21-element covMatrix,
# ordered (d0, phi, omega, z0, tanLambda)
_TRI = np.array([[0, 1, 3, 6, 10],
                 [1, 2, 4, 7, 11],
                 [3, 4, 5, 8, 12],
                 [6, 7, 8, 9, 13],
                 [10, 11, 12, 13, 14]])

# (d0, phi, omega, z0, tanLambda) -> (D, phi0, C, z0, ct)
_SCALE = np.array([1.0, 1.0, -0.5, 1.0, 1.0])


def flip_aleph_edm4hep(d0, phi, omega, z0, tanl, cov21):
    """Apply the ALEPH sign convention (d0 -> -d0, omega -> -omega) to RAW
    edm4hep track states and their packed covariance.

    Only needed when driving the prototype straight off the raw EDM4hep files;
    the stage1 inputs the module is meant to consume already carry the flip
    (flipD0_copy upstream).  Kept here so that the offline validation can
    reproduce the production inputs bit-for-bit.
    """
    cov = np.array(cov21, dtype=float, copy=True)
    if cov.ndim == 1:
        cov = cov[None, :]
    for k in (1, 4, 6, 8, 10, 12):     # cross terms of exactly one flipped par
        cov[:, k] = -cov[:, k]
    return (-np.asarray(d0, float), np.asarray(phi, float),
            -np.asarray(omega, float), np.asarray(z0, float),
            np.asarray(tanl, float), cov)


def edm4hep_to_internal(d0, phi, omega, z0, tanl, cov21):
    """(already ALEPH-flipped, cm-native) edm4hep track states -> (par, cov).

    Returns par (N,5) = (D, phi0, C, z0, ct) and cov (N,5,5), cm/rad units.
    """
    d0 = np.atleast_1d(np.asarray(d0, float))
    par = np.stack([d0,
                    np.atleast_1d(np.asarray(phi, float)),
                    np.atleast_1d(np.asarray(omega, float)),
                    np.atleast_1d(np.asarray(z0, float)),
                    np.atleast_1d(np.asarray(tanl, float))], axis=1)
    par = par * _SCALE
    c21 = np.atleast_2d(np.asarray(cov21, float))
    cov = c21[:, _TRI] * _SCALE[None, :, None] * _SCALE[None, None, :]
    return par, cov


# --------------------------------------------------------------------------
# helix model (vectorised over tracks)
# --------------------------------------------------------------------------

def _sinc(u):
    """sin(u)/u, stable at u -> 0 (series to machine precision below 1e-4)."""
    u = np.asarray(u, float)
    small = np.abs(u) < 1e-4
    out = np.empty_like(u)
    us = u[small]
    out[small] = 1.0 - us ** 2 / 6.0 + us ** 4 / 120.0
    ub = u[~small]
    out[~small] = np.sin(ub) / ub
    return out


def _dsinc(u):
    """d/du [sin(u)/u] = (u cos u - sin u)/u^2, stable at u -> 0."""
    u = np.asarray(u, float)
    small = np.abs(u) < 1e-4
    out = np.empty_like(u)
    us = u[small]
    out[small] = -us / 3.0 + us ** 3 / 30.0
    ub = u[~small]
    out[~small] = (ub * np.cos(ub) - np.sin(ub)) / ub ** 2
    return out


def helix_point(par, L):
    """X(L) for each track.  par (N,5), L (N,) -> (N,3) cm."""
    D, p0, C, z0, ct = par.T
    u = C * L
    f = L * _sinc(u)                      # = sin(CL)/C
    ang = p0 + u
    x = -D * np.sin(p0) + f * np.cos(ang)
    y = D * np.cos(p0) + f * np.sin(ang)
    z = z0 + ct * L
    return np.stack([x, y, z], axis=1)


def helix_dXdL(par, L):
    """a = dX/dL, (N,3).  Transverse part is a unit vector by construction."""
    D, p0, C, z0, ct = par.T
    ang = p0 + 2.0 * C * L
    return np.stack([np.cos(ang), np.sin(ang), ct * np.ones_like(ang)], axis=1)


def helix_dXdpar(par, L):
    """A = dX/da at fixed L, (N,3,5)."""
    D, p0, C, z0, ct = par.T
    u = C * L
    s_u = _sinc(u)
    f = L * s_u                                   # sin(CL)/C
    dfdC = L * L * _dsinc(u)                      # d/dC [sin(CL)/C]
    ang = p0 + u
    ca, sa = np.cos(ang), np.sin(ang)
    cp, sp = np.cos(p0), np.sin(p0)
    N = par.shape[0]
    A = np.zeros((N, 3, 5))
    # d/dD
    A[:, 0, 0] = -sp
    A[:, 1, 0] = cp
    # d/dphi0
    A[:, 0, 1] = -D * cp - f * sa
    A[:, 1, 1] = -D * sp + f * ca
    # d/dC   (at fixed L, z does not depend on C -- this is the whole reason
    #         the L-parametrisation is better conditioned than the s one)
    A[:, 0, 2] = dfdC * ca - f * L * sa
    A[:, 1, 2] = dfdC * sa + f * L * ca
    # d/dz0
    A[:, 2, 3] = 1.0
    # d/dct
    A[:, 2, 4] = L
    return A


# --------------------------------------------------------------------------
# deterministic linear algebra
# --------------------------------------------------------------------------

def reg_inv(M, rcond=1e-12):
    """Symmetric inverse via eigen-decomposition with a relative eigenvalue
    floor.  Deterministic, no pivoting, no recursion (unlike Delphes RegInv,
    which is a recursive block inversion whose branch structure depends on
    exact zeros).  Returns (inverse, condition number, ok_flag)."""
    M = 0.5 * (M + M.T)
    w, V = np.linalg.eigh(M)
    wmax = np.max(np.abs(w))
    if wmax <= 0 or not np.isfinite(wmax):
        return np.full_like(M, np.nan), np.inf, False
    floor = rcond * wmax
    ok = np.all(w > floor)
    wc = np.where(w > floor, w, floor)
    cond = wmax / np.min(wc)
    return (V * (1.0 / wc)) @ V.T, cond, ok


def _cholesky_solve(H, g, rcond):
    """Solve H dx = -g.  Cholesky first (H is a positive-definite Hessian by
    construction); on failure fall back to the regularised eigen-inverse."""
    try:
        c = np.linalg.cholesky(H)
        y = np.linalg.solve(c, -g)
        dx = np.linalg.solve(c.T, y)
        return dx, True
    except np.linalg.LinAlgError:
        Hi, _, _ = reg_inv(H, rcond)
        return -Hi @ g, False


# --------------------------------------------------------------------------
# the fit
# --------------------------------------------------------------------------

def _track_terms(par, cov, x, L, cfg):
    """Per-track projected weight D_i, residual r_i, chi2_i, updated phase L_i.

    The phase is refined first (inner Newton on L at fixed x, minimising the
    W-metric distance) and the projection matrices are then evaluated there.
    """
    N = par.shape[0]
    L = L.copy()
    for _ in range(cfg.phase_iter):
        X = helix_point(par, L)
        A = helix_dXdpar(par, L)
        a = helix_dXdL(par, L)
        Winv = A @ cov @ np.transpose(A, (0, 2, 1))
        dL = np.zeros(N)
        conv = np.zeros(N, bool)
        for i in range(N):
            W, _, _ = reg_inv(Winv[i], cfg.rcond)
            aw = W @ a[i]
            denom = a[i] @ aw
            if denom <= 0 or not np.isfinite(denom):
                conv[i] = True
                continue
            dL[i] = (x - X[i]) @ aw / denom
            conv[i] = abs(dL[i]) < cfg.phase_tol
        L = L + dL
        if np.all(conv):
            break

    X = helix_point(par, L)
    A = helix_dXdpar(par, L)
    a = helix_dXdL(par, L)
    Winv = A @ cov @ np.transpose(A, (0, 2, 1))
    Dm = np.zeros((N, 3, 3))
    chi2 = np.zeros(N)
    cond_max = 0.0
    r = x[None, :] - X
    for i in range(N):
        W, cond, _ = reg_inv(Winv[i], cfg.rcond)
        cond_max = max(cond_max, cond)
        aw = W @ a[i]
        denom = a[i] @ aw
        if denom > 0 and np.isfinite(denom):
            Di = W - np.outer(aw, aw) / denom
        else:
            Di = W
        Dm[i] = 0.5 * (Di + Di.T)
        chi2[i] = r[i] @ Dm[i] @ r[i]
    return Dm, X, chi2, L, cond_max


def _chi2_total(par, cov, x, L, cfg, bs):
    Dm, X, c2, Lnew, cond = _track_terms(par, cov, x, L, cfg)
    tot = float(np.sum(c2))
    c2bs = 0.0
    if bs is not None:
        d = x - bs.center
        c2bs = float(d @ bs.cov_inv @ d)
        tot += c2bs
    return tot, Dm, X, c2, Lnew, cond, c2bs


def _seed_linear(par, cov, bs, cfg):
    """Fast linear seed: project every track at its perigee (L=0) and solve the
    resulting linear system once.  Same idea as Delphes VtxFitNoSteer, but with
    no phase pre-stepping and with the regularised inverse."""
    N = par.shape[0]
    L = np.zeros(N)
    X = helix_point(par, L)
    A = helix_dXdpar(par, L)
    a = helix_dXdL(par, L)
    Winv = A @ cov @ np.transpose(A, (0, 2, 1))
    H = np.zeros((3, 3))
    b = np.zeros(3)
    for i in range(N):
        W, _, _ = reg_inv(Winv[i], cfg.rcond)
        aw = W @ a[i]
        den = a[i] @ aw
        Di = W - np.outer(aw, aw) / den if den > 0 else W
        H += Di
        b += Di @ X[i]
    if bs is not None:
        H = H + bs.cov_inv
        b = b + bs.cov_inv @ bs.center
    Hi, _, ok = reg_inv(H, cfg.rcond)
    return Hi @ b if ok else (bs.center if bs is not None else np.zeros(3))


def _run_from_seed(par, cov, bs, x0, cfg):
    """Damped Gauss-Newton (Levenberg-Marquardt) from one seed."""
    N = par.shape[0]
    x = np.asarray(x0, float).copy()
    L = np.zeros(N)
    chi2, Dm, X, c2trk, L, cond_tr, c2bs = _chi2_total(par, cov, x, L, cfg, bs)
    lam = cfg.lm_lambda0
    nrej = 0
    status = "max_iter"
    it = 0
    for it in range(1, cfg.max_iter + 1):
        H = Dm.sum(axis=0)
        g = np.einsum('nij,nj->i', Dm, x[None, :] - X)
        if bs is not None:
            H = H + bs.cov_inv
            g = g + bs.cov_inv @ (x - bs.center)
        accepted = False
        while True:
            Hd = H + lam * np.diag(np.diag(H))
            dx, _ = _cholesky_solve(Hd, g, cfg.rcond)
            if not np.all(np.isfinite(dx)):
                status = "linalg_fail"
                break
            xn = x + dx
            if np.linalg.norm(xn) > cfg.max_radius:
                # a PV outside the beam pipe/detector volume: not a physical
                # solution, treat like a rejected step and re-damp
                lam *= cfg.lm_up
                nrej += 1
                if lam > cfg.lm_lambda_max or nrej > cfg.max_reject:
                    status = "diverged_radius"
                    break
                continue
            c2n, Dmn, Xn, c2trkn, Ln, cond_n, c2bsn = _chi2_total(
                par, cov, xn, L, cfg, bs)
            dstep = float(dx @ H @ dx)
            dchi2 = abs(chi2 - c2n)
            tolc = cfg.tol_dchi2 + cfg.tol_dchi2_rel * abs(chi2)
            if np.isfinite(c2n) and c2n <= chi2 + 1e-12:
                x, chi2, Dm, X, c2trk, L, cond_tr, c2bs = (
                    xn, c2n, Dmn, Xn, c2trkn, Ln, cond_n, c2bsn)
                lam = max(lam * cfg.lm_down, 1e-14)
                accepted = True
                if dstep < cfg.tol_step and dchi2 < tolc:
                    status = "ok"
                break
            # A REJECTED step that is already negligible in both the parameter
            # metric and the chi2 means we are sitting at the minimum and the
            # only thing left is round-off: that is convergence, not a stall.
            if (np.isfinite(c2n) and dstep < cfg.tol_step
                    and dchi2 < tolc):
                status = "ok"
                break
            lam *= cfg.lm_up
            nrej += 1
            if lam > cfg.lm_lambda_max or nrej > cfg.max_reject:
                status = "lm_stall"
                break
        if status in ("ok", "lm_stall", "linalg_fail", "diverged_radius"):
            break
        if not accepted:
            break

    H = Dm.sum(axis=0)
    if bs is not None:
        H = H + bs.cov_inv
    covx, condH, okH = reg_inv(H, cfg.rcond)
    converged = (status == "ok") and okH and np.all(np.isfinite(x))
    if status == "ok" and not okH:
        status = "singular_hessian"
    return dict(x=x, cov=covx, chi2=chi2, chi2_bs=c2bs, track_chi2=c2trk,
                L=L, n_iter=it, status=status, converged=converged,
                condH=condH, cond_tr=cond_tr, nrej=nrej)


def fit_vertex(par, cov, beamspot=None, x_start=None, config=None):
    """Fit one vertex from track parameters.

    par  (N,5)  internal parametrisation (D, phi0, C, z0, ct), cm/rad
    cov  (N,5,5)
    beamspot  BeamSpot or None (constraint on/off is per call)
    x_start   optional explicit seed (cm); the automatic seed ladder is still
              used as fallback if it fails
    Returns PVFitResult, always fully populated, always with `converged`.
    """
    cfg = config or FitConfig()
    par = np.atleast_2d(np.asarray(par, float))
    cov = np.asarray(cov, float).reshape(-1, 5, 5)
    N = par.shape[0]
    bs = beamspot

    if N < 2 and bs is None:
        return PVFitResult(np.full(3, np.nan), np.full((3, 3), np.nan),
                           np.nan, 0, 0, 0, False, "too_few_tracks",
                           np.zeros(N), np.zeros(N), N, False,
                           np.inf, np.inf, 0, "none")

    seeds = []
    if x_start is not None:
        seeds.append(("user", np.asarray(x_start, float)))
    seeds.append(("linear", _seed_linear(par, cov, bs, cfg)))
    seeds.append(("beamspot", bs.center if bs is not None else np.zeros(3)))
    seeds.append(("perigee_median",
                  np.median(helix_point(par, np.zeros(N)), axis=0)))

    tried, best, best_name = [], None, None
    for name, x0 in seeds:
        if not np.all(np.isfinite(x0)):
            tried.append(name + ":badseed")
            continue
        res = _run_from_seed(par, cov, bs, x0, cfg)
        tried.append("%s:%s" % (name, res["status"]))
        # keep the best CONVERGED result; if none converges keep the first
        if best is None:
            best, best_name = res, name
        if res["converged"]:
            if not best["converged"] or res["chi2"] < best["chi2"] - 1e-9:
                best, best_name = res, name
            break   # first converged seed wins (deterministic ladder order)

    ndf = 2 * N - 3
    ndf_bs = 2 * N if bs is not None else ndf
    return PVFitResult(
        position=best["x"], covariance=best["cov"], chi2=best["chi2"],
        ndf=ndf, ndf_bs=ndf_bs, n_iterations=best["n_iter"],
        converged=bool(best["converged"]), status=best["status"],
        track_chi2=best["track_chi2"], track_phase=best["L"], n_tracks=N,
        beamspot_used=bs is not None, cond_hessian=best["condH"],
        cond_track_max=best["cond_tr"], n_rejected_steps=best["nrej"],
        seed_used=best_name, seeds_tried=tuple(tried),
        chi2_beamspot=best["chi2_bs"])


# --------------------------------------------------------------------------
# iterative primary-track selection (the pvchi2-style pruning)
# --------------------------------------------------------------------------

def select_primary_tracks(par, cov, beamspot=None, chi2_max=5.0,
                          config=None, strict_ties=True, min_tracks=2):
    """Iterative pruning, same semantics as get_PrimaryTracks, but using the
    SAME fitter (and therefore the SAME beam-spot constraint) as the final PV
    fit -- that identity is by construction, not by convention.

    Returns (kept_indices, final_PVFitResult, history).
    `history` records the fit result of every pass, so a non-convergence at any
    stage of the pruning is visible, not just at the end.

    strict_ties=True removes every track whose chi2 equals the maximum (the
    production semantics, but without production's index-shift artefact when
    more than one track ties).
    """
    cfg = config or FitConfig()
    par = np.atleast_2d(np.asarray(par, float))
    cov = np.asarray(cov, float).reshape(-1, 5, 5)
    keep = np.arange(par.shape[0])
    history = []
    if keep.size < min_tracks:
        # trivial case (C++ PVSelResult.trivial): fewer than min_tracks tracks
        # enter the pruning; with a beam spot the fit "converges" at/near the
        # beam-spot centre with no track information. Detectable here as
        # len(history) == 1 and kept.size < min_tracks.
        res = fit_vertex(par, cov, beamspot, config=cfg)
        return keep, res, [res]

    while True:
        res = fit_vertex(par[keep], cov[keep], beamspot, config=cfg)
        history.append(res)
        if not res.converged:
            # do NOT prune on the strength of a fit that did not converge
            return keep, res, history
        c2 = res.track_chi2
        cmax = float(np.max(c2))
        if cmax < chi2_max:
            return keep, res, history
        if strict_ties:
            drop = np.where(c2 >= cmax)[0]
        else:
            drop = np.array([int(np.argmax(c2))])
        if keep.size - drop.size < min_tracks:
            return keep, res, history
        keep = np.delete(keep, drop)
