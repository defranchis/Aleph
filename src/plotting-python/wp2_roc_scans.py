#!/usr/bin/env python3
"""
WP2 cut-design ROC scans on the WP1 truth-matched V0 ntuples.

Evaluates the CUT-type design levers of the SV/V0 revisit handoff on the
existing candidate set (structural levers - single-hypothesis assignment,
global claiming, refit - need the new module and are NOT evaluated here):

  L1  displacement upper bound (junk suppression) + lower bound scan
  L5  Armenteros-Podolanski selection (elliptic Ks band, qt/alpha cuts)
  L6  cosPointing re-tune
  L7  (proxy) candidate chi2 cut

For each scan: signal = candidates truth-matched to the hypothesis
(class 1 for Ks-hyp, class 2 for Lambda-hyp), background = everything else
under the same hypothesis. Reports ROC curves (signal efficiency vs background
rejection) and a proposed working point (max Youden J = eff - (1-rej)) -
PROPOSALS ONLY, adoption requires sign-off.

Usage: python3 wp2_roc_scans.py --input '<glob>' --outdir DIR [--label TEXT]
"""

import argparse
import glob
import json
import os

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.EnableImplicitMT()

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COL = {"ks": "#2a78d6", "lam": "#eb6834"}


def rvec_concat(arrs):
    out = [np.asarray(a) for a in arrs if len(a)]
    return np.concatenate(out) if out else np.array([])


def load(files, prefix):
    df = ROOT.RDataFrame("events", files)
    cols = [f"{prefix}_{v}" for v in
            ["class", "pdg", "invM", "dxyz", "alpha", "qt", "cosPointing", "p"]]
    if prefix == "v0n":
        cols.append("v0n_chi2")
    data = df.AsNumpy(cols)
    # strip the prefix so the rest of the script is prefix-agnostic
    return {c.replace(f"{prefix}_", "v0c_", 1): rvec_concat(data[c]) for c in cols}


def roc(ax, sig_vals, bkg_vals, thresholds, keep_below, label, color):
    """ROC for a threshold cut; keep_below=True means candidates pass if val < thr."""
    n_s, n_b = len(sig_vals), len(bkg_vals)
    effs, rejs = [], []
    for t in thresholds:
        if keep_below:
            eff = np.mean(sig_vals < t) if n_s else 0.
            rej = np.mean(bkg_vals >= t) if n_b else 1.
        else:
            eff = np.mean(sig_vals > t) if n_s else 0.
            rej = np.mean(bkg_vals <= t) if n_b else 1.
        effs.append(eff)
        rejs.append(rej)
    effs, rejs = np.array(effs), np.array(rejs)
    j = effs + rejs - 1.
    best = int(np.argmax(j))
    ax.plot(effs, rejs, "-o", ms=3, lw=1.4, color=color, label=label)
    ax.plot(effs[best], rejs[best], "*", ms=14, color=color)
    return {"thr": float(thresholds[best]), "eff": float(effs[best]),
            "rej": float(rejs[best]), "youdenJ": float(j[best])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--outdir", default="wp2_roc")
    ap.add_argument("--label", default="")
    ap.add_argument("--prefix", default="v0c", choices=["v0c", "v0n"],
                    help="candidate branch prefix: v0c = reference finder, v0n = WP2 module")
    args = ap.parse_args()

    files = sorted(glob.glob(args.input))
    if not files:
        raise SystemExit(f"no inputs match {args.input}")
    os.makedirs(args.outdir, exist_ok=True)
    d = load(files, args.prefix)
    lab = f" [{args.label}]" if args.label else ""
    S = {}

    hyps = [("Ks", 310, 1, COL["ks"]), ("Lambda", 3122, 2, COL["lam"])]

    # ---- L1: displacement upper bound (keep if dxyz < thr) + lower bound ----
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
    thr_up = np.array([10., 20., 30., 50., 75., 100., 150., 200., 300., 500., 1000., 5000.])
    thr_lo = np.array([0., 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2., 3., 5., 8.])
    for name, pdg, cls, col in hyps:
        h = d["v0c_pdg"] == pdg
        sig, bkg = d["v0c_dxyz"][h & (d["v0c_class"] == cls)], d["v0c_dxyz"][h & (d["v0c_class"] != cls)]
        S[f"L1_upper_{name}"] = roc(axs[0], sig, bkg, thr_up, True, f"{name} (dxyz<thr)", col)
        S[f"L1_lower_{name}"] = roc(axs[1], sig, bkg, thr_lo, False, f"{name} (dxyz>thr)", col)
    # v0c_dxyz is in cm (skeptic-verified), despite mm labels in FCCAnalyses comments
    for ax, t in zip(axs, ["upper bound scan [cm]", "lower bound scan [cm]"]):
        ax.set_xlabel("signal efficiency"); ax.set_ylabel("background rejection")
        ax.set_title(f"L1 displacement {t}", fontsize=10)
        ax.grid(alpha=0.25, lw=0.5); ax.legend(frameon=False, fontsize=8)
    fig.suptitle(f"L1: displacement bounds{lab}")
    fig.tight_layout(); fig.savefig(f"{args.outdir}/roc_L1_displacement.png", dpi=140); plt.close(fig)

    # ---- L5: Armenteros - elliptic Ks band distance; Lambda alpha windows ----
    # Ks ellipse: alpha in +-0.83, qt peak 0.206 GeV
    ok = d["v0c_alpha"] > -10
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    h = ok & (d["v0c_pdg"] == 310)
    ell = np.sqrt((d["v0c_alpha"] / 0.83) ** 2 + (d["v0c_qt"] / 0.206) ** 2)
    sig, bkg = ell[h & (d["v0c_class"] == 1)], ell[h & (d["v0c_class"] != 1)]
    # band |ell-1| < thr
    band_sig, band_bkg = np.abs(sig - 1.), np.abs(bkg - 1.)
    S["L5_KsEllipseBand"] = roc(ax, band_sig, band_bkg,
                                np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 1.0]),
                                True, "Ks |AP-ellipse dist|<thr", COL["ks"])
    ax.set_xlabel("signal efficiency"); ax.set_ylabel("background rejection")
    ax.grid(alpha=0.25, lw=0.5); ax.legend(frameon=False, fontsize=8)
    fig.suptitle(f"L5: Armenteros elliptic band (Ks){lab}")
    fig.tight_layout(); fig.savefig(f"{args.outdir}/roc_L5_armenteros.png", dpi=140); plt.close(fig)

    # ---- L6: cosPointing ----
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    thr_cp = 1. - np.geomspace(1e-6, 2e-2, 14)
    for name, pdg, cls, col in hyps:
        h = d["v0c_pdg"] == pdg
        sig = d["v0c_cosPointing"][h & (d["v0c_class"] == cls)]
        bkg = d["v0c_cosPointing"][h & (d["v0c_class"] != cls)]
        S[f"L6_cosPointing_{name}"] = roc(ax, sig, bkg, thr_cp[::-1], False, name, col)
    ax.set_xlabel("signal efficiency"); ax.set_ylabel("background rejection")
    ax.grid(alpha=0.25, lw=0.5); ax.legend(frameon=False, fontsize=8)
    fig.suptitle(f"L6: cosPointing scan (current cut 0.999){lab}")
    fig.tight_layout(); fig.savefig(f"{args.outdir}/roc_L6_cospointing.png", dpi=140); plt.close(fig)

    # ---- L7: vertex chi2 (v0n only - well-conditioned single fit) ----
    if "v0c_chi2" in d:
        fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
        thr_c2 = np.array([0.5, 1., 1.5, 2., 3., 4., 5., 7., 9., 10.])
        for name, pdg, cls, col in hyps:
            h = d["v0c_pdg"] == pdg
            sig = d["v0c_chi2"][h & (d["v0c_class"] == cls)]
            bkg = d["v0c_chi2"][h & (d["v0c_class"] != cls)]
            S[f"L7_chi2_{name}"] = roc(axs[0], sig, bkg, thr_c2, True, name, col)
        axs[0].set_title("vertex $\\chi^2$ < thr", fontsize=10)
        # conversion suppression for Lambda: qt > thr (conversions live at qt~0)
        h = d["v0c_pdg"] == 3122
        sig = d["v0c_qt"][h & (d["v0c_class"] == 2)]
        conv = d["v0c_qt"][h & (d["v0c_class"] == 3)]
        S["L5b_qtVeto_LambdaVsConv"] = roc(axs[1], sig, conv,
                                           np.array([0.003, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08]),
                                           False, "$\\Lambda$ vs conversions ($q_T$>thr)", COL["lam"])
        axs[1].set_title("conversion veto ($q_T$), $\\Lambda$ hyp", fontsize=10)
        for ax in axs:
            ax.set_xlabel("signal efficiency"); ax.set_ylabel("background rejection")
            ax.grid(alpha=0.25, lw=0.5); ax.legend(frameon=False, fontsize=8)
        fig.suptitle(f"L7: chi2 + L5b: conversion veto{lab}")
        fig.tight_layout(); fig.savefig(f"{args.outdir}/roc_L7_chi2_conv.png", dpi=140); plt.close(fig)

    with open(f"{args.outdir}/summary.json", "w") as f:
        json.dump(S, f, indent=2)
    with open(f"{args.outdir}/index.html", "w") as f:
        f.write(f"<html><body><h2>WP2 ROC cut scans{lab}</h2>\n"
                "<p>Working points = max Youden J (stars). PROPOSALS ONLY - "
                "structural levers (single-hypothesis, global claiming, single fit) "
                "not evaluated here.</p>\n")
        f.write(f"<pre>{json.dumps(S, indent=2)}</pre>\n")
        for p in ["roc_L1_displacement.png", "roc_L5_armenteros.png", "roc_L6_cospointing.png"]:
            f.write(f'<div><img src="{p}" style="max-width:95%"></div>\n')
        f.write("</body></html>\n")
    print(json.dumps(S, indent=2))


if __name__ == "__main__":
    main()
