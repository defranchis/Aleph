#!/usr/bin/env python3
"""
WP1 baseline figures of merit for the CURRENT V0 finder (SV/V0 revisit).

Reads stage1 --truthV0 output ntuples (tester or full production) and produces:
  1. efficiency vs |p|, flight distance, cos(theta)  (Ks and Lambda separately)
  2. purity vs invariant mass                        (per hypothesis)
  3. decomposed mass spectra stacked by truth class  (Ks window, Lambda region, full log)
  4. duplicate accounting                            (multi-hypothesis pairs, shared tracks)
  5. kinematic sanity: dxyz by truth class, Armenteros-Podolanski plane

Truth classes (v0c_class): 0 combinatorial, 1 true Ks, 2 true Lambda,
3 gamma-conversion, 4 half-match (one true daughter + stolen/junk track).

Usage:
  python3 wp1_v0_baseline.py --input '<glob>' [--outdir DIR] [--label TEXT]

Writes PNGs + summary.json + index.html into --outdir.
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

# validated categorical palette (dataviz reference, light mode) + neutral for background
COL = {
    "ks":     "#2a78d6",  # true Ks
    "lam":    "#eb6834",  # true Lambda
    "conv":   "#1baf7a",  # gamma conversion
    "half":   "#eda100",  # half-match
    "comb":   "#9a9a94",  # combinatorial (neutral)
}
CLASS_ORDER = [("comb", 0), ("half", 4), ("conv", 3), ("lam", 2), ("ks", 1)]
CLASS_LABEL = {
    "ks": "true $K_S$", "lam": "true $\\Lambda$", "conv": "$\\gamma$ conversion",
    "half": "half-match", "comb": "combinatorial",
}

TREE = "events"


def rvec_concat(arrs):
    """Concatenate a column of per-event RVecs into one flat numpy array."""
    out = [np.asarray(a) for a in arrs if len(a)]
    return np.concatenate(out) if out else np.array([])


def load(files):
    df = ROOT.RDataFrame(TREE, files)
    cols = [
        "truev0_pdg", "truev0_p", "truev0_costheta", "truev0_fd", "truev0_dpv",
        "truev0_nmatched", "truev0_found_any", "truev0_found_correct",
        "v0c_class", "v0c_pairmult", "v0c_trackshared",
        "v0c_alpha", "v0c_qt", "v0c_pdg", "v0c_invM", "v0c_dxyz", "v0c_p",
        "n_v0_event",
    ]
    data = df.AsNumpy(cols)
    flat = {c: rvec_concat(data[c]) for c in cols if c != "n_v0_event"}
    flat["n_events"] = len(data["n_v0_event"])
    flat["n_v0_event_sum"] = int(np.sum(data["n_v0_event"]))
    return flat


def binomial_eff(num, den):
    num = num.astype(float); den = den.astype(float)
    eff = np.divide(num, den, out=np.zeros_like(num, dtype=float), where=den > 0)
    err = np.sqrt(np.maximum(eff * (1 - eff), 1e-12) / np.maximum(den, 1))
    err[den == 0] = 0.
    return eff, err


def eff_plot(ax, x_true, sel_den, sel_num_any, sel_num_cor, bins, xlabel, logx=False):
    den, _ = np.histogram(x_true[sel_den], bins=bins)
    num_a, _ = np.histogram(x_true[sel_den & sel_num_any], bins=bins)
    num_c, _ = np.histogram(x_true[sel_den & sel_num_cor], bins=bins)
    ctr = 0.5 * (bins[1:] + bins[:-1])
    for num, color, lab in [(num_a, "#2a78d6", "found (any hypothesis)"),
                            (num_c, "#eb6834", "found (correct hypothesis)")]:
        eff, err = binomial_eff(num, den)
        ax.errorbar(ctr, eff, yerr=err, fmt="o", ms=4, lw=1.4, color=color, label=lab)
    ax.set_xlabel(xlabel); ax.set_ylabel("efficiency")
    ax.set_ylim(0, 1.05)
    if logx:
        ax.set_xscale("log")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25, lw=0.5)


def stacked_mass(ax, d, hyp_pdg, bins, logy=False):
    sel_h = d["v0c_pdg"] == hyp_pdg
    bottoms = np.zeros(len(bins) - 1)
    ctr = 0.5 * (bins[1:] + bins[:-1])
    width = np.diff(bins)
    for key, cls in CLASS_ORDER:
        h, _ = np.histogram(d["v0c_invM"][sel_h & (d["v0c_class"] == cls)], bins=bins)
        ax.bar(ctr, h, width=width, bottom=bottoms, color=COL[key],
               label=CLASS_LABEL[key], edgecolor="white", linewidth=0.3)
        bottoms += h
    # overlay: multi-hypothesis duplicates
    hdup, _ = np.histogram(d["v0c_invM"][sel_h & (d["v0c_pairmult"] > 1)], bins=bins)
    ax.step(bins, np.append(hdup, hdup[-1]), where="post", color="black",
            lw=1.2, ls="--", label="pair booked $>1\\times$")
    ax.set_xlabel("$m_{inv}$ [GeV]"); ax.set_ylabel("candidates / bin")
    if logy:
        ax.set_yscale("log"); ax.set_ylim(bottom=0.5)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.25, lw=0.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="input ROOT file glob")
    ap.add_argument("--outdir", default="wp1_plots")
    ap.add_argument("--label", default="", help="extra label (e.g. flavour) for titles")
    args = ap.parse_args()

    files = sorted(glob.glob(args.input))
    if not files:
        raise SystemExit(f"no inputs match {args.input}")
    os.makedirs(args.outdir, exist_ok=True)
    d = load(files)
    lab = f" [{args.label}]" if args.label else ""

    S = {}  # summary numbers
    S["n_events"] = d["n_events"]
    S["n_candidates"] = int(len(d["v0c_class"]))
    assert S["n_candidates"] == d["n_v0_event_sum"], "candidate count mismatch vs n_v0_event"

    # ---------------- efficiency + acceptance ----------------
    for name, sel_pdg in [("Ks", np.abs(d["truev0_pdg"]) == 310),
                          ("Lambda", np.abs(d["truev0_pdg"]) == 3122)]:
        n_all = int(np.sum(sel_pdg))
        acc = sel_pdg & (d["truev0_nmatched"] == 2)
        n_acc = int(np.sum(acc))
        n_any = int(np.sum(acc & (d["truev0_found_any"] == 1)))
        n_cor = int(np.sum(acc & (d["truev0_found_correct"] == 1)))
        S[name] = {
            "n_true": n_all, "n_bothTrackMatched": n_acc,
            "acceptance": n_acc / n_all if n_all else 0.,
            "eff_any": n_any / n_acc if n_acc else 0.,
            "eff_correct": n_cor / n_acc if n_acc else 0.,
        }
        fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
        sel_any = d["truev0_found_any"] == 1
        sel_cor = d["truev0_found_correct"] == 1
        eff_plot(axs[0], d["truev0_p"], acc, sel_any, sel_cor,
                 np.linspace(0, 40, 17), "true $|p|$ [GeV]")
        eff_plot(axs[1], d["truev0_dpv"], acc, sel_any, sel_cor,
                 np.geomspace(0.05, 150, 15), "true decay distance from PV [cm]", logx=True)
        eff_plot(axs[2], d["truev0_costheta"], acc, sel_any, sel_cor,
                 np.linspace(-1, 1, 17), "true $\\cos\\theta$")
        fig.suptitle(f"{name} reconstruction efficiency (both daughters track-matched){lab}")
        fig.tight_layout()
        fig.savefig(f"{args.outdir}/eff_{name}.png", dpi=140)
        plt.close(fig)

    # ---------------- purity vs mass ----------------
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, (hyp, cls_sig, rng, name) in zip(
            axs, [(310, 1, (0.3, 0.75), "$K_S$ hypothesis"),
                  (3122, 2, (1.05, 1.35), "$\\Lambda$ hypothesis")]):
        bins = np.linspace(*rng, 26)
        sel_h = d["v0c_pdg"] == hyp
        den, _ = np.histogram(d["v0c_invM"][sel_h], bins=bins)
        num, _ = np.histogram(d["v0c_invM"][sel_h & (d["v0c_class"] == cls_sig)], bins=bins)
        eff, err = binomial_eff(num, den)
        ctr = 0.5 * (bins[1:] + bins[:-1])
        ax.errorbar(ctr, eff, yerr=err, fmt="o", ms=4, color="#2a78d6")
        ax.set_xlabel("$m_{inv}$ [GeV]"); ax.set_ylabel("truth-matched fraction")
        ax.set_ylim(0, 1.05); ax.set_title(name, fontsize=10)
        ax.grid(alpha=0.25, lw=0.5)
    fig.suptitle(f"candidate purity vs mass{lab}")
    fig.tight_layout()
    fig.savefig(f"{args.outdir}/purity_vs_mass.png", dpi=140)
    plt.close(fig)

    # per-hypothesis integrated purity
    for hyp, cls_sig, name in [(310, 1, "Ks"), (3122, 2, "Lambda")]:
        sel_h = d["v0c_pdg"] == hyp
        n_h = int(np.sum(sel_h))
        S[f"purity_{name}"] = {
            "n_candidates": n_h,
            "frac_signal": float(np.mean(d["v0c_class"][sel_h] == cls_sig)) if n_h else 0.,
            "frac_comb": float(np.mean(d["v0c_class"][sel_h] == 0)) if n_h else 0.,
            "frac_conv": float(np.mean(d["v0c_class"][sel_h] == 3)) if n_h else 0.,
            "frac_half": float(np.mean(d["v0c_class"][sel_h] == 4)) if n_h else 0.,
        }

    # ---------------- decomposed mass spectra ----------------
    fig, axs = plt.subplots(1, 3, figsize=(14, 4))
    stacked_mass(axs[0], d, 310, np.linspace(0.35, 0.65, 61))
    axs[0].set_title("$K_S$ hypothesis, signal window", fontsize=10)
    stacked_mass(axs[1], d, 3122, np.linspace(1.05, 1.35, 61))
    axs[1].set_title("$\\Lambda$ hypothesis, signal window", fontsize=10)
    stacked_mass(axs[2], d, 310, np.linspace(0.1, 1.45, 91), logy=True)
    axs[2].set_title("$K_S$ hypothesis, full range", fontsize=10)
    fig.suptitle(f"reco mass spectra decomposed by truth class{lab}")
    fig.tight_layout()
    fig.savefig(f"{args.outdir}/mass_decomposed.png", dpi=140)
    plt.close(fig)

    # ---------------- duplicate accounting ----------------
    n_c = len(d["v0c_class"])
    S["duplicates"] = {
        "frac_pair_multibooked": float(np.mean(d["v0c_pairmult"] > 1)) if n_c else 0.,
        "frac_track_shared": float(np.mean(d["v0c_trackshared"] == 1)) if n_c else 0.,
    }

    # ---------------- kinematic sanity ----------------
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    bins = np.geomspace(0.05, 2e4, 60)
    for key, cls in CLASS_ORDER:
        sel = d["v0c_class"] == cls
        axs[0].hist(np.clip(d["v0c_dxyz"][sel], bins[0], bins[-1]), bins=bins,
                    histtype="step", lw=1.5, color=COL[key], label=CLASS_LABEL[key])
    axs[0].set_xscale("log"); axs[0].set_yscale("log")
    # branch values are cm (skeptic-verified 2026-07-22), despite mm labels in FCCAnalyses comments
    axs[0].set_xlabel("candidate $d_{xyz}$ [cm]"); axs[0].set_ylabel("candidates")
    axs[0].legend(frameon=False, fontsize=8); axs[0].grid(alpha=0.25, lw=0.5)

    ok = d["v0c_alpha"] > -10
    for key, cls in CLASS_ORDER:
        sel = ok & (d["v0c_class"] == cls)
        axs[1].scatter(d["v0c_alpha"][sel], d["v0c_qt"][sel], s=3, alpha=0.5,
                       color=COL[key], label=CLASS_LABEL[key], rasterized=True)
    axs[1].set_xlabel("$\\alpha = (p_L^+ - p_L^-)/(p_L^+ + p_L^-)$")
    axs[1].set_ylabel("$q_T$ [GeV]")
    axs[1].set_xlim(-1.2, 1.2); axs[1].set_ylim(0, 0.35)
    axs[1].legend(frameon=False, fontsize=8, markerscale=3)
    axs[1].grid(alpha=0.25, lw=0.5)
    fig.suptitle(f"kinematic sanity: displacement and Armenteros-Podolanski{lab}")
    fig.tight_layout()
    fig.savefig(f"{args.outdir}/kinematics.png", dpi=140)
    plt.close(fig)

    # ---------------- summary ----------------
    with open(f"{args.outdir}/summary.json", "w") as f:
        json.dump(S, f, indent=2)

    pngs = ["eff_Ks.png", "eff_Lambda.png", "purity_vs_mass.png",
            "mass_decomposed.png", "kinematics.png"]
    with open(f"{args.outdir}/index.html", "w") as f:
        f.write(f"<html><body><h2>WP1 V0 baseline{lab}</h2>\n")
        f.write(f"<pre>{json.dumps(S, indent=2)}</pre>\n")
        for p in pngs:
            f.write(f'<div><img src="{p}" style="max-width:95%"></div>\n')
        f.write("</body></html>\n")

    print(json.dumps(S, indent=2))
    print(f"wrote {args.outdir}/")


if __name__ == "__main__":
    main()
