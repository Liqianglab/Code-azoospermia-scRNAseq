#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified developmental divergence (mixing entropy) plots for external validation cohorts.

This script reproduces the *continuous* (cell-level) divergence plots:
- GSE235321: Normal vs NOA1/NOA2 (stage-guided continuous pseudotime + sliding-window entropy)
- GSE149512: platform-corrected PCs (continuous pseudotime + sliding-window entropy)

Inputs (CSV, already exported by upstream pipeline):
- GSE235321_externalValidation_cellLevel_pseudotime_entropy_moduleScores.csv
  required columns: Pseudotime2, MixEntropy, Group, PlotStage
- GSE149512_crossPlatformIntegrated_cellLevel_pseudotime_entropy.csv
  required columns: pseudotime_scaled, entropy, Subtype

Outputs:
- ExternalValidation_GSE235321_developmental_divergence_entropy_FINAL_CONTINUOUS.pdf/png
- ExternalValidation_GSE149512_developmental_divergence_entropy_FINAL_CONTINUOUS.pdf/png
- Smoothed curves as CSV (optional): *_smoothedCurve_FINAL_CONTINUOUS.csv

Notes:
- Font: uses 'sans-serif' with fallback list [Arial, Liberation Sans, DejaVu Sans].
  If Arial is not installed, matplotlib will fall back to Liberation/DejaVu; you can swap to Arial in Illustrator.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from statsmodels.nonparametric.smoothers_lowess import lowess
import logging

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]

PALETTE = {
    "Normal": "#1f77b4",
    "NOA1": "#ff7f0e",
    "NOA2": "#2ca02c",
    "Adult": "#1f77b4",
    "AZFa": "#ff7f0e",
    "iNOA": "#2ca02c",
    "KS": "#8c564b",
}

def lowess_smooth(x, y, frac=0.08):
    d = pd.DataFrame({"x": x, "y": y}).sort_values("x")
    sm = lowess(d["y"].values, d["x"].values, frac=frac, return_sorted=True)
    return sm[:, 0], sm[:, 1]

def plot_divergence_continuous(
    df,
    xcol,
    ycol,
    groupcol,
    title,
    out_pdf,
    out_png,
    stage_bands=None,
    arrest_x=None,
    legend_title="Group",
):
    fig = plt.figure(figsize=(12, 6))  # 2:1
    ax = plt.gca()

    # stage background
    if stage_bands is not None:
        y0, y1 = -0.05, 1.05
        for label, x0, x1, color, alpha in stage_bands:
            ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, color=color, lw=0, alpha=alpha))
            ax.text((x0 + x1) / 2, 1.02, label, ha="center", va="top",
                    fontsize=14, fontweight="bold")

    # scatter
    for g in df[groupcol].unique():
        sub = df[df[groupcol] == g]
        ax.scatter(sub[xcol], sub[ycol], s=10, alpha=0.65,
                   color=PALETTE.get(g, None), label=g, edgecolors="none")

    # smooth
    sx, sy = lowess_smooth(df[xcol].values, df[ycol].values, frac=0.08)
    ax.plot(sx, sy, color="black", lw=4)

    # arrest
    if arrest_x is not None:
        ax.axvline(arrest_x, color="red", linestyle="--", lw=2)
        ax.text(arrest_x + 0.1, 0.03, "Arrest / split point",
                color="red", fontsize=13, fontstyle="italic")

    ax.set_title(title, fontsize=26, fontweight="bold", pad=14)
    ax.set_xlabel("Pseudotime (spermatogenesis progression)", fontsize=16)
    ax.set_ylabel("Group-mixing entropy (bits)", fontsize=16)
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle="--", alpha=0.25)
    ax.legend(title=legend_title, frameon=True, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

def compute_stage_ranges_from_plotstage(df, xcol="Pseudotime2", stagecol="PlotStage", order=None):
    if order is None:
        order = list(pd.unique(df[stagecol]))
    ranges = []
    for s in order:
        vals = df.loc[df[stagecol] == s, xcol]
        if len(vals) == 0:
            continue
        ranges.append((s, float(vals.min()), float(vals.max())))
    return ranges

def main():
    # ---------- GSE235321 ----------
    f235 = "GSE235321_externalValidation_cellLevel_pseudotime_entropy_moduleScores.csv"
    df235 = pd.read_csv(f235)

    # stage bands derived from PlotStage
    order235 = ["SPG (Stem)", "Early SPC", "Late SPC", "Round Sperm", "Elongated Sperm"]
    sr = compute_stage_ranges_from_plotstage(df235, xcol="Pseudotime2", stagecol="PlotStage", order=order235)

    BG = {
        "SPG (Stem)": ("#f9d5d3", 0.35),
        "Early SPC": ("#dbeaf7", 0.35),
        "Late SPC": ("#e8d9f5", 0.35),
        "Round Sperm": ("#dff0d8", 0.35),
        "Elongated Sperm": ("#d3eccc", 0.25),
    }
    bands235 = []
    for label, x0, x1 in sr:
        color, alpha = BG.get(label, ("#eeeeee", 0.2))
        # nicer multiline labels to match other plot
        nice = label.replace(" (Stem)", "\n(Stem)").replace(" ", "\n") if label != "SPG (Stem)" else "SPG\n(Stem)"
        nice = "Early\nSPC" if label == "Early SPC" else nice
        nice = "Late\nSPC" if label == "Late SPC" else nice
        nice = "Round\nSperm" if label == "Round Sperm" else nice
        nice = "Elongated\nSperm" if label == "Elongated Sperm" else nice
        bands235.append((nice, x0, x1, color, alpha))

    # arrest point = minimum of smoothed curve within Late SPC band
    sx, sy = lowess_smooth(df235["Pseudotime2"].values, df235["MixEntropy"].values, frac=0.08)
    late = df235[df235["PlotStage"] == "Late SPC"]["Pseudotime2"]
    late0, late1 = float(late.min()), float(late.max())
    m = (sx >= late0) & (sx <= late1)
    arrest235 = float(sx[m][np.argmin(sy[m])])

    out235_pdf = "ExternalValidation_GSE235321_developmental_divergence_entropy_FINAL_CONTINUOUS.pdf"
    out235_png = "ExternalValidation_GSE235321_developmental_divergence_entropy_FINAL_CONTINUOUS.png"
    plot_divergence_continuous(
        df235, "Pseudotime2", "MixEntropy", "Group",
        "Developmental divergence (GSE235321)", out235_pdf, out235_png,
        stage_bands=bands235, arrest_x=arrest235, legend_title="Group"
    )

    # save smoothed curve (optional)
    pd.DataFrame({"pseudotime": sx, "entropy_smooth": sy}).to_csv(
        "GSE235321_developmental_divergence_entropy_smoothedCurve_FINAL_CONTINUOUS.csv", index=False
    )

    # ---------- GSE149512 ----------
    f149 = "GSE149512_crossPlatformIntegrated_cellLevel_pseudotime_entropy.csv"
    df149 = pd.read_csv(f149)

    bands149 = [
        ("SPG\n(Stem)", 0.0, 3.2, "#f9d5d3", 0.35),
        ("Early\nSPC", 3.2, 4.3, "#dbeaf7", 0.35),
        ("Late\nSPC", 4.3, 5.5, "#e8d9f5", 0.35),
        ("Round\nSperm", 5.5, 7.6, "#dff0d8", 0.35),
        ("Elongated\nSperm", 7.6, 10.0, "#d3eccc", 0.25),
    ]

    sx2, sy2 = lowess_smooth(df149["pseudotime_scaled"].values, df149["entropy"].values, frac=0.08)
    # search around the early-to-late meiotic transition
    m2 = (sx2 >= 3.6) & (sx2 <= 4.8)
    arrest149 = float(sx2[m2][np.argmin(sy2[m2])])

    out149_pdf = "ExternalValidation_GSE149512_developmental_divergence_entropy_FINAL_CONTINUOUS.pdf"
    out149_png = "ExternalValidation_GSE149512_developmental_divergence_entropy_FINAL_CONTINUOUS.png"
    plot_divergence_continuous(
        df149, "pseudotime_scaled", "entropy", "Subtype",
        "Developmental Divergence (Platform-corrected PCs)", out149_pdf, out149_png,
        stage_bands=bands149, arrest_x=arrest149, legend_title="Subtype"
    )

    pd.DataFrame({"pseudotime": sx2, "entropy_smooth": sy2}).to_csv(
        "GSE149512_developmental_divergence_entropy_smoothedCurve_FINAL_CONTINUOUS.csv", index=False
    )

    print("Done. Wrote:", out235_pdf, out235_png, out149_pdf, out149_png)

if __name__ == "__main__":
    main()
