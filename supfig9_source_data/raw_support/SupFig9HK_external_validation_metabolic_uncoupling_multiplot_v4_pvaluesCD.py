#!/usr/bin/env python3
"""
External validation: donor-level lactate coupling and metabolic uncoupling (GSE149512)
Adds pairwise p-values (Mann–Whitney U, two-sided) for panels C/D (vs Normal).

Inputs
------
- GSE149512_donorLevel_export_uptake_uncoupling_spermatidFraction.csv
  Required columns:
    sample_code, Class, export_ST_z, uptake_germ_z, uncoupling_z, spermatid_fraction

Outputs
-------
- ExternalValidation_GSE149512_metabolic_uncoupling_MULTIPANEL_v4_pvaluesCD.pdf
- ExternalValidation_GSE149512_metabolic_uncoupling_MULTIPANEL_v4_pvaluesCD.png
- GSE149512_metabolic_uncoupling_panelsCD_pvalues.csv
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu

# Embed fonts as editable text in PDF
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

# Prefer Arial if present; otherwise matplotlib will fallback
mpl.rcParams["font.family"] = ["Arial", "DejaVu Sans"]

ORDER = ["Normal", "AZFa", "iNOA", "KS"]
COLORS = {
    "Normal": "#1b9e77",  # green
    "AZFa": "#d95f02",    # orange
    "iNOA": "#7570b3",    # purple
    "KS": "#a6761d",      # brown
}

def add_bracket(ax: plt.Axes, x1: float, x2: float, y: float, h: float, text: str, fontsize: int = 10) -> None:
    """Draw a significance bracket and a text label."""
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", lw=1)
    ax.text((x1 + x2) / 2, y + h, text, ha="center", va="bottom", fontsize=fontsize)

def mann_whitney_p(a: np.ndarray, b: np.ndarray) -> float:
    """Two-sided Mann–Whitney U p-value."""
    if len(a) == 0 or len(b) == 0:
        return np.nan
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)

def main() -> None:
    in_csv = "GSE149512_donorLevel_export_uptake_uncoupling_spermatidFraction.csv"
    if not os.path.exists(in_csv):
        raise FileNotFoundError(f"Missing input: {in_csv}")

    df = pd.read_csv(in_csv)

    # Correlations for panels A/B
    rhoA, pA = spearmanr(df["export_ST_z"], df["uptake_germ_z"])
    rhoB, pB = spearmanr(df["uncoupling_z"], df["spermatid_fraction"])

    # Pairwise p-values vs Normal for panels C/D
    rows = []
    normal_export = df.loc[df["Class"] == "Normal", "export_ST_z"].dropna().to_numpy()
    normal_uptake = df.loc[df["Class"] == "Normal", "uptake_germ_z"].dropna().to_numpy()

    p_export = {}
    p_uptake = {}
    for g in ORDER[1:]:
        g_export = df.loc[df["Class"] == g, "export_ST_z"].dropna().to_numpy()
        g_uptake = df.loc[df["Class"] == g, "uptake_germ_z"].dropna().to_numpy()

        p_export[g] = mann_whitney_p(normal_export, g_export)
        p_uptake[g] = mann_whitney_p(normal_uptake, g_uptake)

        rows.append({
            "metric": "ST_lactate_export_z",
            "comparison": f"{g} vs Normal",
            "p_value": p_export[g],
            "n_normal": len(normal_export),
            "n_group": len(g_export),
        })
        rows.append({
            "metric": "Germ_lactate_uptake_oxidation_z",
            "comparison": f"{g} vs Normal",
            "p_value": p_uptake[g],
            "n_normal": len(normal_uptake),
            "n_group": len(g_uptake),
        })

    pd.DataFrame(rows).to_csv("GSE149512_metabolic_uncoupling_panelsCD_pvalues.csv", index=False)

    # ---- Plot (2x2) ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axA, axB, axC, axD = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    # Panel A: export vs uptake (Spearman)
    for g in ORDER:
        sub = df[df["Class"] == g]
        axA.scatter(
            sub["export_ST_z"], sub["uptake_germ_z"],
            s=90, color=COLORS[g], edgecolor="black", linewidth=0.6, label=g
        )
        for _, r in sub.iterrows():
            axA.text(r["export_ST_z"] + 0.03, r["uptake_germ_z"] + 0.03, str(r["sample_code"]), fontsize=9)

    lims = [
        min(df["export_ST_z"].min(), df["uptake_germ_z"].min()) - 0.2,
        max(df["export_ST_z"].max(), df["uptake_germ_z"].max()) + 0.2,
    ]
    axA.plot(lims, lims, linestyle="--", color="black", lw=1)
    axA.set_xlim(lims)
    axA.set_ylim(lims)
    axA.set_title("Coupling between lactate export and uptake", fontsize=12, fontweight="bold")
    axA.set_xlabel("ST lactate export (z)")
    axA.set_ylabel("Germ lactate uptake/oxidation (z)")
    axA.text(
        0.02, 0.96,
        f"Spearman ρ={rhoA:.2f}, P={pA:.3f}",
        transform=axA.transAxes, ha="left", va="top", fontsize=10
    )
    axA.grid(True, linestyle="--", alpha=0.25)

    # Panel B: uncoupling vs spermatid fraction (Spearman + regression)
    for g in ORDER:
        sub = df[df["Class"] == g]
        axB.scatter(
            sub["uncoupling_z"], sub["spermatid_fraction"],
            s=90, color=COLORS[g], edgecolor="black", linewidth=0.6, label=g
        )
        for _, r in sub.iterrows():
            axB.text(r["uncoupling_z"] + 0.03, r["spermatid_fraction"] + 0.002, str(r["sample_code"]), fontsize=9)

    x = df["uncoupling_z"].to_numpy()
    y = df["spermatid_fraction"].to_numpy()
    coef = np.polyfit(x, y, 1)
    xline = np.linspace(x.min(), x.max(), 100)
    axB.plot(xline, coef[0] * xline + coef[1], color="black", lw=2)

    axB.set_title("Uncoupling predicts post-meiotic output", fontsize=12, fontweight="bold")
    axB.set_xlabel("Uncoupling index (export_z − uptake_z)")
    axB.set_ylabel("Spermatid fraction (post-meiotic / germ)")
    axB.text(
        0.02, 0.96,
        f"Spearman ρ={rhoB:.2f}, P={pB:.4f}",
        transform=axB.transAxes, ha="left", va="top", fontsize=10
    )
    axB.grid(True, linestyle="--", alpha=0.25)

    # Legend
    handles, labels = axB.get_legend_handles_labels()
    leg = axB.legend(handles, labels, title="Group", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True)
    for h in leg.legend_handles:
        try:
            h.set_edgecolor("black")
            h.set_linewidth(0.6)
        except Exception:
            pass

    # Panel C: export by group + p-values (vs Normal)
    metric = "export_ST_z"
    data = [df[df["Class"] == g][metric].dropna().to_numpy() for g in ORDER]
    axC.boxplot(
        data, positions=np.arange(1, len(ORDER) + 1), widths=0.55, patch_artist=True,
        medianprops=dict(color="black", linewidth=2),
        boxprops=dict(facecolor="white", color="black"),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
    )
    for i, g in enumerate(ORDER, start=1):
        vals = df[df["Class"] == g][metric].dropna().to_numpy()
        axC.scatter(
            np.full_like(vals, i, dtype=float) + np.random.uniform(-0.05, 0.05, size=len(vals)),
            vals, s=90, color=COLORS[g], edgecolor="black", linewidth=0.6, zorder=3
        )
    axC.set_xticks(range(1, len(ORDER) + 1))
    axC.set_xticklabels(ORDER, rotation=20, ha="right")
    axC.set_title("ST lactate export by group", fontsize=12, fontweight="bold")
    axC.set_ylabel("ST export score (z)")
    axC.grid(True, axis="y", linestyle="--", alpha=0.3)

    ymin, ymax = axC.get_ylim()
    yr = ymax - ymin
    base = ymax + 0.05 * yr
    step = 0.07 * yr
    h = 0.02 * yr
    for j, g in enumerate(ORDER[1:], start=2):
        add_bracket(axC, 1, j, base + (j - 2) * step, h, f"P={p_export[g]:.3f}")
    axC.set_ylim(ymin, base + (len(ORDER) - 2) * step + 0.12 * yr)

    # Panel D: uptake by group + p-values (vs Normal)
    metric = "uptake_germ_z"
    data = [df[df["Class"] == g][metric].dropna().to_numpy() for g in ORDER]
    axD.boxplot(
        data, positions=np.arange(1, len(ORDER) + 1), widths=0.55, patch_artist=True,
        medianprops=dict(color="black", linewidth=2),
        boxprops=dict(facecolor="white", color="black"),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
    )
    for i, g in enumerate(ORDER, start=1):
        vals = df[df["Class"] == g][metric].dropna().to_numpy()
        axD.scatter(
            np.full_like(vals, i, dtype=float) + np.random.uniform(-0.05, 0.05, size=len(vals)),
            vals, s=90, color=COLORS[g], edgecolor="black", linewidth=0.6, zorder=3
        )
    axD.set_xticks(range(1, len(ORDER) + 1))
    axD.set_xticklabels(ORDER, rotation=20, ha="right")
    axD.set_title("Germ lactate uptake/oxidation by group", fontsize=12, fontweight="bold")
    axD.set_ylabel("Germ uptake/oxidation score (z)")
    axD.grid(True, axis="y", linestyle="--", alpha=0.3)

    ymin, ymax = axD.get_ylim()
    yr = ymax - ymin
    base = ymax + 0.05 * yr
    step = 0.07 * yr
    h = 0.02 * yr
    for j, g in enumerate(ORDER[1:], start=2):
        add_bracket(axD, 1, j, base + (j - 2) * step, h, f"P={p_uptake[g]:.3f}")
    axD.set_ylim(ymin, base + (len(ORDER) - 2) * step + 0.12 * yr)

    fig.suptitle("GSE149512: donor-level lactate coupling and metabolic uncoupling", fontsize=16, fontweight="bold", y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    out_pdf = "ExternalValidation_GSE149512_metabolic_uncoupling_MULTIPANEL_v4_pvaluesCD.pdf"
    out_png = "ExternalValidation_GSE149512_metabolic_uncoupling_MULTIPANEL_v4_pvaluesCD.png"
    plt.savefig(out_pdf)
    plt.savefig(out_png, dpi=300)
    plt.close(fig)

    print("Saved:", out_pdf)
    print("Saved:", out_png)
    print("Saved:", "GSE149512_metabolic_uncoupling_panelsCD_pvalues.csv")

if __name__ == "__main__":
    main()
