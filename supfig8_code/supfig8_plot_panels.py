#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

BASE = Path(__file__).resolve().parents[1]  # .../source data
SRC_DIR = BASE / "supfig8_source_data"
PLOT_DIR = BASE / "supfig8_code" / "supfig8_plots"

GROUP_COLORS = {
    "Ctrl": "#1f77b4",
    "OA": "#d62728",
    "AZFc_Del": "#9467bd",
    "iNOA_B": "#ff7f0e",
    "iNOA_S": "#2ca02c",
    "KS": "#8c564b",
}


def ensure_dir() -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)


def save(fig: plt.Figure, name: str) -> None:
    fig.savefig(PLOT_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(PLOT_DIR / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def plot_panel_a() -> None:
    mat = pd.read_csv(SRC_DIR / "SupFig8A_somatic_scores_heatmap_z_matrix.csv")
    grp = pd.read_csv(SRC_DIR / "SupFig8A_group_bar_annotation.csv")

    row_names = mat["score_display"].tolist()
    col_names = [c for c in mat.columns if c != "score_display"]
    z = mat[col_names].to_numpy(dtype=float)

    color_list = [GROUP_COLORS.get(g, "#999999") for g in grp["Group"]]
    group_idx = np.arange(len(color_list))[None, :]

    fig = plt.figure(figsize=(12, 4.8))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[0.25, 1.0], hspace=0.08)

    ax0 = fig.add_subplot(gs[0, 0])
    ax0.imshow(group_idx, aspect="auto", cmap=ListedColormap(color_list))
    ax0.set_yticks([])
    ax0.set_xticks([])
    ax0.set_title("Somatic niche scores (group bar + row-wise z-score)", fontsize=12, pad=4)

    ax1 = fig.add_subplot(gs[1, 0])
    im = ax1.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax1.set_yticks(np.arange(len(row_names)))
    ax1.set_yticklabels(row_names, fontsize=9)
    ax1.set_xticks(np.arange(len(col_names)))
    ax1.set_xticklabels(col_names, rotation=45, ha="right", fontsize=8)
    ax1.set_xlabel("Donor / sample")

    cbar = fig.colorbar(im, ax=ax1, fraction=0.02, pad=0.01)
    cbar.set_label("Z-score", fontsize=9)

    handles = [
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=c, markeredgecolor="none", markersize=8, label=g)
        for g, c in GROUP_COLORS.items()
    ]
    ax0.legend(handles=handles, loc="upper right", ncol=3, frameon=False, fontsize=8)

    save(fig, "SupFig8A_heatmap")


def plot_panel_b() -> None:
    df = pd.read_csv(SRC_DIR / "SupFig8B_pca_points.csv")

    fig, ax = plt.subplots(figsize=(5.8, 4.6))
    for g, sub in df.groupby("Group", observed=False):
        ax.scatter(sub["PC1_v2"], sub["PC2_v2"], s=42, color=GROUP_COLORS.get(g, "#999999"), label=g, edgecolor="none")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Panel B: PCA (PC1_v2 vs PC2_v2)", fontsize=12)
    save(fig, "SupFig8B_pca")


def _scatter_with_fit(path_csv: Path, path_summary: Path, panel: str, ylabel: str, out_name: str) -> None:
    df = pd.read_csv(path_csv)
    sm = pd.read_csv(path_summary).iloc[0]

    fig, ax = plt.subplots(figsize=(5.8, 4.6))
    for g, sub in df.groupby("Group", observed=False):
        ax.scatter(sub.iloc[:, 2], sub.iloc[:, 3], s=42, color=GROUP_COLORS.get(g, "#999999"), label=g, edgecolor="none")

    x = df.iloc[:, 2].to_numpy(dtype=float)
    slope = float(sm["line_slope"])
    intercept = float(sm["line_intercept"])
    xx = np.linspace(np.nanmin(x), np.nanmax(x), 100)
    yy = slope * xx + intercept
    ax.plot(xx, yy, color="black", linewidth=1)

    rho = float(sm["spearman_rho"])
    pval = float(sm["p_value"])
    ax.text(
        0.98,
        0.06,
        f"Spearman rho = {rho:.2f}\nP = {pval:.6g}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
    )

    ax.set_xlabel("Niche collapse axis (PC1v2)")
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(f"Panel {panel}", fontsize=12)
    save(fig, out_name)


def plot_panel_c() -> None:
    _scatter_with_fit(
        SRC_DIR / "SupFig8C_niche_vs_germ_maturity.csv",
        SRC_DIR / "SupFig8C_spearman_summary.csv",
        panel="C",
        ylabel="Germ maturity index",
        out_name="SupFig8C_scatter",
    )


def plot_panel_d() -> None:
    _scatter_with_fit(
        SRC_DIR / "SupFig8D_niche_vs_germ_latest_stage.csv",
        SRC_DIR / "SupFig8D_spearman_summary.csv",
        panel="D",
        ylabel="Germline latest stage (0-7)",
        out_name="SupFig8D_scatter",
    )


def plot_panel_e() -> None:
    mat = pd.read_csv(SRC_DIR / "SupFig8E_score_correlation_spearman_matrix.csv", index_col=0)
    txt = pd.read_csv(SRC_DIR / "SupFig8E_score_correlation_text_matrix.csv", index_col=0)

    val = mat.to_numpy(dtype=float)
    labels = mat.index.tolist()

    fig, ax = plt.subplots(figsize=(8.2, 7.2))
    im = ax.imshow(val, cmap="RdBu_r", vmin=-1, vmax=1)

    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)

    for i in range(len(labels)):
        for j in range(len(labels)):
            s = txt.iloc[i, j]
            if isinstance(s, str) and s:
                ax.text(j, i, s, ha="center", va="center", fontsize=7, color="black")

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Spearman rho", fontsize=9)
    ax.set_title("Panel E: Score correlations (donor-level)", fontsize=12)
    save(fig, "SupFig8E_correlation_heatmap")


def main() -> None:
    ensure_dir()
    plot_panel_a()
    plot_panel_b()
    plot_panel_c()
    plot_panel_d()
    plot_panel_e()


if __name__ == "__main__":
    main()
