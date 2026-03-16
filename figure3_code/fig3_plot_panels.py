#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
SRC = BASE / "figure3_source_data"
OUT = BASE / "figure3_code" / "figure3_plots"

GROUP_ORDER = ["Ctrl", "OA", "AZFc_Del", "iNOA_B", "iNOA_S", "KS"]
GROUP_COLORS = {
    "Ctrl": "#1B9E77",
    "OA": "#E41A1C",
    "AZFc_Del": "#D95F02",
    "iNOA_B": "#7570B3",
    "iNOA_S": "#1F78B4",
    "KS": "#A6761D",
}


def ensure_out() -> None:
    OUT.mkdir(parents=True, exist_ok=True)


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def _ordered(items: list[str], preferred: list[str]) -> list[str]:
    return [x for x in preferred if x in items] + [x for x in items if x not in preferred]


def panel_a(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig3A_somatic_cell_ratio_by_group.csv")
    pivot = (
        df.pivot(index="group", columns="somatic_celltype", values="ratio_in_group")
        .fillna(0)
        .sort_index()
    )
    groups = _ordered(list(pivot.index), GROUP_ORDER)
    pivot = pivot.loc[groups]

    bottom = np.zeros(len(pivot))
    for celltype in pivot.columns:
        vals = pivot[celltype].to_numpy(dtype=float)
        ax.bar(
            np.arange(len(groups)),
            vals,
            bottom=bottom,
            width=0.8,
            label=celltype,
            linewidth=0.3,
            edgecolor="white",
        )
        bottom += vals

    ax.set_xticks(np.arange(len(groups)))
    ax.set_xticklabels(groups, rotation=35, ha="right")
    ax.set_ylabel("Fraction")
    ax.set_title("A Somatic composition by group", loc="left", fontsize=11, fontweight="bold")
    ax.legend(fontsize=7, frameon=False, ncol=2)


def panel_b(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig3B_radar_matrix_group_by_celltype.csv")
    groups = _ordered(df["disease"].astype(str).tolist(), GROUP_ORDER)
    mat = df.set_index("disease").loc[groups]
    val = mat.to_numpy(dtype=float)

    im = ax.imshow(val, aspect="auto", cmap="YlGnBu")
    ax.set_xticks(np.arange(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(groups)))
    ax.set_yticklabels(groups, fontsize=8)
    ax.set_title("B Somatic cell counts matrix", loc="left", fontsize=11, fontweight="bold")
    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Cell count", fontsize=8)


def panel_c(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig3C_bubble_top10_pathways_by_group.csv")
    if "panel_celltype" in df.columns and not df["panel_celltype"].empty:
        target = df["panel_celltype"].astype(str).value_counts().index[0]
        df = df[df["panel_celltype"].astype(str) == target].copy()
    else:
        target = "all"

    groups = _ordered(sorted(df["group"].astype(str).unique().tolist()), GROUP_ORDER)
    pathways = (
        df.groupby("pathway_readable", as_index=False)["neglog10_fdr"].max()
        .sort_values("neglog10_fdr", ascending=False)
        .head(12)["pathway_readable"]
        .tolist()
    )
    sub = df[df["pathway_readable"].isin(pathways)].copy()
    pathway_order = [p for p in pathways if p in sub["pathway_readable"].tolist()]

    x_map = {g: i for i, g in enumerate(groups)}
    y_map = {p: i for i, p in enumerate(pathway_order)}

    x = sub["group"].map(x_map)
    y = sub["pathway_readable"].map(y_map)
    size = (pd.to_numeric(sub["neglog10_fdr"], errors="coerce").fillna(0).clip(lower=0) + 0.1) * 28
    color = pd.to_numeric(sub["NES"], errors="coerce").fillna(0)

    sc = ax.scatter(x, y, s=size, c=color, cmap="RdBu_r", vmin=-3, vmax=3, edgecolors="black", linewidths=0.2)
    ax.set_xticks(np.arange(len(groups)))
    ax.set_xticklabels(groups, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(pathway_order)))
    ax.set_yticklabels(pathway_order, fontsize=7)
    ax.set_title(f"C Top pathways bubble ({target})", loc="left", fontsize=11, fontweight="bold")
    cbar = plt.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("NES", fontsize=8)


def _box_with_points(ax: plt.Axes, df: pd.DataFrame, y_col: str, title: str) -> None:
    groups = _ordered(sorted(df["group"].astype(str).unique().tolist()), GROUP_ORDER)
    arr = [pd.to_numeric(df.loc[df["group"] == g, y_col], errors="coerce").dropna().to_numpy() for g in groups]

    bp = ax.boxplot(arr, patch_artist=True, showfliers=False)
    for i, box in enumerate(bp["boxes"]):
        g = groups[i]
        box.set_facecolor(GROUP_COLORS.get(g, "#cccccc"))
        box.set_alpha(0.45)
        box.set_edgecolor("#333333")

    rng = np.random.default_rng(0)
    for i, g in enumerate(groups):
        vals = pd.to_numeric(df.loc[df["group"] == g, y_col], errors="coerce").dropna()
        if vals.empty:
            continue
        take = vals.sample(min(800, len(vals)), random_state=1)
        x = i + 1 + rng.uniform(-0.12, 0.12, size=len(take))
        ax.scatter(x, take.to_numpy(), s=4, alpha=0.25, color=GROUP_COLORS.get(g, "#777777"), edgecolors="none")

    ax.set_xticks(np.arange(1, len(groups) + 1))
    ax.set_xticklabels(groups, rotation=35, ha="right", fontsize=8)
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")


def panel_d(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig3D_cytokine_score_cells.csv")
    _box_with_points(ax, df, "score", "D Cytokine score")
    ax.set_ylabel("Score")


def panel_e(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig3E_sasp_score_cells.csv")
    _box_with_points(ax, df, "score", "E SASP score")
    ax.set_ylabel("Score")


def panel_f(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig3F_heatmap_matrix_wide.csv")
    cols = [c for c in GROUP_ORDER if c in df.columns]
    if not cols:
        cols = [c for c in df.columns if c != "gene"]
    mat = df[cols].to_numpy(dtype=float)

    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(cols, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df["gene"].astype(str).tolist(), fontsize=7)
    ax.set_title("F Ubiquitination heatmap", loc="left", fontsize=11, fontweight="bold")
    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Z-score", fontsize=8)


def main() -> None:
    ensure_out()

    fig, axes = plt.subplots(2, 3, figsize=(19, 11))
    panel_a(axes[0, 0])
    panel_b(axes[0, 1])
    panel_c(axes[0, 2])
    panel_d(axes[1, 0])
    panel_e(axes[1, 1])
    panel_f(axes[1, 2])

    for ax in axes.ravel():
        ax.grid(True, axis="y", linestyle="--", alpha=0.25)

    fig.suptitle("Figure 3 (source-data reconstruction)", fontsize=18, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    save(fig, "Figure3_main")


if __name__ == "__main__":
    main()
