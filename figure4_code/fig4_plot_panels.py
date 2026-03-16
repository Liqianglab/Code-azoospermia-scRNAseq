#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
SRC = BASE / "figure4_source_data"
OUT = BASE / "figure4_code" / "figure4_plots"

GROUP_ORDER = ["Ctrl", "OA", "AZFc_Del", "iNOA_B", "iNOA_S", "KS"]
GROUP_COLORS = {
    "Ctrl": "#1B9E77",
    "OA": "#E41A1C",
    "AZFc_Del": "#D95F02",
    "iNOA_B": "#7570B3",
    "iNOA_S": "#1F78B4",
    "KS": "#A6761D",
}
STAGE_COLORS = {"Stage_a": "#4daf4a", "Stage_b": "#377eb8", "Stage_c": "#984ea3"}


def ensure_out() -> None:
    OUT.mkdir(parents=True, exist_ok=True)


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def _ordered(items: list[str], preferred: list[str]) -> list[str]:
    return [x for x in preferred if x in items] + [x for x in items if x not in preferred]


def panel_a(ax: plt.Axes) -> None:
    cells = pd.read_csv(SRC / "Fig4A_ST_cells_trajectory.csv")
    lines = pd.read_csv(SRC / "Fig4A_trajectory_line_points.csv")

    stages = [s for s in ["Stage_a", "Stage_b", "Stage_c"] if s in cells["stage_label_candidate"].astype(str).unique()]
    rng = np.random.default_rng(0)

    for stage in stages:
        sub = cells[cells["stage_label_candidate"].astype(str) == stage]
        if len(sub) > 2500:
            sub = sub.iloc[rng.choice(len(sub), 2500, replace=False)]
        ax.scatter(
            sub["traj_x"],
            sub["traj_y"],
            s=5,
            alpha=0.45,
            color=STAGE_COLORS.get(stage, "#999999"),
            edgecolors="none",
            label=stage,
        )

    if {"line_id", "traj_x", "traj_y"}.issubset(lines.columns):
        for _, sub in lines.groupby("line_id", observed=False):
            ax.plot(sub["traj_x"], sub["traj_y"], color="black", linewidth=0.8, alpha=0.45)

    ax.set_title("A ST trajectory", loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel("traj_x")
    ax.set_ylabel("traj_y")
    ax.legend(frameon=False, fontsize=8)


def panel_b(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig4B_heatmap_matrix_wide.csv")
    cols = [c for c in ["Stage_a", "Stage_b", "Stage_c", "ST1", "ST2", "ST3"] if c in df.columns]
    mat = df[cols].to_numpy(dtype=float)
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)

    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(cols, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df["gene"].astype(str).tolist(), fontsize=7)
    ax.set_title("B ST-stage heatmap", loc="left", fontsize=11, fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Z-score", fontsize=8)


def panel_d(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig4D_stage_ratio_by_group.csv")
    pivot = (
        df.pivot(index="group", columns="stage_label_candidate", values="ratio_in_group")
        .fillna(0)
        .sort_index()
    )
    groups = _ordered(list(pivot.index), GROUP_ORDER)
    pivot = pivot.loc[groups]

    stages = [s for s in ["Stage_a", "Stage_b", "Stage_c"] if s in pivot.columns] + [
        s for s in pivot.columns if s not in {"Stage_a", "Stage_b", "Stage_c"}
    ]

    bottom = np.zeros(len(groups))
    for st in stages:
        vals = pivot[st].to_numpy(dtype=float)
        ax.bar(
            np.arange(len(groups)),
            vals,
            bottom=bottom,
            width=0.8,
            label=st,
            color=STAGE_COLORS.get(st, "#aaaaaa"),
            edgecolor="white",
            linewidth=0.3,
        )
        bottom += vals

    ax.set_xticks(np.arange(len(groups)))
    ax.set_xticklabels(groups, rotation=35, ha="right")
    ax.set_ylabel("Fraction")
    ax.set_title("D Stage ratio by group", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8)


def panel_f(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig4F_regulon_activity_matrix_wide.csv")
    z_cols = [c for c in df.columns if c.endswith("_activity_zscore")]
    group_cols = []
    for g in GROUP_ORDER:
        col = f"{g}_activity_zscore"
        if col in z_cols:
            group_cols.append(col)
    for c in z_cols:
        if c not in group_cols:
            group_cols.append(c)

    mat = df[group_cols].to_numpy(dtype=float)
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)

    xlabels = [c.replace("_activity_zscore", "") for c in group_cols]
    ax.set_xticks(np.arange(len(xlabels)))
    ax.set_xticklabels(xlabels, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df["regulon_label"].astype(str).tolist(), fontsize=7)
    ax.set_title("F Regulon activity (z-score)", loc="left", fontsize=11, fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Z-score", fontsize=8)


def panel_g(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig4G_BTB_score_cells.csv")
    groups = _ordered(sorted(df["group"].astype(str).unique().tolist()), GROUP_ORDER)
    arr = [pd.to_numeric(df.loc[df["group"] == g, "btb_score"], errors="coerce").dropna().to_numpy() for g in groups]

    bp = ax.boxplot(arr, patch_artist=True, showfliers=False)
    for i, box in enumerate(bp["boxes"]):
        g = groups[i]
        box.set_facecolor(GROUP_COLORS.get(g, "#cccccc"))
        box.set_alpha(0.5)
        box.set_edgecolor("#333333")

    rng = np.random.default_rng(1)
    for i, g in enumerate(groups):
        vals = pd.to_numeric(df.loc[df["group"] == g, "btb_score"], errors="coerce").dropna()
        if vals.empty:
            continue
        take = vals.sample(min(900, len(vals)), random_state=2)
        x = i + 1 + rng.uniform(-0.12, 0.12, size=len(take))
        ax.scatter(x, take.to_numpy(), s=4, alpha=0.22, color=GROUP_COLORS.get(g, "#777777"), edgecolors="none")

    ax.set_xticks(np.arange(1, len(groups) + 1))
    ax.set_xticklabels(groups, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("BTB score")
    ax.set_title("G BTB integrity score", loc="left", fontsize=11, fontweight="bold")


def panel_note(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(
        0.02,
        0.9,
        "Panels C/E/H/I are image or marker overlays\nfrom source-data notes and are not redrawn\nin this public plotting script.",
        ha="left",
        va="top",
        fontsize=10,
    )


def main() -> None:
    ensure_out()

    fig, axes = plt.subplots(2, 3, figsize=(19, 11))
    panel_a(axes[0, 0])
    panel_b(axes[0, 1])
    panel_d(axes[0, 2])
    panel_f(axes[1, 0])
    panel_g(axes[1, 1])
    panel_note(axes[1, 2])

    for ax in axes.ravel():
        if ax.has_data():
            ax.grid(True, axis="y", linestyle="--", alpha=0.25)

    fig.suptitle("Figure 4 (source-data reconstruction)", fontsize=18, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    save(fig, "Figure4_main")


if __name__ == "__main__":
    main()
