#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
SRC = BASE / "figure5_source_data"
OUT = BASE / "figure5_code" / "figure5_plots"

GROUP_ORDER = ["Ctrl", "OA", "AZFc_Del", "iNOA_B", "iNOA_S", "KS"]
GROUP_COLORS = {
    "Ctrl": "#1B9E77",
    "OA": "#E41A1C",
    "AZFc_Del": "#D95F02",
    "iNOA_B": "#7570B3",
    "iNOA_S": "#1F78B4",
    "KS": "#A6761D",
}
STAGE_COLORS = {"Stage_a": "#66c2a5", "Stage_b": "#fc8d62", "Stage_c": "#8da0cb"}


def ensure_out() -> None:
    OUT.mkdir(parents=True, exist_ok=True)


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def _ordered(items: list[str], preferred: list[str]) -> list[str]:
    return [x for x in preferred if x in items] + [x for x in items if x not in preferred]


def panel_a(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig5A_LC_cells_umap_stage_group.csv")
    stages = [s for s in ["Stage_a", "Stage_b", "Stage_c"] if s in df["stage"].astype(str).unique()]
    rng = np.random.default_rng(0)

    for stage in stages:
        sub = df[df["stage"].astype(str) == stage]
        if len(sub) > 2600:
            sub = sub.iloc[rng.choice(len(sub), 2600, replace=False)]
        ax.scatter(
            sub["UMAP1"],
            sub["UMAP2"],
            s=5,
            alpha=0.45,
            color=STAGE_COLORS.get(stage, "#999999"),
            edgecolors="none",
            label=stage,
        )

    ax.set_title("A LC trajectory UMAP", loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(frameon=False, fontsize=8)


def panel_b(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig5B_heatmap_matrix_wide.csv")
    cols = [c for c in ["Stage_a", "Stage_b", "Stage_c", "LC1", "LC2", "LC3"] if c in df.columns]
    mat = df[cols].to_numpy(dtype=float)

    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(cols, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df["gene"].astype(str).tolist(), fontsize=7)
    ax.set_title("B LC-stage heatmap", loc="left", fontsize=11, fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Z-score", fontsize=8)


def panel_d(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig5D_stage_ratio_by_group.csv")
    pivot = df.pivot(index="group", columns="stage", values="ratio").fillna(0)
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
            color=STAGE_COLORS.get(st, "#aaaaaa"),
            width=0.8,
            edgecolor="white",
            linewidth=0.3,
            label=st,
        )
        bottom += vals

    ax.set_xticks(np.arange(len(groups)))
    ax.set_xticklabels(groups, rotation=35, ha="right")
    ax.set_ylabel("Fraction")
    ax.set_title("D LC stage ratio", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8)


def panel_f(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig5F_GO_down_vs_ctrl_top5_by_group.csv")
    sub = df.copy().sort_values(["group", "rank_in_group"])
    sub["label"] = sub["group"].astype(str) + " | " + sub["Description"].astype(str)
    top = sub.head(20)

    y = np.arange(len(top))
    val = pd.to_numeric(top["minus_log10_p_adjust"], errors="coerce").fillna(0).to_numpy()
    colors = [GROUP_COLORS.get(g, "#999999") for g in top["group"].astype(str)]

    ax.barh(y, val, color=colors, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(top["label"].tolist(), fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("-log10(p.adjust)")
    ax.set_title("F GO terms (down vs Ctrl)", loc="left", fontsize=11, fontweight="bold")


def panel_i(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig5I_hormone_values.csv")
    metric_order = list(df["Metric"].astype(str).drop_duplicates())
    group_order = ["Control", "OA", "NOA"]

    width = 0.22
    pos_all = []
    labels = []
    for mi, metric in enumerate(metric_order):
        center = mi * 1.2
        for gi, group in enumerate(group_order):
            vals = pd.to_numeric(
                df.loc[(df["Metric"] == metric) & (df["Group"] == group), "Value"],
                errors="coerce",
            ).dropna()
            if vals.empty:
                continue
            pos = center + (gi - 1) * width
            bp = ax.boxplot([vals.to_numpy()], positions=[pos], widths=width * 0.9, patch_artist=True, showfliers=False)
            color = {"Control": "#1B9E77", "OA": "#E41A1C", "NOA": "#7570B3"}.get(group, "#999999")
            bp["boxes"][0].set_facecolor(color)
            bp["boxes"][0].set_alpha(0.5)
            bp["boxes"][0].set_edgecolor("#333333")
            pos_all.append(pos)
            labels.append((metric, group))

    ax.set_xticks([i * 1.2 for i in range(len(metric_order))])
    ax.set_xticklabels(metric_order, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("Value")
    ax.set_title("I Hormone distribution", loc="left", fontsize=11, fontweight="bold")

    handles = [
        plt.Line2D([0], [0], color="none", marker="s", markerfacecolor="#1B9E77", markersize=8, label="Control"),
        plt.Line2D([0], [0], color="none", marker="s", markerfacecolor="#E41A1C", markersize=8, label="OA"),
        plt.Line2D([0], [0], color="none", marker="s", markerfacecolor="#7570B3", markersize=8, label="NOA"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper right")


def panel_note(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(
        0.02,
        0.9,
        "Panels C/E/G/H use marker overlays or image panels\nfrom source-data notes and are not redrawn\nin this public plotting script.",
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
    panel_i(axes[1, 1])
    panel_note(axes[1, 2])

    for ax in axes.ravel():
        if ax.has_data():
            ax.grid(True, axis="y", linestyle="--", alpha=0.25)

    fig.suptitle("Figure 5 (source-data reconstruction)", fontsize=18, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    save(fig, "Figure5_main")


if __name__ == "__main__":
    main()
