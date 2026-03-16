#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
SRC = BASE / "fig7_source_data"
OUT = BASE / "fig7_code" / "fig7_plots"

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
    df = pd.read_csv(SRC / "Fig7A_scmeta_diff_selected_pathways_long.csv")
    df = df.copy()
    if "source_table" in df.columns and "with_OA" in set(df["source_table"].astype(str)):
        df = df[df["source_table"].astype(str) == "with_OA"].copy()

    mat = (
        df.groupby(["pathway", "group"], as_index=False)["score"].mean()
        .pivot(index="pathway", columns="group", values="score")
        .fillna(0)
    )
    cols = _ordered(list(mat.columns), GROUP_ORDER)
    mat = mat[cols]

    im = ax.imshow(mat.to_numpy(dtype=float), aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(cols, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(mat.index)))
    ax.set_yticklabels(mat.index.tolist(), fontsize=7)
    ax.set_title("A scMetabolism selected pathways", loc="left", fontsize=11, fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Score", fontsize=8)


def panel_b(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig7B_diagram_gene_log2fc.csv")
    sub = df.copy().sort_values("log2fc", ascending=True)
    if len(sub) > 26:
        idx = np.r_[np.arange(13), np.arange(len(sub) - 13, len(sub))]
        sub = sub.iloc[idx].sort_values("log2fc", ascending=True)

    y = np.arange(len(sub))
    colors = ["#377eb8" if v < 0 else "#e41a1c" for v in pd.to_numeric(sub["log2fc"], errors="coerce").fillna(0)]
    ax.barh(y, sub["log2fc"], color=colors, alpha=0.85)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["gene"].astype(str).tolist(), fontsize=7)
    ax.set_xlabel("log2FC (disease vs ctrl)")
    ax.set_title("B Diagram genes", loc="left", fontsize=11, fontweight="bold")


def _box_with_points(ax: plt.Axes, df: pd.DataFrame, title: str) -> None:
    groups = _ordered(sorted(df["group"].astype(str).unique().tolist()), GROUP_ORDER)
    arr = [pd.to_numeric(df.loc[df["group"] == g, "score"], errors="coerce").dropna().to_numpy() for g in groups]

    bp = ax.boxplot(arr, patch_artist=True, showfliers=False)
    for i, box in enumerate(bp["boxes"]):
        g = groups[i]
        box.set_facecolor(GROUP_COLORS.get(g, "#cccccc"))
        box.set_alpha(0.5)
        box.set_edgecolor("#333333")

    rng = np.random.default_rng(0)
    for i, g in enumerate(groups):
        vals = pd.to_numeric(df.loc[df["group"] == g, "score"], errors="coerce").dropna()
        if vals.empty:
            continue
        take = vals.sample(min(900, len(vals)), random_state=1)
        x = i + 1 + rng.uniform(-0.12, 0.12, size=len(take))
        ax.scatter(x, take.to_numpy(), s=4, alpha=0.2, color=GROUP_COLORS.get(g, "#777777"), edgecolors="none")

    ax.set_xticks(np.arange(1, len(groups) + 1))
    ax.set_xticklabels(groups, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("Score")
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")


def panel_c(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig7C_ST_lactate_glycolysis_cells.csv")
    _box_with_points(ax, df, "C ST lactate-glycolysis score")


def panel_d(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig7D_late_primary_SPCs_oxphos_cells.csv")
    _box_with_points(ax, df, "D Late primary SPC OXPHOS score")


def panel_e(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig7E_dotplot_summary_all.csv")
    groups = _ordered(sorted(df["group"].astype(str).unique().tolist()), GROUP_ORDER)
    label = df["cell_type"].astype(str) + " | " + df["gene"].astype(str)
    df = df.assign(label=label)

    labels = (
        df.groupby("label", as_index=False)["mean_log1p"]
        .mean()
        .sort_values("mean_log1p", ascending=False)["label"]
        .tolist()
    )
    if len(labels) > 18:
        labels = labels[:18]
    sub = df[df["label"].isin(labels)].copy()

    x_map = {g: i for i, g in enumerate(groups)}
    y_map = {lb: i for i, lb in enumerate(labels)}
    x = sub["group"].map(x_map)
    y = sub["label"].map(y_map)
    size = pd.to_numeric(sub["pct_expressed"], errors="coerce").fillna(0).clip(lower=0) * 2.2
    color = pd.to_numeric(sub["mean_log1p"], errors="coerce").fillna(0)

    sc = ax.scatter(x, y, s=size, c=color, cmap="viridis", edgecolors="black", linewidths=0.2)
    ax.set_xticks(np.arange(len(groups)))
    ax.set_xticklabels(groups, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_title("E Lactate shuttle dotplot summary", loc="left", fontsize=11, fontweight="bold")

    cbar = plt.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Mean log1p expression", fontsize=8)


def panel_note(ax: plt.Axes) -> None:
    ax.axis("off")
    note = pd.read_csv(SRC / "Fig7G_note.csv")
    text = str(note.iloc[0].to_dict()) if not note.empty else "Panel F/G are image-note panels."
    ax.text(0.02, 0.9, "F/G notes:\n" + text[:600], ha="left", va="top", fontsize=9)


def main() -> None:
    ensure_out()

    fig, axes = plt.subplots(2, 3, figsize=(19, 11))
    panel_a(axes[0, 0])
    panel_b(axes[0, 1])
    panel_c(axes[0, 2])
    panel_d(axes[1, 0])
    panel_e(axes[1, 1])
    panel_note(axes[1, 2])

    for ax in axes.ravel():
        if ax.has_data():
            ax.grid(True, axis="y", linestyle="--", alpha=0.25)

    fig.suptitle("Figure 7 (source-data reconstruction)", fontsize=18, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    save(fig, "Figure7_main")


if __name__ == "__main__":
    main()
