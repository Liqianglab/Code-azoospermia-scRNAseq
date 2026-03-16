#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
SRC = BASE / "fig8_source_data"
OUT = BASE / "fig8_code" / "fig8_plots"

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
    df = pd.read_csv(SRC / "Fig8A_ctrl_oa_umap_cells.csv")
    sub = df[df["group"].astype(str).isin(["Ctrl", "OA"])].copy()
    if len(sub) > 22000:
        sub = sub.sample(22000, random_state=1)

    for g in ["Ctrl", "OA"]:
        d = sub[sub["group"] == g]
        ax.scatter(d["UMAP1"], d["UMAP2"], s=4, alpha=0.35, color=GROUP_COLORS[g], edgecolors="none", label=g)

    ax.set_title("A Ctrl vs OA germ UMAP", loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(frameon=False, fontsize=8)


def panel_b(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig8B_germ_composition_percent_long.csv")
    pivot = df.pivot(index="group", columns="cluster", values="percent").fillna(0)
    groups = _ordered(list(pivot.index), ["Ctrl", "OA"])
    pivot = pivot.loc[groups]

    bottom = np.zeros(len(groups))
    for cluster in pivot.columns:
        vals = pivot[cluster].to_numpy(dtype=float)
        ax.bar(np.arange(len(groups)), vals, bottom=bottom, width=0.82, label=cluster, edgecolor="white", linewidth=0.3)
        bottom += vals

    ax.set_xticks(np.arange(len(groups)))
    ax.set_xticklabels(groups)
    ax.set_ylabel("Fraction")
    ax.set_title("B Germ composition", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=7, ncol=2)


def panel_cd(ax: plt.Axes) -> None:
    up = pd.read_csv(SRC / "Fig8C_up_in_OA_genes_from_diff_long.csv")
    down = pd.read_csv(SRC / "Fig8D_down_in_OA_genes_from_diff_long.csv")
    up_cnt = up.groupby("subtype", as_index=False).size().rename(columns={"size": "n_genes"})
    up_cnt["direction"] = "Up in OA"
    dn_cnt = down.groupby("subtype", as_index=False).size().rename(columns={"size": "n_genes"})
    dn_cnt["direction"] = "Down in OA"
    df = pd.concat([up_cnt, dn_cnt], ignore_index=True)

    subtypes = sorted(df["subtype"].astype(str).unique().tolist())
    x = np.arange(len(subtypes))
    w = 0.38

    up_vals = [int(df[(df["subtype"] == s) & (df["direction"] == "Up in OA")]["n_genes"].sum()) for s in subtypes]
    dn_vals = [int(df[(df["subtype"] == s) & (df["direction"] == "Down in OA")]["n_genes"].sum()) for s in subtypes]

    ax.bar(x - w / 2, up_vals, width=w, color="#d73027", label="Up in OA")
    ax.bar(x + w / 2, dn_vals, width=w, color="#4575b4", label="Down in OA")

    ax.set_xticks(x)
    ax.set_xticklabels(subtypes, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("DEG count")
    ax.set_title("C-D OA differential genes by subtype", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8)


def panel_e(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig8E_GO_selected_combined.csv")
    df = df[df["found"].astype(str).str.lower().isin(["true", "1"])].copy()
    if df.empty:
        ax.axis("off")
        return

    df = df.sort_values(["direction_in_OA", "term_order"])
    y = np.arange(len(df))
    vals = pd.to_numeric(df["neg_log10_pvalue"], errors="coerce").fillna(0).to_numpy()
    colors = ["#4575b4" if "Down" in d else "#d73027" for d in df["direction_in_OA"].astype(str)]

    ax.barh(y, vals, color=colors, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(df["Description"].astype(str).tolist(), fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("-log10(pvalue)")
    ax.set_title("E Selected GO terms", loc="left", fontsize=11, fontweight="bold")


def panel_g(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig8G_candidate_ctrl_oa_expression_zscore.csv")
    if "gene_order" in df.columns:
        df = df.sort_values("gene_order")

    mat = df[["Ctrl_z", "OA_z"]].to_numpy(dtype=float)
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Ctrl", "OA"])
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df["gene"].astype(str).tolist(), fontsize=7)
    ax.set_title("G Candidate gene z-score", loc="left", fontsize=11, fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Z-score", fontsize=8)


def panel_j(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig8J_score_cells_combined_long.csv")
    groups = _ordered(sorted(df["group"].astype(str).unique().tolist()), GROUP_ORDER)
    types = sorted(df["score_type"].astype(str).unique().tolist())

    width = 0.38
    xbase = np.arange(len(groups))

    for ti, st in enumerate(types):
        arr = [
            pd.to_numeric(df.loc[(df["group"] == g) & (df["score_type"] == st), "score"], errors="coerce")
            .dropna()
            .to_numpy()
            for g in groups
        ]
        pos = xbase + (ti - (len(types) - 1) / 2) * width
        bp = ax.boxplot(arr, positions=pos, widths=width * 0.85, patch_artist=True, showfliers=False)
        color = ["#984ea3", "#4daf4a", "#377eb8", "#ff7f00"][ti % 4]
        for box in bp["boxes"]:
            box.set_facecolor(color)
            box.set_alpha(0.35)
            box.set_edgecolor("#333333")

    ax.set_xticks(xbase)
    ax.set_xticklabels(groups, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("Score")
    ax.set_title("J LC mature/immature scores", loc="left", fontsize=11, fontweight="bold")

    handles = [
        plt.Line2D([0], [0], color="none", marker="s", markerfacecolor=["#984ea3", "#4daf4a", "#377eb8", "#ff7f00"][i % 4], markersize=8, label=t)
        for i, t in enumerate(types)
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper right")


def panel_i(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Fig8I_somatic_composition_percent_long.csv")
    pivot = df.pivot(index="group", columns="panel_cluster", values="percent").fillna(0)
    groups = _ordered(list(pivot.index), ["Ctrl", "OA"])
    pivot = pivot.loc[groups]

    bottom = np.zeros(len(groups))
    for cluster in pivot.columns:
        vals = pivot[cluster].to_numpy(dtype=float)
        ax.bar(np.arange(len(groups)), vals, bottom=bottom, width=0.82, label=cluster, edgecolor="white", linewidth=0.3)
        bottom += vals

    ax.set_xticks(np.arange(len(groups)))
    ax.set_xticklabels(groups)
    ax.set_ylabel("Fraction")
    ax.set_title("I Somatic composition", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=7)


def panel_note(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(
        0.02,
        0.9,
        "Panels F/H are source-note and auxiliary\nviews in this public reconstruction script.",
        ha="left",
        va="top",
        fontsize=10,
    )


def main() -> None:
    ensure_out()

    fig, axes = plt.subplots(3, 3, figsize=(20, 15))
    panel_a(axes[0, 0])
    panel_b(axes[0, 1])
    panel_cd(axes[0, 2])
    panel_e(axes[1, 0])
    panel_g(axes[1, 1])
    panel_i(axes[1, 2])
    panel_j(axes[2, 0])
    panel_note(axes[2, 1])
    axes[2, 2].axis("off")

    for ax in axes.ravel():
        if ax.has_data():
            ax.grid(True, axis="y", linestyle="--", alpha=0.25)

    fig.suptitle("Figure 8 (source-data reconstruction)", fontsize=18, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    save(fig, "Figure8_main")


if __name__ == "__main__":
    main()
