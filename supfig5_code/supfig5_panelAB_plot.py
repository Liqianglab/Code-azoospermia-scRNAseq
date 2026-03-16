#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import gridspec


SCORE_ORDER = [
    ("Luteinizing_hormone", "Luteinizing hormone Score"),
    ("Steroid_metabolism", "Steroid Metabolism Score"),
    ("Glucocorticoid_receptor_pathway", "Glucocorticoid receptor pathway score"),
]
STAGE_ORDER = ["LC_a", "LC_b", "LC_c"]
STAGE_COLORS = ["#1f77b4", "#2ca02c", "#ff4d4d"]
GROUP_ORDER = ["AZFc", "iNOA_B", "iNOA_S", "KS"]
DIRECTION_ORDER = ["Up", "Down"]
GROUP_COLORS = {
    "AZFc": "#1f77b4",
    "iNOA_B": "#d62728",
    "iNOA_S": "#2ca02c",
    "KS": "#9467bd",
}


def read_csv_rows(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def to_float(value):
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def plot_panel_a(ax_list, panel_a_rows):
    grouped = defaultdict(list)
    for row in panel_a_rows:
        score_id = str(row.get("score_id", "")).strip()
        stage = str(row.get("stage", "")).strip()
        value = to_float(row.get("score_value"))
        if value is None:
            value = to_float(row.get("score_proxy"))
        if score_id and stage and value is not None:
            grouped[(score_id, stage)].append(value)

    for ax, (score_id, score_label) in zip(ax_list, SCORE_ORDER):
        data = [grouped.get((score_id, stage), []) for stage in STAGE_ORDER]

        violin = ax.violinplot(
            data,
            positions=[1, 2, 3],
            widths=0.8,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )

        for body, color in zip(violin["bodies"], STAGE_COLORS):
            body.set_facecolor(color)
            body.set_edgecolor("black")
            body.set_alpha(0.9)

        box = ax.boxplot(
            data,
            positions=[1, 2, 3],
            widths=0.13,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.0},
            whiskerprops={"linewidth": 0.8, "color": "black"},
            capprops={"linewidth": 0.8, "color": "black"},
            boxprops={"linewidth": 0.8, "color": "black"},
        )
        for patch in box["boxes"]:
            patch.set_facecolor("white")
            patch.set_alpha(0.85)

        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(STAGE_ORDER, fontsize=9)
        ax.set_title(score_label, fontsize=11)
        ax.set_ylabel("Score", fontsize=9)
        ax.tick_params(axis="y", labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)


def plot_panel_b(ax_grid, panel_b_top_rows):
    grouped = defaultdict(list)
    for row in panel_b_top_rows:
        direction = str(row.get("direction", "")).strip()
        group = str(row.get("group", "")).strip()
        rank = int(float(row.get("rank_in_group", "0") or 0))
        desc = str(row.get("Description", "")).strip()
        logp = to_float(row.get("minus_log10_pvalue"))
        if direction and group and rank > 0 and desc and logp is not None:
            grouped[(direction, group)].append(
                {
                    "rank": rank,
                    "Description": desc,
                    "minus_log10_pvalue": logp,
                }
            )

    for row_idx, group in enumerate(GROUP_ORDER):
        for col_idx, direction in enumerate(DIRECTION_ORDER):
            ax = ax_grid[row_idx][col_idx]
            rows = grouped.get((direction, group), [])
            rows = sorted(rows, key=lambda r: r["rank"])[:15]

            # Put most significant term at top
            labels = [r["Description"] for r in rows][::-1]
            values = [r["minus_log10_pvalue"] for r in rows][::-1]

            y = list(range(len(labels)))
            ax.barh(y, values, color=GROUP_COLORS.get(group, "#777777"), edgecolor="none")
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=7)
            ax.tick_params(axis="x", labelsize=8)
            ax.set_xlabel("-log10(p-value)", fontsize=9)
            ax.set_title(group, fontsize=11)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            if row_idx == 0:
                ax.text(
                    0.0,
                    1.18,
                    direction,
                    transform=ax.transAxes,
                    ha="left",
                    va="bottom",
                    fontsize=13,
                    fontweight="bold",
                )


def main():
    code_dir = Path(__file__).resolve().parent
    source_dir = code_dir.parent / "supfig5_source_data"
    out_dir = code_dir / "supfig5_plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    panel_a_raw_path = source_dir / "SupFig5A_score_cells.csv"
    panel_a_proxy_path = source_dir / "SupFig5A_score_proxy_cells.csv"
    panel_a_rows = read_csv_rows(panel_a_raw_path if panel_a_raw_path.exists() else panel_a_proxy_path)
    panel_b_top_rows = read_csv_rows(source_dir / "SupFig5B_GO_updown_top15_by_group.csv")

    # Panel A only
    fig_a, axes_a = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    plot_panel_a(axes_a, panel_a_rows)
    fig_a.savefig(out_dir / "SupFig5A_violin.svg", dpi=300)
    fig_a.savefig(out_dir / "SupFig5A_violin.png", dpi=300)
    fig_a.savefig(out_dir / "SupFig5A_violin_proxy.svg", dpi=300)
    fig_a.savefig(out_dir / "SupFig5A_violin_proxy.png", dpi=300)
    plt.close(fig_a)

    # Panel B only
    fig_b, axes_b = plt.subplots(4, 2, figsize=(14, 18), constrained_layout=True)
    plot_panel_b(axes_b, panel_b_top_rows)
    fig_b.savefig(out_dir / "SupFig5B_GO_updown_top15.svg", dpi=300)
    fig_b.savefig(out_dir / "SupFig5B_GO_updown_top15.png", dpi=300)
    plt.close(fig_b)

    # Combined proxy layout
    fig = plt.figure(figsize=(14, 22), constrained_layout=False)
    gs = gridspec.GridSpec(5, 2, figure=fig, height_ratios=[1.1, 1, 1, 1, 1], hspace=0.55, wspace=0.95)

    gs_top = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0, :], wspace=0.4)
    axes_top = [fig.add_subplot(gs_top[0, i]) for i in range(3)]
    plot_panel_a(axes_top, panel_a_rows)

    axes_bottom = []
    for i in range(4):
        row_axes = []
        for j in range(2):
            row_axes.append(fig.add_subplot(gs[i + 1, j]))
        axes_bottom.append(row_axes)
    plot_panel_b(axes_bottom, panel_b_top_rows)

    fig.suptitle("Supplementary Figure 5 (reconstruction)", fontsize=16, y=0.995)
    fig.savefig(out_dir / "SupFig5AB_layout.svg", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "SupFig5AB_layout.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "SupFig5AB_proxy_layout.svg", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "SupFig5AB_proxy_layout.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("SupFig5 plots generated in:", out_dir)


if __name__ == "__main__":
    main()
