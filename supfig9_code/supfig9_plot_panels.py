#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

BASE = Path(__file__).resolve().parents[1]
SRC = BASE / "supfig9_source_data"
OUT = BASE / "supfig9_code" / "supfig9_plots"

FIG_BG = "#eeeeee"
AX_BG = "#f2f2f2"

COL_A = {"Normal": "#4C78A8", "NOA1": "#F58518", "NOA2": "#54A24B"}
COL_BG = {"Control": "#1B9E77", "AZFa": "#E41A1C", "iNOA": "#6A3D9A", "KS": "#A6761D"}
COL_HJK = {"Normal": "#1B9E77", "AZFa": "#D95F02", "iNOA": "#7570B3", "KS": "#8C6D31"}

GROUP_ORDER_A = ["Normal", "NOA1", "NOA2"]
GROUP_ORDER_BG = ["Control", "AZFa", "iNOA", "KS"]
GROUP_ORDER_HJK = ["Normal", "AZFa", "iNOA", "KS"]

PANEL_BG_META = {
    "B": ("ST maturity", "ST_maturity"),
    "C": ("LC immaturity", "LCimmature_z"),
    "D": ("Somatic cytokine signaling", "Cytokine_z"),
    "E": ("Somatic SASP", "SASP"),
    "F": ("ST BTB integrity", "BTB"),
    "G": ("Germ apoptosis", "GermApoptosis_z"),
}
PANEL_BG_ORDER = ["B", "D", "F", "C", "E", "G"]


def ensure_out() -> None:
    OUT.mkdir(parents=True, exist_ok=True)


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def read_csv_if_exists(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def style_axis(ax: plt.Axes, hide_top_right: bool = True) -> None:
    ax.set_facecolor(AX_BG)
    ax.grid(True, linestyle="--", alpha=0.25, color="#bdbdbd")
    ax.spines["left"].set_color("#555555")
    ax.spines["bottom"].set_color("#555555")
    if hide_top_right:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    else:
        ax.spines["top"].set_color("#555555")
        ax.spines["right"].set_color("#555555")


def _panel_letter(fig: plt.Figure, ax: plt.Axes, letter: str, xpad: float = -0.018, ypad: float = 0.005) -> None:
    pos = ax.get_position()
    fig.text(pos.x0 + xpad, pos.y1 + ypad, letter, fontsize=16, ha="left", va="bottom")


def _corr_text(panel: str) -> str:
    csum = read_csv_if_exists(SRC / "SupFig9HJ_correlation_summary.csv")
    if csum.empty:
        return ""
    row = csum[(csum["panel"] == panel) & (csum["mode"] == "figure_annotation")]
    if row.empty:
        row = csum[(csum["panel"] == panel) & (csum["mode"] == "donor_level")]
    if row.empty:
        return ""
    rho = float(row["spearman_rho"].iloc[0])
    pval = float(row["p_value"].iloc[0])
    ptxt = f"{pval:.3g}" if pval < 0.01 else f"{pval:.4g}"
    return f"Spearman rho={rho:.2f}, P={ptxt}"


def draw_panel_a(ax: plt.Axes) -> bool:
    df = read_csv_if_exists(SRC / "SupFig9A_celllevel_pseudotime_entropy.csv")
    sm = read_csv_if_exists(SRC / "SupFig9A_smoothed_curve.csv")
    st = read_csv_if_exists(SRC / "SupFig9A_stage_ranges.csv")
    ap = read_csv_if_exists(SRC / "SupFig9A_arrest_point.csv")
    if df.empty:
        ax.axis("off")
        return False

    style_axis(ax, hide_top_right=False)

    stage_bg = {
        "SPG (Stem)": ("#f2d7d6", 0.45),
        "Early SPC": ("#d9e7f6", 0.45),
        "Late SPC": ("#e5dcef", 0.45),
        "Round Sperm": ("#dcead7", 0.45),
        "Elongated Sperm": ("#cedfc6", 0.38),
    }
    stage_label = {
        "SPG (Stem)": "SPG",
        "Early SPC": "Early SPC",
        "Late SPC": "Late SPC",
        "Round Sperm": "Round Sperm",
        "Elongated Sperm": "Elongated Sperm",
    }

    if not st.empty and {"plot_stage", "x_min", "x_max"}.issubset(st.columns):
        for _, r in st.sort_values("x_min").iterrows():
            stage = str(r["plot_stage"])
            color, alpha = stage_bg.get(stage, ("#dddddd", 0.25))
            x0 = float(r["x_min"])
            x1 = float(r["x_max"])
            ax.axvspan(x0, x1, color=color, alpha=alpha, lw=0)
            text_y = 1.045
            text_size = 12
            if stage == "Round Sperm":
                text_y = 0.92
                text_size = 14
            elif stage == "Elongated Sperm":
                text_y = 0.80
                text_size = 13
            ax.text((x0 + x1) / 2, text_y, stage_label.get(stage, stage), ha="center", va="top", fontsize=text_size, fontweight="bold")

    for g in GROUP_ORDER_A:
        sub = df[df["group"] == g].copy()
        if sub.empty:
            continue
        sub = sub.sort_values(["pseudotime", "mix_entropy"])
        ax.scatter(
            sub["pseudotime"],
            sub["mix_entropy"],
            s=16,
            alpha=0.7,
            color=COL_A.get(g, "#777777"),
            edgecolors="none",
            zorder=2,
        )

    if not sm.empty and {"pseudotime", "entropy_smooth"}.issubset(sm.columns):
        curve = sm.sort_values("pseudotime")
        ax.plot(curve["pseudotime"], curve["entropy_smooth"], color="black", linewidth=3.0, zorder=4)

    if not ap.empty and "arrest_point_pseudotime" in ap.columns:
        x = float(ap["arrest_point_pseudotime"].iloc[0])
        if np.isfinite(x):
            ax.axvline(x, color="#d62728", linestyle="--", linewidth=2.8, zorder=3)
            ax.text(x - 1.02, 0.00, "Arrest point", color="#d62728", fontsize=15, fontweight="bold")

    ax.set_xlim(0, 10.4)
    ax.set_ylim(-0.08, 1.08)
    ax.set_xlabel("Pseudotime(spermatogenesis progression)", fontsize=16, fontweight="bold", labelpad=4)
    ax.set_ylabel("Group-mixing entropy (bits)", fontsize=16)
    ax.set_title("Developmental divergence(GSE235321)", fontsize=17, fontweight="bold", pad=14)
    ax.tick_params(labelsize=12)

    handles = [
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=COL_A[g], markeredgecolor="none", markersize=5, label=g)
        for g in GROUP_ORDER_A
    ]
    ax.legend(handles=handles, title="Group", frameon=True, fontsize=10, title_fontsize=12, loc="upper right")
    return True


def violin_with_donor(
    ax: plt.Axes,
    cells: pd.DataFrame,
    donor: pd.DataFrame,
    group_order: list[str],
    title: str,
    pvals: pd.DataFrame | None = None,
) -> None:
    style_axis(ax)
    arrays: list[np.ndarray] = []
    present: list[str] = []
    pos_map = {g: i + 1 for i, g in enumerate(group_order)}

    for g in group_order:
        vals = pd.to_numeric(cells.loc[cells["group"] == g, "value"], errors="coerce").dropna().to_numpy()
        if len(vals) == 0:
            continue
        arrays.append(vals)
        present.append(g)

    if arrays:
        parts = ax.violinplot(
            arrays,
            positions=[pos_map[g] for g in present],
            widths=0.82,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )
        for i, body in enumerate(parts["bodies"]):
            g = present[i]
            body.set_facecolor(COL_BG.get(g, "#999999"))
            body.set_edgecolor("#333333")
            body.set_linewidth(1.2)
            body.set_alpha(0.86)

        ax.boxplot(
            arrays,
            positions=[pos_map[g] for g in present],
            widths=0.22,
            patch_artist=True,
            showfliers=False,
            boxprops={"facecolor": "white", "edgecolor": "#222222", "linewidth": 1.1},
            medianprops={"color": "#222222", "linewidth": 1.2},
            whiskerprops={"color": "#222222", "linewidth": 1.0},
            capprops={"color": "#222222", "linewidth": 1.0},
        )

    rng = np.random.default_rng(0)
    for g in group_order:
        med = donor.loc[donor["group"] == g, "median"] if "group" in donor.columns else pd.Series(dtype=float)
        if med.empty:
            continue
        x = pos_map[g] + rng.uniform(-0.07, 0.07, size=len(med))
        ax.scatter(x, med.to_numpy(), s=28, color=COL_BG.get(g, "#999999"), edgecolor="#333333", linewidth=0.7, zorder=4)

    if pvals is not None and not pvals.empty:
        ymin, ymax = ax.get_ylim()
        span = (ymax - ymin) if ymax > ymin else 1.0
        ytxt = ymax + 0.04 * span
        for _, r in pvals.iterrows():
            comp = str(r.get("comparison", ""))
            star = str(r.get("star", ""))
            if " vs Control" not in comp:
                continue
            g = comp.replace(" vs Control", "")
            if g in pos_map:
                ax.text(pos_map[g], ytxt, star, ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.set_ylim(ymin, ymax + 0.14 * span)

    ax.set_xticks([pos_map[g] for g in group_order])
    ax.set_xticklabels(group_order, rotation=38, ha="right", fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.grid(axis="y", linestyle="--", alpha=0.25, color="#bdbdbd")


def draw_panel_bg(ax: plt.Axes, letter: str) -> bool:
    if letter not in PANEL_BG_META:
        ax.axis("off")
        return False
    title, metric = PANEL_BG_META[letter]
    cell_path = SRC / f"SupFig9{letter}_{metric}_cells.csv"
    donor_path = SRC / f"SupFig9{letter}_{metric}_donor_summary.csv"
    pval_path = SRC / f"SupFig9{letter}_{metric}_pvalues.csv"
    if not cell_path.exists() or not donor_path.exists():
        ax.axis("off")
        return False
    cells = pd.read_csv(cell_path)
    donor = pd.read_csv(donor_path)
    pvals = pd.read_csv(pval_path) if pval_path.exists() else pd.DataFrame()
    violin_with_donor(ax, cells, donor, GROUP_ORDER_BG, title, pvals)
    ax.tick_params(labelsize=9)
    return True


def scatter_with_labels(
    ax: plt.Axes,
    df: pd.DataFrame,
    xcol: str,
    ycol: str,
    title: str,
    xlabel: str,
    ylabel: str,
    corr_text: str = "",
) -> None:
    style_axis(ax)
    x = pd.to_numeric(df[xcol], errors="coerce")
    y = pd.to_numeric(df[ycol], errors="coerce")
    df = df.assign(_x=x, _y=y).dropna(subset=["_x", "_y"])

    for g in GROUP_ORDER_HJK:
        sub = df[df["group_plot"] == g]
        if sub.empty:
            continue
        ax.scatter(sub["_x"], sub["_y"], s=52, color=COL_HJK.get(g, "#999999"), edgecolor="#333333", linewidth=0.7, zorder=3)

    if len(df) >= 2:
        coef = np.polyfit(df["_x"], df["_y"], 1)
        xx = np.linspace(float(df["_x"].min()), float(df["_x"].max()), 120)
        yy = coef[0] * xx + coef[1]
        ax.plot(xx, yy, color="#333333", linestyle="--", linewidth=1.2, zorder=2)

    xspan = float(df["_x"].max() - df["_x"].min()) if len(df) else 1.0
    yspan = float(df["_y"].max() - df["_y"].min()) if len(df) else 1.0
    dx = 0.012 * (xspan if xspan > 0 else 1.0)
    dy = 0.015 * (yspan if yspan > 0 else 1.0)
    for _, r in df.iterrows():
        label = str(r.get("sample_code", ""))
        ax.text(float(r["_x"]) + dx, float(r["_y"]) + dy, label, fontsize=10)

    if corr_text:
        ax.text(0.98, 0.97, corr_text, transform=ax.transAxes, ha="right", va="top", fontsize=10)

    ax.set_title(title, fontsize=11)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.tick_params(labelsize=10)


def box_with_points(
    ax: plt.Axes,
    df: pd.DataFrame,
    group_col: str,
    value_col: str,
    group_order: list[str],
    title: str,
    ylabel: str,
) -> None:
    style_axis(ax)
    arrays: list[np.ndarray] = []
    labels: list[str] = []
    for g in group_order:
        vals = pd.to_numeric(df.loc[df[group_col] == g, value_col], errors="coerce").dropna().to_numpy()
        if len(vals) == 0:
            continue
        arrays.append(vals)
        labels.append(g)

    if arrays:
        pos = np.arange(1, len(labels) + 1)
        ax.boxplot(
            arrays,
            positions=pos,
            widths=0.5,
            patch_artist=True,
            showfliers=False,
            boxprops={"facecolor": "#dddddd", "edgecolor": "#777777", "linewidth": 1.1},
            medianprops={"color": "#555555", "linewidth": 1.2},
            whiskerprops={"color": "#777777", "linewidth": 1.0},
            capprops={"color": "#777777", "linewidth": 1.0},
        )
        rng = np.random.default_rng(0)
        for i, g in enumerate(labels):
            vals = pd.to_numeric(df.loc[df[group_col] == g, value_col], errors="coerce").dropna().to_numpy()
            x = np.full(len(vals), pos[i]) + rng.uniform(-0.08, 0.08, size=len(vals))
            ax.scatter(x, vals, s=35, color=COL_HJK.get(g, "#999999"), edgecolor="#333333", linewidth=0.7, zorder=3)
        ax.set_xticks(pos)
        ax.set_xticklabels(labels, rotation=18, ha="right")

    ax.set_title(title, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.tick_params(labelsize=10)
    ax.grid(axis="y", linestyle="--", alpha=0.25, color="#bdbdbd")


def _load_hijk() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    h = read_csv_if_exists(SRC / "SupFig9H_coupling_scatter_donor_points.csv")
    i_df = read_csv_if_exists(SRC / "SupFig9I_st_lactate_export_by_group.csv")
    j = read_csv_if_exists(SRC / "SupFig9J_uncoupling_vs_postmeiotic.csv")
    k_df = read_csv_if_exists(SRC / "SupFig9K_germ_uptake_by_group.csv")
    return h, i_df, j, k_df


def draw_panel_h(ax: plt.Axes, corr_text: str = "") -> bool:
    h, _, _, _ = _load_hijk()
    if h.empty:
        ax.axis("off")
        return False
    scatter_with_labels(
        ax,
        h,
        "st_export_z",
        "germ_uptake_z",
        "Coupling between lactate export and uptake",
        "ST lactate export (z)",
        "Germ lactate uptake/oxidation (z)",
        corr_text=corr_text,
    )
    return True


def draw_panel_i(ax: plt.Axes) -> bool:
    _, i_df, _, _ = _load_hijk()
    if i_df.empty:
        ax.axis("off")
        return False
    box_with_points(ax, i_df, "group_plot", "st_export_z", GROUP_ORDER_HJK, "ST lactate export by group", "ST export score (z)")
    return True


def draw_panel_j(ax: plt.Axes, corr_text: str = "") -> bool:
    _, _, j, _ = _load_hijk()
    if j.empty:
        ax.axis("off")
        return False
    scatter_with_labels(
        ax,
        j,
        "uncoupling_index",
        "postmeiotic_fraction",
        "Uncoupling predicts post-meiotic output",
        "Uncoupling index (export_z - uptake_z)",
        "Spermatid fraction (post-meiotic / germ)",
        corr_text=corr_text,
    )
    return True


def draw_panel_k(ax: plt.Axes) -> bool:
    _, _, _, k_df = _load_hijk()
    if k_df.empty:
        ax.axis("off")
        return False
    box_with_points(
        ax,
        k_df,
        "group_plot",
        "germ_uptake_z",
        GROUP_ORDER_HJK,
        "Germ lactate uptake/oxidation by group",
        "Germ uptake/oxidation score (z)",
    )
    return True


def plot_panel_a() -> None:
    fig, ax = plt.subplots(figsize=(12.5, 6.3), facecolor=FIG_BG)
    draw_panel_a(ax)
    save(fig, "SupFig9A_divergence")


def plot_panels_b_to_g() -> None:
    fig, axes = plt.subplots(2, 3, figsize=(14.5, 8.7), facecolor=FIG_BG)
    panel_axes = [
        ("B", axes[0, 0]),
        ("D", axes[0, 1]),
        ("F", axes[0, 2]),
        ("C", axes[1, 0]),
        ("E", axes[1, 1]),
        ("G", axes[1, 2]),
    ]
    for letter, ax in panel_axes:
        draw_panel_bg(ax, letter)
    fig.suptitle("Gene score(GSE149512)", fontsize=18, fontweight="bold", y=0.99)
    fig.subplots_adjust(top=0.90, hspace=0.32, wspace=0.12)
    save(fig, "SupFig9BG_violin_grid")


def plot_panels_h_to_k() -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.8), facecolor=FIG_BG)
    draw_panel_h(ax, corr_text=_corr_text("H"))
    save(fig, "SupFig9H_scatter")

    fig, ax = plt.subplots(figsize=(5.8, 4.6), facecolor=FIG_BG)
    draw_panel_i(ax)
    save(fig, "SupFig9I_box")

    fig, ax = plt.subplots(figsize=(6.4, 4.8), facecolor=FIG_BG)
    draw_panel_j(ax, corr_text=_corr_text("J"))
    save(fig, "SupFig9J_scatter")

    fig, ax = plt.subplots(figsize=(5.8, 4.6), facecolor=FIG_BG)
    draw_panel_k(ax)
    save(fig, "SupFig9K_box")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), facecolor=FIG_BG)
    draw_panel_h(axes[0, 0], corr_text=_corr_text("H"))
    draw_panel_j(axes[0, 1], corr_text=_corr_text("J"))
    draw_panel_i(axes[1, 0])
    draw_panel_k(axes[1, 1])
    fig.subplots_adjust(hspace=0.32, wspace=0.20)
    save(fig, "SupFig9HIJK_combined")


def plot_supfig9_main() -> None:
    fig = plt.figure(figsize=(14.6, 20.2), facecolor=FIG_BG)
    gs = fig.add_gridspec(
        nrows=5,
        ncols=12,
        height_ratios=[3.5, 2.0, 2.0, 2.25, 2.25],
        hspace=0.42,
        wspace=0.35,
    )

    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0:4])
    ax_d = fig.add_subplot(gs[1, 4:8])
    ax_f = fig.add_subplot(gs[1, 8:12])
    ax_c = fig.add_subplot(gs[2, 0:4])
    ax_e = fig.add_subplot(gs[2, 4:8])
    ax_g = fig.add_subplot(gs[2, 8:12])
    ax_h = fig.add_subplot(gs[3, 0:5])
    ax_leg = fig.add_subplot(gs[3, 5:7])
    ax_j = fig.add_subplot(gs[3, 7:12])
    ax_i = fig.add_subplot(gs[4, 0:6])
    ax_k = fig.add_subplot(gs[4, 6:12])

    draw_panel_a(ax_a)
    draw_panel_bg(ax_b, "B")
    draw_panel_bg(ax_d, "D")
    draw_panel_bg(ax_f, "F")
    draw_panel_bg(ax_c, "C")
    draw_panel_bg(ax_e, "E")
    draw_panel_bg(ax_g, "G")
    draw_panel_h(ax_h, corr_text=_corr_text("H"))
    draw_panel_j(ax_j, corr_text=_corr_text("J"))
    draw_panel_i(ax_i)
    draw_panel_k(ax_k)

    ax_leg.set_facecolor(FIG_BG)
    ax_leg.axis("off")
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markerfacecolor=COL_HJK[g],
            markeredgecolor="#333333",
            markeredgewidth=0.8,
            markersize=9,
            label=g,
        )
        for g in GROUP_ORDER_HJK
    ]
    ax_leg.legend(handles=handles, title="Group", loc="center", frameon=True, fontsize=12, title_fontsize=14)

    fig.text(0.03, 0.992, "Supplementary Figure 9", ha="left", va="top", fontsize=44)
    mid_bg_x = (ax_b.get_position().x0 + ax_f.get_position().x1) / 2
    fig.text(mid_bg_x, ax_b.get_position().y1 + 0.006, "Gene score(GSE149512)", ha="center", va="bottom", fontsize=18, fontweight="bold")

    _panel_letter(fig, ax_a, "A")
    _panel_letter(fig, ax_b, "B")
    _panel_letter(fig, ax_c, "C")
    _panel_letter(fig, ax_d, "D")
    _panel_letter(fig, ax_e, "E")
    _panel_letter(fig, ax_f, "F")
    _panel_letter(fig, ax_g, "G")
    _panel_letter(fig, ax_h, "H")
    _panel_letter(fig, ax_i, "I")
    _panel_letter(fig, ax_j, "J")
    _panel_letter(fig, ax_k, "K")

    save(fig, "SupFig9_main")


def main() -> None:
    ensure_out()
    plot_panel_a()
    plot_panels_b_to_g()
    plot_panels_h_to_k()
    plot_supfig9_main()


if __name__ == "__main__":
    main()
