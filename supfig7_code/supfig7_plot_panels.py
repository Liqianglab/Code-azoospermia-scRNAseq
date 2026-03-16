#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
SRC = BASE / "supfig7_source_data"
OUT = BASE / "supfig7_code" / "supfig7_plots"


def ensure_out() -> None:
    OUT.mkdir(parents=True, exist_ok=True)


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def panel_ab(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "SupFig7AB_panel_numbers_wide.csv")
    long_rows: list[dict[str, object]] = []
    for _, r in df.iterrows():
        panel = str(r["panel"])
        long_rows.append({"panel": panel, "celltype": "Shared", "count": int(r["shared_count"])})
        for key in ["ECs_specific_count", "LCs_specific_count", "Lym_specific_count", "Myoids_specific_count", "STs_specific_count"]:
            long_rows.append({"panel": panel, "celltype": key.replace("_specific_count", ""), "count": int(r[key])})
    long = pd.DataFrame(long_rows)

    x_labels = ["A", "B"]
    celltypes = ["Shared", "ECs", "LCs", "Lym", "Myoids", "STs"]
    x = np.arange(len(x_labels))
    width = 0.12

    for i, ct in enumerate(celltypes):
        vals = [
            int(long[(long["panel"] == p) & (long["celltype"] == ct)]["count"].sum())
            for p in x_labels
        ]
        ax.bar(x + (i - (len(celltypes) - 1) / 2) * width, vals, width=width, label=ct)

    ax.set_xticks(x)
    ax.set_xticklabels(["A Up in OA", "B Down in OA"], fontsize=8)
    ax.set_ylabel("Gene count")
    ax.set_title("A-B OA overlap counts", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=7, ncol=3)


def _plot_go(ax: plt.Axes, csv_name: str, title: str) -> None:
    df = pd.read_csv(SRC / csv_name)
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
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")


def panel_e(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "SupFig7E_cilia13_dotplot_ctrl_oa_long.csv")

    groups = ["Ctrl", "OA"]
    genes = (
        df.groupby("gene", as_index=False)["normalised_expression_value"]
        .mean()
        .sort_values("normalised_expression_value", ascending=False)["gene"]
        .tolist()
    )
    if len(genes) > 13:
        genes = genes[:13]
    sub = df[df["gene"].isin(genes)].copy()

    x_map = {g: i for i, g in enumerate(groups)}
    y_map = {g: i for i, g in enumerate(genes)}

    x = sub["gname"].map(x_map)
    y = sub["gene"].map(y_map)
    size = pd.to_numeric(sub["percentage"], errors="coerce").fillna(0).clip(lower=0) * 800
    color = pd.to_numeric(sub["normalised_expression_value"], errors="coerce").fillna(0)

    sc = ax.scatter(x, y, s=size, c=color, cmap="viridis", edgecolors="black", linewidths=0.2)
    ax.set_xticks(np.arange(len(groups)))
    ax.set_xticklabels(groups)
    ax.set_yticks(np.arange(len(genes)))
    ax.set_yticklabels(genes, fontsize=7)
    ax.set_title("E Cilia13 dotplot", loc="left", fontsize=11, fontweight="bold")
    cbar = plt.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Normalised expression", fontsize=8)


def panel_fg_note(ax: plt.Axes) -> None:
    f_note = pd.read_csv(SRC / "SupFig7F_source_note.csv").iloc[0].to_dict()
    g_note = pd.read_csv(SRC / "SupFig7G_source_note.csv").iloc[0].to_dict()
    txt = (
        "F note: " + str(f_note.get("note", "")) + "\n\n" +
        "G note: " + str(g_note.get("note", ""))
    )
    ax.axis("off")
    ax.text(0.02, 0.95, txt, ha="left", va="top", fontsize=9)


def main() -> None:
    ensure_out()

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    panel_ab(axes[0, 0])
    _plot_go(axes[0, 1], "SupFig7C_GO_selected_combined.csv", "C GO selected terms")
    _plot_go(axes[0, 2], "SupFig7D_GO_selected_combined.csv", "D GO selected terms")
    panel_e(axes[1, 0])
    panel_fg_note(axes[1, 1])
    axes[1, 2].axis("off")

    for ax in axes.ravel():
        if ax.has_data():
            ax.grid(True, axis="y", linestyle="--", alpha=0.25)

    fig.suptitle("Supplementary Figure 7 (source-data reconstruction)", fontsize=18, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    save(fig, "SupFig7_main")


if __name__ == "__main__":
    main()
