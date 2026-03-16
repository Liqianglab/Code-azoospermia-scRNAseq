#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]  # .../source data
IN_FILE = BASE / "supfig6_source_data" / "SupFig6A_scmeta_heatmap_top30.csv"
OUT_DIR = Path(__file__).resolve().parent / "supfig6_plots"


def main() -> None:
    if not IN_FILE.exists():
        raise FileNotFoundError(f"Input not found: {IN_FILE}")

    df = pd.read_csv(IN_FILE)
    if "pathway" not in df.columns:
        raise ValueError("SupFig6A_scmeta_heatmap_top30.csv must contain 'pathway' column")

    mat = df.set_index("pathway")
    values = mat.to_numpy(dtype=float)

    vmin = float(np.nanmin(values)) if np.isfinite(values).any() else -2.0
    vmax = float(np.nanmax(values)) if np.isfinite(values).any() else 2.0
    lim = max(abs(vmin), abs(vmax), 2.0)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fig_w = max(8.0, 0.75 * mat.shape[1] + 3.5)
    fig_h = max(10.0, 0.24 * mat.shape[0] + 2.2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=150)

    im = ax.imshow(values, aspect="auto", cmap="YlGnBu_r", vmin=-lim, vmax=lim)

    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels(mat.index.tolist(), fontsize=9)
    ax.set_xlabel("celltype")
    ax.set_ylabel("pathway")
    ax.set_title("SupFig6A: scMetabolism heatmap (top pathways)")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Expression")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "SupFig6A_heatmap_top30.png", bbox_inches="tight")
    fig.savefig(OUT_DIR / "SupFig6A_heatmap_top30.svg", bbox_inches="tight")


if __name__ == "__main__":
    main()
