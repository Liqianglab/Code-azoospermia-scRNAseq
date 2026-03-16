from __future__ import annotations

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
SRC = BASE / "fig6_source_data"


def plot_dotplot(df: pd.DataFrame, out: Path, title: str) -> None:
    groups = list(df["gname"].drop_duplicates())
    genes = list(df["Gene name"].drop_duplicates())

    x_map = {g: i for i, g in enumerate(groups)}
    y_map = {g: i for i, g in enumerate(genes)}

    x = df["gname"].map(x_map)
    y = df["Gene name"].map(y_map)
    sizes = (df["Percentage"].astype(float) * 900).clip(lower=10)
    colors = df["Normalised expression value"].astype(float)

    fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=200)
    sc = ax.scatter(x, y, s=sizes, c=colors, cmap="Reds", edgecolors="none")

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=45, ha="right")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes)
    ax.invert_yaxis()
    ax.set_title(title)

    cbar = plt.colorbar(sc, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Average expression")

    plt.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)


def main() -> None:
    c_df = pd.read_csv(SRC / "Fig6C_dotplot_LCa.csv")
    d_df = pd.read_csv(SRC / "Fig6D_dotplot_ST_no_OA.csv")

    plot_dotplot(c_df, SRC / "Fig6C_dotplot_LCa_replot.pdf", "Fig6C LC_a")
    plot_dotplot(d_df, SRC / "Fig6D_dotplot_ST_replot.pdf", "Fig6D ST")


if __name__ == "__main__":
    main()
