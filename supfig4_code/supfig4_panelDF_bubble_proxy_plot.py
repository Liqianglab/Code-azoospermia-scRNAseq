#!/usr/bin/env python3
import csv
import html
from pathlib import Path

STALE_PLOT_FILES = [
    "SupFig4D_ST1vsST3_proxy_bubble.svg",
    "SupFig4F_ST3vsST1_proxy_bubble.svg",
]


def parse_float(text):
    value = str(text).strip()
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def blend(color1, color2, ratio):
    ratio = max(0.0, min(1.0, ratio))
    return tuple(round(color1[idx] + (color2[idx] - color1[idx]) * ratio) for idx in range(3))


def color_from_nes(nes, lim):
    blue = (59, 76, 192)
    white = (247, 247, 247)
    red = (180, 4, 38)
    if lim <= 0:
        rgb = white
    elif nes >= 0:
        rgb = blend(white, red, nes / lim)
    else:
        rgb = blend(white, blue, abs(nes) / lim)
    return "#{:02X}{:02X}{:02X}".format(*rgb)


def load_rows(path):
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            nes = parse_float(raw.get("NES", ""))
            neglog = parse_float(raw.get("neglog10_fdr", ""))
            if nes is None or neglog is None:
                continue
            rows.append(
                {
                    "pathway": (
                        raw.get("pathway_readable", "").strip()
                        or raw.get("Description", "").replace("KEGG_", "").replace("REACTOME_", "").replace("_", " ")
                    ),
                    "nes": nes,
                    "neglog10_fdr": neglog,
                    "padj": raw.get("p.adjust", ""),
                }
            )
    return rows


def render_panel(rows, title, panel_label, out_path):
    rows = list(rows)
    left = 500
    top = 90
    row_h = 30
    plot_w = 420
    right = 230
    bottom = 90
    width = left + plot_w + right
    height = top + row_h * len(rows) + bottom

    max_abs_nes = max(abs(row["nes"]) for row in rows) if rows else 1.0
    max_size = max(row["neglog10_fdr"] for row in rows) if rows else 1.0
    x_min = -max_abs_nes * 1.05
    x_max = max_abs_nes * 1.05

    def x_scale(value):
        if x_max == x_min:
            return left + plot_w / 2
        return left + (value - x_min) / (x_max - x_min) * plot_w

    def radius(value):
        if max_size <= 0:
            return 3.0
        return 3.0 + 9.0 * (value / max_size)

    lines = [f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}'>"]
    lines.append(
        f"<text x='20' y='40' font-family='Arial' font-size='31' font-weight='bold'>{html.escape(panel_label)}</text>"
    )
    lines.append(
        f"<text x='{left}' y='42' font-family='Arial' font-size='22' font-weight='bold'>{html.escape(title)}</text>"
    )
    lines.append(
        "<text x='20' y='70' font-family='Arial' font-size='16'>Proxy plot from source-data GSEA table</text>"
    )

    zero_x = x_scale(0.0)
    lines.append(
        f"<line x1='{zero_x:.2f}' y1='{top - 8}' x2='{zero_x:.2f}' y2='{top + row_h * len(rows)}' stroke='#777777' stroke-width='1.2' stroke-dasharray='4,4'/>"
    )
    lines.append(
        f"<rect x='{left}' y='{top}' width='{plot_w}' height='{row_h * len(rows)}' fill='none' stroke='#333333' stroke-width='1.2'/>"
    )

    tick_values = [round(x_min, 2), 0.0, round(x_max, 2)]
    for tick in tick_values:
        x = x_scale(tick)
        lines.append(
            f"<line x1='{x:.2f}' y1='{top + row_h * len(rows)}' x2='{x:.2f}' y2='{top + row_h * len(rows) + 6}' stroke='#333333' stroke-width='1'/>"
        )
        lines.append(
            f"<text x='{x:.2f}' y='{top + row_h * len(rows) + 22}' text-anchor='middle' font-family='Arial' font-size='12'>{tick:.2f}</text>"
        )
    lines.append(
        f"<text x='{left + plot_w / 2:.2f}' y='{top + row_h * len(rows) + 48}' text-anchor='middle' font-family='Arial' font-size='13'>NES</text>"
    )

    for idx, row in enumerate(rows):
        y = top + (idx + 0.5) * row_h
        lines.append(
            f"<line x1='{left}' y1='{top + (idx + 1) * row_h:.2f}' x2='{left + plot_w}' y2='{top + (idx + 1) * row_h:.2f}' stroke='#F0F0F0' stroke-width='1'/>"
        )
        lines.append(
            f"<text x='{left - 12}' y='{y + 4:.2f}' text-anchor='end' font-family='Arial' font-size='12'>{html.escape(row['pathway'])}</text>"
        )
        cx = x_scale(row["nes"])
        fill = color_from_nes(row["nes"], max_abs_nes)
        r = radius(row["neglog10_fdr"])
        lines.append(
            f"<circle cx='{cx:.2f}' cy='{y:.2f}' r='{r:.2f}' fill='{fill}' fill-opacity='0.95' stroke='none'/>"
        )

    legend_x = left + plot_w + 90
    legend_y = top + 8
    lines.append(
        f"<text x='{legend_x - 4}' y='{legend_y - 10}' font-family='Arial' font-size='14' font-weight='bold'>-log10(FDR)</text>"
    )
    for i, size in enumerate([1, 2, 3, 4]):
        if size > max_size and size > 1:
            continue
        cy = legend_y + 20 + i * 34
        r = radius(size)
        lines.append(f"<circle cx='{legend_x + 10}' cy='{cy:.2f}' r='{r:.2f}' fill='#1A1A1A'/>")
        lines.append(f"<text x='{legend_x + 28}' y='{cy + 4:.2f}' font-family='Arial' font-size='12'>{size}</text>")

    lines.append(
        f"<text x='{legend_x - 4}' y='{legend_y + 170}' font-family='Arial' font-size='14' font-weight='bold'>Color</text>"
    )
    lines.append(
        f"<text x='{legend_x - 4}' y='{legend_y + 188}' font-family='Arial' font-size='12'>NES: blue (-) to red (+)</text>"
    )
    lines.append("</svg>")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def main():
    script_dir = Path(__file__).resolve().parent
    data_dir = script_dir.parent / "supfig4_source_data"
    out_dir = script_dir / "supfig4_plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in STALE_PLOT_FILES:
        stale_plot = out_dir / name
        if stale_plot.exists():
            stale_plot.unlink()

    panel_d_rows = load_rows(data_dir / "SupFig4D_ST1_disease_vs_ctrl_top10.csv")
    panel_f_rows = load_rows(data_dir / "SupFig4F_ST3_disease_vs_ctrl_top10.csv")

    render_panel(
        panel_d_rows,
        "ST1 disease vs ctrl: top enriched pathways",
        "SupFig4D (proxy)",
        out_dir / "SupFig4D_ST1_disease_vs_ctrl_proxy_bubble.svg",
    )
    render_panel(
        panel_f_rows,
        "ST3 disease vs ctrl: top enriched pathways",
        "SupFig4F (proxy)",
        out_dir / "SupFig4F_ST3_disease_vs_ctrl_proxy_bubble.svg",
    )
    print(f"Wrote SVG plots to: {out_dir}")


if __name__ == "__main__":
    main()
