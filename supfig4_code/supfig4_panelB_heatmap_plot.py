#!/usr/bin/env python3
import csv
import html
import math
from pathlib import Path


def parse_float(text):
    value = str(text).strip()
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def row_zscore(values):
    numeric = [value for value in values if value is not None]
    if len(numeric) < 2:
        return [0.0 if value is not None else None for value in values]
    mean_val = sum(numeric) / len(numeric)
    variance = sum((value - mean_val) ** 2 for value in numeric) / len(numeric)
    std = math.sqrt(variance)
    if std == 0:
        return [0.0 if value is not None else None for value in values]
    return [None if value is None else (value - mean_val) / std for value in values]


def blend(color1, color2, ratio):
    ratio = max(0.0, min(1.0, ratio))
    return tuple(round(color1[idx] + (color2[idx] - color1[idx]) * ratio) for idx in range(3))


def color_from_zscore(zscore, lim=2.0):
    if zscore is None:
        return "#E6E6E6"
    blue = (49, 54, 149)
    white = (247, 247, 247)
    red = (165, 0, 38)
    value = max(-lim, min(lim, zscore))
    if value >= 0:
        rgb = blend(white, red, value / lim)
    else:
        rgb = blend(white, blue, abs(value) / lim)
    return "#{:02X}{:02X}{:02X}".format(*rgb)


def read_matrix(path):
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            gene = raw.get("target_gene", "").strip()
            values = [
                parse_float(raw.get("Stage_a", "")),
                parse_float(raw.get("Stage_b", "")),
                parse_float(raw.get("Stage_c", "")),
            ]
            rows.append(
                {
                    "gene": gene,
                    "available": raw.get("available_in_aucell", "0").strip() == "1",
                    "values": values,
                    "zscores": row_zscore(values),
                }
            )
    return rows


def render_svg(rows, out_path):
    columns = ["Stage_a", "Stage_b", "Stage_c"]
    left = 220
    top = 120
    cell_w = 120
    cell_h = 38
    right = 220
    bottom = 110
    width = left + len(columns) * cell_w + right
    height = top + len(rows) * cell_h + bottom

    plot_w = len(columns) * cell_w
    plot_h = len(rows) * cell_h
    legend_x = left + plot_w + 80
    legend_y = top + 6
    grad_h = 220
    grad_w = 20

    lines = [f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}'>"]
    lines.append(
        "<defs><linearGradient id='grad' x1='0%' y1='100%' x2='0%' y2='0%'>"
        "<stop offset='0%' stop-color='#313695'/>"
        "<stop offset='50%' stop-color='#F7F7F7'/>"
        "<stop offset='100%' stop-color='#A50026'/>"
        "</linearGradient></defs>"
    )
    lines.append(
        "<text x='20' y='42' font-family='Arial' font-size='30' font-weight='bold'>Supplementary Figure 4B (proxy)</text>"
    )
    lines.append(
        "<text x='20' y='75' font-family='Arial' font-size='18'>Regulon activity by ST stage (row z-score, available rows from AUC matrix)</text>"
    )

    for idx, col in enumerate(columns):
        x = left + (idx + 0.5) * cell_w
        lines.append(
            f"<text x='{x:.2f}' y='{top - 14}' text-anchor='middle' font-family='Arial' font-size='15'>{html.escape(col)}</text>"
        )

    for ridx, row in enumerate(rows):
        y = top + (ridx + 0.5) * cell_h + 5
        label = row["gene"] + ("" if row["available"] else " (missing)")
        lines.append(
            f"<text x='{left - 12}' y='{y:.2f}' text-anchor='end' font-family='Arial' font-size='14'>{html.escape(label)}</text>"
        )

        for cidx in range(len(columns)):
            x = left + cidx * cell_w
            z = row["zscores"][cidx]
            fill = color_from_zscore(z)
            lines.append(
                f"<rect x='{x:.2f}' y='{top + ridx * cell_h:.2f}' width='{cell_w}' height='{cell_h}' fill='{fill}' stroke='#D0D0D0' stroke-width='1'/>"
            )

    lines.append(
        f"<rect x='{left}' y='{top}' width='{plot_w}' height='{plot_h}' fill='none' stroke='#333333' stroke-width='1.2'/>"
    )

    lines.append(
        f"<text x='{legend_x - 4}' y='{legend_y - 10}' font-family='Arial' font-size='14' font-weight='bold'>Matrix (z-score)</text>"
    )
    lines.append(
        f"<rect x='{legend_x}' y='{legend_y}' width='{grad_w}' height='{grad_h}' fill='url(#grad)' stroke='#333333' stroke-width='0.8'/>"
    )
    lines.append(
        f"<text x='{legend_x + grad_w + 8}' y='{legend_y + 6}' font-family='Arial' font-size='12'>+2</text>"
    )
    lines.append(
        f"<text x='{legend_x + grad_w + 8}' y='{legend_y + grad_h / 2 + 4}' font-family='Arial' font-size='12'>0</text>"
    )
    lines.append(
        f"<text x='{legend_x + grad_w + 8}' y='{legend_y + grad_h + 4}' font-family='Arial' font-size='12'>-2</text>"
    )
    lines.append(
        f"<rect x='{legend_x}' y='{legend_y + grad_h + 30}' width='18' height='18' fill='#E6E6E6' stroke='#BEBEBE'/>"
    )
    lines.append(
        f"<text x='{legend_x + 28}' y='{legend_y + grad_h + 44}' font-family='Arial' font-size='12'>missing regulon</text>"
    )

    lines.append("</svg>")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def main():
    script_dir = Path(__file__).resolve().parent
    data_path = script_dir.parent / "supfig4_source_data" / "SupFig4B_regulon_stage_matrix_target.csv"
    out_dir = script_dir / "supfig4_plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "SupFig4B_regulon_stage_heatmap.svg"

    rows = read_matrix(data_path)
    render_svg(rows, out_path)
    print(f"Wrote SVG plot to: {out_path}")


if __name__ == "__main__":
    main()
