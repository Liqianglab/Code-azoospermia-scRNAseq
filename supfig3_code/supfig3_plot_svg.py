#!/usr/bin/env python3
import csv
import html
from pathlib import Path


def read_panel_matrix(path):
    rows = []
    comparisons = []
    pathways = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row["NES"] = float(row["NES"])
            row["neglog10_fdr"] = float(row["neglog10_fdr"])
            row["is_missing"] = int(row["is_missing"])
            rows.append(row)
            if row["comparison"] not in comparisons:
                comparisons.append(row["comparison"])
            if row["pathway_readable"] not in pathways:
                pathways.append(row["pathway_readable"])
    return {"rows": rows, "comparisons": comparisons, "pathways": pathways}


def blend(c1, c2, t):
    t = max(0.0, min(1.0, t))
    return tuple(round(c1[i] + (c2[i] - c1[i]) * t) for i in range(3))


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


def radius_from_size(value, max_value):
    if max_value <= 0:
        return 2.0
    return 2.0 + 8.0 * (value / max_value)


def panel_geometry(panel):
    left = 360
    top = 90
    cell_w = 130
    cell_h = 28
    right = 230
    bottom = 130

    width = left + cell_w * len(panel["comparisons"]) + right
    height = top + cell_h * len(panel["pathways"]) + bottom
    return {
        "left": left,
        "top": top,
        "cell_w": cell_w,
        "cell_h": cell_h,
        "right": right,
        "bottom": bottom,
        "width": width,
        "height": height,
    }


def draw_panel(svg, panel, geom, x0, y0, title, panel_label, gradient_id):
    comp_index = {c: i for i, c in enumerate(panel["comparisons"])}
    path_index = {p: i for i, p in enumerate(panel["pathways"])}
    max_abs_nes = max(abs(r["NES"]) for r in panel["rows"] if r["is_missing"] == 0)
    max_size = max(r["neglog10_fdr"] for r in panel["rows"])

    left = geom["left"]
    top = geom["top"]
    cell_w = geom["cell_w"]
    cell_h = geom["cell_h"]
    plot_w = cell_w * len(panel["comparisons"])
    plot_h = cell_h * len(panel["pathways"])

    svg.append(
        f"<defs><linearGradient id='{gradient_id}' x1='0%' y1='100%' x2='0%' y2='0%'>"
        f"<stop offset='0%' stop-color='#3B4CC0'/>"
        f"<stop offset='50%' stop-color='#F7F7F7'/>"
        f"<stop offset='100%' stop-color='#B40426'/></linearGradient></defs>"
    )

    svg.append(
        f"<text x='{x0 + 24}' y='{y0 + 42}' font-size='34' font-family='Arial' font-weight='bold'>{html.escape(panel_label)}</text>"
    )
    svg.append(
        f"<text x='{x0 + left}' y='{y0 + 46}' font-size='24' font-family='Arial' font-weight='bold'>{html.escape(title)}</text>"
    )

    for i in range(len(panel["pathways"]) + 1):
        y_line = y0 + top + i * cell_h
        svg.append(
            f"<line x1='{x0 + left}' y1='{y_line:.2f}' x2='{x0 + left + plot_w}' y2='{y_line:.2f}' "
            f"stroke='#E6E6E6' stroke-width='1'/>"
        )

    svg.append(
        f"<rect x='{x0 + left}' y='{y0 + top}' width='{plot_w}' height='{plot_h}' fill='none' stroke='#333333' stroke-width='1.2'/>"
    )

    for pathway in panel["pathways"]:
        i = path_index[pathway]
        y = y0 + top + (i + 0.5) * cell_h + 5
        svg.append(
            f"<text x='{x0 + left - 10}' y='{y:.2f}' text-anchor='end' font-size='15' font-family='Arial'>{html.escape(pathway)}</text>"
        )

    xlab_y = y0 + top + plot_h + 36
    for comp in panel["comparisons"]:
        j = comp_index[comp]
        x = x0 + left + (j + 0.5) * cell_w
        svg.append(
            f"<text transform='translate({x:.2f},{xlab_y:.2f}) rotate(35)' text-anchor='start' "
            f"font-size='15' font-family='Arial'>{html.escape(comp)}</text>"
        )

    for row in panel["rows"]:
        i = path_index[row["pathway_readable"]]
        j = comp_index[row["comparison"]]
        x = x0 + left + (j + 0.5) * cell_w
        y = y0 + top + (i + 0.5) * cell_h
        r = radius_from_size(row["neglog10_fdr"], max_size)
        if row["is_missing"] == 1:
            color = "#E6E6E6"
            opacity = "0.80"
        else:
            color = color_from_nes(row["NES"], max_abs_nes)
            opacity = "0.95"
        svg.append(
            f"<circle cx='{x:.2f}' cy='{y:.2f}' r='{r:.2f}' fill='{color}' fill-opacity='{opacity}' stroke='none'/>"
        )

    legend_x = x0 + left + plot_w + 90
    legend_y = y0 + top + 12
    grad_w = 18
    grad_h = 190
    svg.append(
        f"<text x='{legend_x - 3}' y='{legend_y - 10}' font-size='15' font-family='Arial' font-weight='bold'>NES</text>"
    )
    svg.append(
        f"<rect x='{legend_x}' y='{legend_y}' width='{grad_w}' height='{grad_h}' fill='url(#{gradient_id})' stroke='#333333' stroke-width='0.8'/>"
    )
    svg.append(
        f"<text x='{legend_x + grad_w + 8}' y='{legend_y + 5}' font-size='12' font-family='Arial'>{max_abs_nes:.2f}</text>"
    )
    svg.append(
        f"<text x='{legend_x + grad_w + 8}' y='{legend_y + grad_h / 2 + 4}' font-size='12' font-family='Arial'>0</text>"
    )
    svg.append(
        f"<text x='{legend_x + grad_w + 8}' y='{legend_y + grad_h + 4}' font-size='12' font-family='Arial'>-{max_abs_nes:.2f}</text>"
    )

    size_title_y = legend_y + grad_h + 44
    svg.append(
        f"<text x='{legend_x - 3}' y='{size_title_y}' font-size='15' font-family='Arial' font-weight='bold'>-log10(FDR)</text>"
    )
    size_values = [1, 2, 3, 4, 6]
    y_cursor = size_title_y + 20
    for value in size_values:
        if value > max_size and value != 1:
            continue
        radius = radius_from_size(value, max_size)
        cy = y_cursor + radius
        svg.append(f"<circle cx='{legend_x + 12}' cy='{cy:.2f}' r='{radius:.2f}' fill='#1A1A1A'/>")
        svg.append(
            f"<text x='{legend_x + 34}' y='{cy + 4:.2f}' font-size='12' font-family='Arial'>{value}</text>"
        )
        y_cursor = cy + radius + 8


def compose_svg(panelA, panelB):
    geomA = panel_geometry(panelA)
    geomB = panel_geometry(panelB)
    canvas_w = max(geomA["width"], geomB["width"]) + 30
    gap = 56
    canvas_h = geomA["height"] + geomB["height"] + gap + 20

    xA = (canvas_w - geomA["width"]) / 2
    yA = 0
    xB = (canvas_w - geomB["width"]) / 2
    yB = geomA["height"] + gap

    svg = [f"<svg xmlns='http://www.w3.org/2000/svg' width='{canvas_w:.0f}' height='{canvas_h:.0f}'>"]
    draw_panel(
        svg,
        panelA,
        geomA,
        xA,
        yA,
        "Late primary spermatocytes: Hallmark enrichment (disease vs Ctrl)",
        "A",
        "gradA",
    )
    draw_panel(
        svg,
        panelB,
        geomB,
        xB,
        yB,
        "Round spermatids: Hallmark enrichment (disease vs Ctrl)",
        "B",
        "gradB",
    )
    svg.append("</svg>")
    return "\n".join(svg)


def compose_single(panel, title, panel_label, gradient_id):
    geom = panel_geometry(panel)
    canvas_w = geom["width"] + 20
    canvas_h = geom["height"] + 10
    x0 = (canvas_w - geom["width"]) / 2
    y0 = 0
    svg = [f"<svg xmlns='http://www.w3.org/2000/svg' width='{canvas_w:.0f}' height='{canvas_h:.0f}'>"]
    draw_panel(svg, panel, geom, x0, y0, title, panel_label, gradient_id)
    svg.append("</svg>")
    return "\n".join(svg)


def main():
    script_dir = Path(__file__).resolve().parent
    data_dir = script_dir.parent / "supfig3_source_data"
    out_dir = script_dir / "supfig3_plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    panelA = read_panel_matrix(data_dir / "SupFig3A_plot_matrix.csv")
    panelB = read_panel_matrix(data_dir / "SupFig3B_plot_matrix.csv")

    singleA = compose_single(
        panelA,
        "Late primary spermatocytes: Hallmark enrichment (disease vs Ctrl)",
        "A",
        "gradA_single",
    )
    singleB = compose_single(
        panelB,
        "Round spermatids: Hallmark enrichment (disease vs Ctrl)",
        "B",
        "gradB_single",
    )
    combined = compose_svg(panelA, panelB)

    (out_dir / "SupFig3A.svg").write_text(singleA, encoding="utf-8")
    (out_dir / "SupFig3B.svg").write_text(singleB, encoding="utf-8")
    (out_dir / "SupFig3_combined.svg").write_text(combined, encoding="utf-8")
    print(f"Wrote SVG plots to: {out_dir}")


if __name__ == "__main__":
    main()
