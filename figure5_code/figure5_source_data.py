#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

NS = {
    "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}
CELL_REF_RE = re.compile(r"^([A-Z]+)([0-9]+)$")

GROUP_ORDER = ["Ctrl", "AZFc_Del", "iNOA_B", "iNOA_S", "KS"]
GROUP_ORDER_PANEL_E = ["Ctrl", "AZFc_Del", "iNOA_S", "iNOA_B", "KS"]
GROUP_ORDER_PANEL_G = ["Ctrl", "KS", "AZFc_Del", "iNOA_B", "iNOA_S"]
GROUP_ORDER_HORMONE = ["Control", "OA", "NOA"]
STAGE_ORDER = ["Stage_a", "Stage_b", "Stage_c"]
STAGE_MAP = {"LC1": "Stage_a", "LC2": "Stage_b", "LC3": "Stage_c"}

PANEL_C_GENE_COLUMNS = {
    "NOTCH2": "NOTCH2 normalised expression value",
    "PDGFRB": "PDGFRB normalised expression value",
    "CYP17A1": "CYP17A1 normalised expression value",
    "INSL3": "INSL3 normalised expression value",
}

PANEL_G_GENE_ORDER = [
    "CYP51A1",
    "EBP",
    "ABCA5",
    "LIMA1",
    "SCARB1",
    "DHCR7",
    "TM7SF2",
    "DHCR24",
    "FDX1",
    "CYB5R3",
    "CYP11A1",
    "STAR",
    "CYP17A1",
]

PANEL_I_METRICS = [
    "Testosterone",
    "Luteinizing Hormone",
    "Thyroid-Stimulating Hormone",
]


def to_float(value):
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.upper() == "NA":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def quantile(values, q: float):
    if not values:
        return ""
    sorted_vals = sorted(values)
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = (len(sorted_vals) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return sorted_vals[lo]
    frac = pos - lo
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * frac


def normalize_group(token: str) -> str:
    t = str(token).strip()
    upper = t.upper()
    if upper in {"CTRL", "CONTROL"}:
        return "Ctrl"
    if upper == "AZFC" or "AZFC" in upper:
        return "AZFc_Del"
    if upper == "INOA_B":
        return "iNOA_B"
    if upper in {"INOA_S", "INOS_S"}:
        return "iNOA_S"
    if upper == "KS":
        return "KS"
    if upper == "OA":
        return "OA"
    if upper == "NOA":
        return "NOA"
    return t


def infer_group_from_cell_id(cell_id: str) -> str:
    token = str(cell_id).strip()
    if token.startswith("AZFc_Del"):
        return "AZFc_Del"
    if token.startswith("iNOA_B"):
        return "iNOA_B"
    if token.startswith("iNOA_S"):
        return "iNOA_S"
    if token.startswith("Ctrl"):
        return "Ctrl"
    if token.startswith("KS"):
        return "KS"
    if token.startswith("OA"):
        return "OA"
    match = re.match(r"^(Ctrl|OA|AZFc_Del|iNOA_B|iNOA_S|KS)", token)
    if match:
        return normalize_group(match.group(1))
    return ""


def read_csv_rows(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def column_index(col_letters: str) -> int:
    idx = 0
    for char in col_letters:
        idx = idx * 26 + (ord(char) - ord("A") + 1)
    return idx - 1


def read_xlsx_sheet1_rows(path: Path):
    with zipfile.ZipFile(path) as xlsx:
        shared_strings = []
        if "xl/sharedStrings.xml" in xlsx.namelist():
            shared_root = ET.fromstring(xlsx.read("xl/sharedStrings.xml"))
            for si in shared_root.findall("a:si", NS):
                shared_strings.append(
                    "".join((text_node.text or "") for text_node in si.findall(".//a:t", NS))
                )

        sheet_xml = xlsx.read("xl/worksheets/sheet1.xml")
        sheet_root = ET.fromstring(sheet_xml)

        rows = []
        for row_elem in sheet_root.findall(".//a:sheetData/a:row", NS):
            row_map = {}
            max_col = -1
            for cell_elem in row_elem.findall("a:c", NS):
                ref = cell_elem.get("r", "")
                ref_match = CELL_REF_RE.match(ref)
                col_idx = column_index(ref_match.group(1)) if ref_match else len(row_map)

                cell_type = cell_elem.get("t")
                text_value = ""
                value_elem = cell_elem.find("a:v", NS)

                if cell_type == "s" and value_elem is not None:
                    index = int(value_elem.text or "0")
                    if 0 <= index < len(shared_strings):
                        text_value = shared_strings[index]
                elif cell_type == "inlineStr":
                    inline_text = cell_elem.find("a:is/a:t", NS)
                    if inline_text is not None and inline_text.text is not None:
                        text_value = inline_text.text
                elif value_elem is not None and value_elem.text is not None:
                    text_value = value_elem.text

                row_map[col_idx] = text_value
                if col_idx > max_col:
                    max_col = col_idx

            if max_col >= 0:
                row_values = [row_map.get(idx, "") for idx in range(max_col + 1)]
            else:
                row_values = []
            rows.append(row_values)

        return rows


def sort_group(value: str):
    return GROUP_ORDER.index(value) if value in GROUP_ORDER else 999


def sort_stage(value: str):
    return STAGE_ORDER.index(value) if value in STAGE_ORDER else 999


def build_panel_a_c_d_e(expression_rows, out_dir: Path):
    lc_rows = []
    for row in expression_rows:
        cell = str(row.get("cell", "")).strip()
        group = infer_group_from_cell_id(cell)
        cluster = str(row.get("Major cell types", "")).strip()
        stage = STAGE_MAP.get(cluster, "")

        lc_rows.append(
            {
                "cell": cell,
                "group": group,
                "cluster": cluster,
                "stage": stage,
                "UMAP1": to_float(row.get("UMAP1")),
                "UMAP2": to_float(row.get("UMAP2")),
                "NOTCH2": to_float(row.get(PANEL_C_GENE_COLUMNS["NOTCH2"])),
                "PDGFRB": to_float(row.get(PANEL_C_GENE_COLUMNS["PDGFRB"])),
                "CYP17A1": to_float(row.get(PANEL_C_GENE_COLUMNS["CYP17A1"])),
                "INSL3": to_float(row.get(PANEL_C_GENE_COLUMNS["INSL3"])),
                "MAFB": to_float(row.get("MAFB normalised expression value")),
                "STAR": to_float(row.get("STAR normalised expression value")),
            }
        )

    lc_rows.sort(key=lambda r: (sort_group(r["group"]), r["cell"]))

    write_csv(
        out_dir / "Fig5A_LC_cells_umap_stage_group.csv",
        lc_rows,
        [
            "cell",
            "group",
            "cluster",
            "stage",
            "UMAP1",
            "UMAP2",
            "NOTCH2",
            "PDGFRB",
            "CYP17A1",
            "INSL3",
            "MAFB",
            "STAR",
        ],
    )

    # A counts
    by_group = Counter(row["group"] for row in lc_rows)
    by_stage = Counter(row["stage"] for row in lc_rows)
    by_group_stage = Counter((row["group"], row["stage"]) for row in lc_rows)

    rows_group = [
        {"group": g, "n_cells": by_group[g]}
        for g in sorted(by_group.keys(), key=sort_group)
        if g
    ]
    write_csv(out_dir / "Fig5A_LC_counts_by_group.csv", rows_group, ["group", "n_cells"])

    rows_stage = [
        {"stage": s, "n_cells": by_stage[s]}
        for s in sorted(by_stage.keys(), key=sort_stage)
        if s
    ]
    write_csv(out_dir / "Fig5A_LC_counts_by_stage.csv", rows_stage, ["stage", "n_cells"])

    rows_group_stage = []
    for group in sorted(by_group.keys(), key=sort_group):
        if not group:
            continue
        total = by_group[group]
        for stage in STAGE_ORDER:
            count = by_group_stage.get((group, stage), 0)
            rows_group_stage.append(
                {
                    "group": group,
                    "stage": stage,
                    "count": count,
                    "ratio": (count / total) if total else "",
                    "group_total": total,
                }
            )
    write_csv(
        out_dir / "Fig5A_LC_stage_counts_by_group.csv",
        rows_group_stage,
        ["group", "stage", "count", "ratio", "group_total"],
    )

    write_csv(
        out_dir / "Fig5A_stage_candidate_mapping.csv",
        [
            {"LC_cluster": "LC1", "Stage": "Stage_a"},
            {"LC_cluster": "LC2", "Stage": "Stage_b"},
            {"LC_cluster": "LC3", "Stage": "Stage_c"},
        ],
        ["LC_cluster", "Stage"],
    )

    write_csv(
        out_dir / "Fig5A_pseudotime_note.csv",
        [
            {
                "note": "Panel A right pseudotime trajectory is available as PDF panel in raw_support/Fig5A_monocle2_facet_group.pdf; numeric pseudotime coordinate table was not recovered in this folder.",
                "raw_panel": "raw_support/Fig5A_monocle2_facet_group.pdf",
            }
        ],
        ["note", "raw_panel"],
    )

    # Panel C long + summary
    marker_rows = []
    for row in lc_rows:
        for gene in ["NOTCH2", "PDGFRB", "CYP17A1", "INSL3"]:
            marker_rows.append(
                {
                    "cell": row["cell"],
                    "group": row["group"],
                    "cluster": row["cluster"],
                    "stage": row["stage"],
                    "gene": gene,
                    "expression": row[gene],
                }
            )

    marker_rows.sort(
        key=lambda r: (
            ["NOTCH2", "PDGFRB", "CYP17A1", "INSL3"].index(r["gene"]),
            sort_stage(r["stage"]),
            sort_group(r["group"]),
            r["cell"],
        )
    )

    write_csv(
        out_dir / "Fig5C_stage_marker_cells.csv",
        marker_rows,
        ["cell", "group", "cluster", "stage", "gene", "expression"],
    )

    summary = []
    grouped = defaultdict(list)
    for row in marker_rows:
        val = to_float(row["expression"])
        if val is not None:
            grouped[(row["gene"], row["stage"])].append(val)

    for gene in ["NOTCH2", "PDGFRB", "CYP17A1", "INSL3"]:
        for stage in STAGE_ORDER:
            vals = grouped.get((gene, stage), [])
            summary.append(
                {
                    "gene": gene,
                    "stage": stage,
                    "n_cells": len(vals),
                    "mean": (sum(vals) / len(vals)) if vals else "",
                    "median": quantile(vals, 0.5),
                    "q1": quantile(vals, 0.25),
                    "q3": quantile(vals, 0.75),
                }
            )

    write_csv(
        out_dir / "Fig5C_stage_marker_summary.csv",
        summary,
        ["gene", "stage", "n_cells", "mean", "median", "q1", "q3"],
    )

    # Panel D (reuses A group-stage counts)
    rows_d_counts = []
    rows_d_ratio = []
    for group in GROUP_ORDER:
        total = by_group.get(group, 0)
        if total == 0:
            continue
        for stage in STAGE_ORDER:
            count = by_group_stage.get((group, stage), 0)
            rows_d_counts.append(
                {
                    "group": group,
                    "stage": stage,
                    "count": count,
                    "group_total": total,
                }
            )
            rows_d_ratio.append(
                {
                    "group": group,
                    "stage": stage,
                    "ratio": (count / total) if total else "",
                }
            )

    write_csv(
        out_dir / "Fig5D_stage_counts_by_group.csv",
        rows_d_counts,
        ["group", "stage", "count", "group_total"],
    )
    write_csv(
        out_dir / "Fig5D_stage_ratio_by_group.csv",
        rows_d_ratio,
        ["group", "stage", "ratio"],
    )

    # Panel E (UMAP per-group highlight)
    rows_e = []
    for facet_group in GROUP_ORDER_PANEL_E:
        for row in lc_rows:
            rows_e.append(
                {
                    "facet_group": facet_group,
                    "cell": row["cell"],
                    "group": row["group"],
                    "stage": row["stage"],
                    "UMAP1": row["UMAP1"],
                    "UMAP2": row["UMAP2"],
                    "is_highlight": 1 if row["group"] == facet_group else 0,
                }
            )

    write_csv(
        out_dir / "Fig5E_group_highlight_umap_cells.csv",
        rows_e,
        ["facet_group", "cell", "group", "stage", "UMAP1", "UMAP2", "is_highlight"],
    )

    e_summary = []
    for facet_group in GROUP_ORDER_PANEL_E:
        total = sum(1 for row in lc_rows)
        hi = sum(1 for row in lc_rows if row["group"] == facet_group)
        e_summary.append(
            {
                "facet_group": facet_group,
                "highlight_cells": hi,
                "total_cells": total,
                "highlight_ratio": (hi / total) if total else "",
            }
        )
    write_csv(
        out_dir / "Fig5E_group_highlight_summary.csv",
        e_summary,
        ["facet_group", "highlight_cells", "total_cells", "highlight_ratio"],
    )

    write_csv(
        out_dir / "Fig5E_trajectory_note.csv",
        [
            {
                "note": "Panel E lower pseudotime trajectories were assembled from monocle panel PDF in raw_support/Fig5A_monocle2_facet_group.pdf; numeric trajectory coordinates by group were not recovered.",
                "raw_panel": "raw_support/Fig5A_monocle2_facet_group.pdf",
            }
        ],
        ["note", "raw_panel"],
    )


def build_panel_b(raw_dir: Path, out_dir: Path):
    rows = read_xlsx_sheet1_rows(raw_dir / "Fig5B_heatmap_export.xlsx")
    if not rows:
        return

    header = rows[0]
    cluster_cols = []
    for idx, name in enumerate(header):
        token = str(name).strip()
        if idx == 0:
            continue
        if token in STAGE_MAP:
            cluster_cols.append((idx, token))

    wide_rows = []
    long_rows = []
    gene_order = []

    for i, row in enumerate(rows[1:], start=1):
        if not row:
            continue
        gene = str(row[0]).strip() if len(row) > 0 else ""
        if gene == "":
            continue

        out_row = {"gene": gene}
        gene_order.append({"order": i, "gene": gene})
        for col_idx, cluster in cluster_cols:
            value = to_float(row[col_idx] if col_idx < len(row) else "")
            stage = STAGE_MAP.get(cluster, "")
            out_row[cluster] = value
            out_row[stage] = value
            long_rows.append(
                {
                    "gene": gene,
                    "cluster": cluster,
                    "stage": stage,
                    "scaled_expression": value,
                }
            )
        wide_rows.append(out_row)

    wide_cols = ["gene"] + [c for _, c in cluster_cols] + STAGE_ORDER
    write_csv(out_dir / "Fig5B_heatmap_matrix_wide.csv", wide_rows, wide_cols)
    write_csv(
        out_dir / "Fig5B_heatmap_matrix_long.csv",
        long_rows,
        ["gene", "cluster", "stage", "scaled_expression"],
    )
    write_csv(out_dir / "Fig5B_heatmap_gene_order.csv", gene_order, ["order", "gene"])


def parse_go_enrichment_sheet(rows):
    if not rows:
        return []
    header = rows[0]
    idx = {name: i for i, name in enumerate(header)}

    parsed = []
    for row in rows[1:]:
        if not row:
            continue

        def get_col(name, default=""):
            i = idx.get(name)
            if i is None or i >= len(row):
                return default
            return row[i]

        group_raw = str(get_col("group", "")).strip()
        group = normalize_group(group_raw)
        # This file has an extra final column (Up/down) that is not present in header row.
        regulation = str(row[11]).strip() if len(row) > 11 else ""
        count = to_float(get_col("Count", ""))
        p_adjust = to_float(get_col("p.adjust", ""))

        parsed.append(
            {
                "ONTOLOGY": get_col("ONTOLOGY", ""),
                "ID": get_col("ID", ""),
                "Description": get_col("Description", ""),
                "GeneRatio": get_col("GeneRatio", ""),
                "BgRatio": get_col("BgRatio", ""),
                "pvalue": to_float(get_col("pvalue", "")),
                "p.adjust": p_adjust,
                "qvalue": to_float(get_col("qvalue", "")),
                "geneID": get_col("geneID", ""),
                "Count": int(count) if count is not None else "",
                "group_raw": group_raw,
                "group": group,
                "regulation": regulation,
                "minus_log10_p_adjust": (-math.log10(p_adjust)) if (p_adjust and p_adjust > 0) else "",
            }
        )

    return parsed


def build_panel_f(raw_dir: Path, out_dir: Path):
    rows = read_xlsx_sheet1_rows(raw_dir / "Fig5F_LC_GOALL_enrichment_sig.xlsx")
    parsed = parse_go_enrichment_sheet(rows)

    down_rows = [r for r in parsed if str(r["regulation"]).lower() == "down"]
    down_rows.sort(key=lambda r: (sort_group(r["group"]), r["p.adjust"] if r["p.adjust"] != "" else 999))

    write_csv(
        out_dir / "Fig5F_GO_down_vs_ctrl_long.csv",
        down_rows,
        [
            "group_raw",
            "group",
            "regulation",
            "Description",
            "ID",
            "ONTOLOGY",
            "Count",
            "GeneRatio",
            "BgRatio",
            "pvalue",
            "p.adjust",
            "qvalue",
            "minus_log10_p_adjust",
            "geneID",
        ],
    )

    top_rows = []
    for group in ["AZFc_Del", "iNOA_B", "iNOA_S", "KS"]:
        candidates = [r for r in down_rows if r["group"] == group]
        candidates.sort(key=lambda r: (r["p.adjust"] if r["p.adjust"] != "" else 999, -r["Count"]))
        # Panel F in the assembled figure emphasizes metabolism/steroid-related terms.
        focused = [
            r
            for r in candidates
            if (
                "metabolic" in str(r["Description"]).lower()
                or "biosynthetic" in str(r["Description"]).lower()
            )
        ]
        selected = focused[:5]
        if len(selected) < 5:
            selected_ids = {id(row) for row in selected}
            for row in candidates:
                if id(row) in selected_ids:
                    continue
                selected.append(row)
                if len(selected) == 5:
                    break

        for rank, row in enumerate(selected, start=1):
            top_rows.append(
                {
                    "group_raw": row["group_raw"],
                    "group": row["group"],
                    "rank_in_group": rank,
                    "Description": row["Description"],
                    "Count": row["Count"],
                    "minus_log10_p_adjust": row["minus_log10_p_adjust"],
                    "p.adjust": row["p.adjust"],
                    "geneID": row["geneID"],
                }
            )

    write_csv(
        out_dir / "Fig5F_GO_down_vs_ctrl_top5_by_group.csv",
        top_rows,
        [
            "group_raw",
            "group",
            "rank_in_group",
            "Description",
            "Count",
            "minus_log10_p_adjust",
            "p.adjust",
            "geneID",
        ],
    )


def build_panel_g(raw_dir: Path, out_dir: Path):
    matrix_rows = read_csv_rows(raw_dir / "Fig5G_ALL-LC-matrix.csv")
    if not matrix_rows:
        return

    # First column name is empty in this file; csv.DictReader records it as ''.
    gene_col = ""
    if gene_col not in matrix_rows[0]:
        for key in matrix_rows[0].keys():
            if key.strip() == "":
                gene_col = key
                break

    matrix_by_gene = {}
    for row in matrix_rows:
        gene = str(row.get(gene_col, "")).strip()
        if gene:
            matrix_by_gene[gene] = row

    long_rows = []
    for gene in PANEL_G_GENE_ORDER:
        row = matrix_by_gene.get(gene)
        for group in GROUP_ORDER_PANEL_G:
            avg_expr = ""
            if row is not None:
                avg_expr = to_float(row.get(group, ""))
            long_rows.append(
                {
                    "gene": gene,
                    "group": group,
                    "average_expression": avg_expr,
                    "percent_expressed": "",
                }
            )

    write_csv(
        out_dir / "Fig5G_dotplot_average_expression_long.csv",
        long_rows,
        ["gene", "group", "average_expression", "percent_expressed"],
    )

    write_csv(
        out_dir / "Fig5G_dotplot_gene_order.csv",
        [{"order": i + 1, "gene": gene} for i, gene in enumerate(PANEL_G_GENE_ORDER)],
        ["order", "gene"],
    )

    write_csv(
        out_dir / "Fig5G_dotplot_note.csv",
        [
            {
                "note": "Average expression values were recovered from ALL-LC-matrix.csv. Percent-expressed values used for bubble size in panel G were not found as an explicit numeric table in this folder; panel visual exists in raw_support/Fig5G_dotplot_panel.pdf.",
                "raw_panel": "raw_support/Fig5G_dotplot_panel.pdf",
                "average_matrix": "raw_support/Fig5G_ALL-LC-matrix.csv",
            }
        ],
        ["note", "raw_panel", "average_matrix"],
    )


def build_panel_h_i(raw_dir: Path, out_dir: Path):
    rows = read_csv_rows(raw_dir / "Fig5I_Violin_Input_generated.csv")

    # Panel H counts from age rows
    age_counts = Counter()
    for row in rows:
        if str(row.get("Metric", "")).strip() == "age":
            g = str(row.get("Group", "")).strip()
            age_counts[g] += 1

    display_map = {"Control": "Normal", "OA": "OA", "NOA": "iNOA"}
    h_rows = []
    for g in ["Control", "OA", "NOA"]:
        h_rows.append(
            {
                "group_raw": g,
                "display_group": display_map.get(g, g),
                "n_patients": age_counts.get(g, 0),
            }
        )

    write_csv(out_dir / "Fig5H_cohort_counts.csv", h_rows, ["group_raw", "display_group", "n_patients"])

    # Panel I hormone values + summaries
    hormone_rows = []
    for row in rows:
        metric = str(row.get("Metric", "")).strip()
        group = str(row.get("Group", "")).strip()
        value = to_float(row.get("Value"))
        if metric in PANEL_I_METRICS and group in GROUP_ORDER_HORMONE and value is not None:
            hormone_rows.append({"Group": group, "Metric": metric, "Value": value})

    hormone_rows.sort(
        key=lambda r: (
            PANEL_I_METRICS.index(r["Metric"]),
            GROUP_ORDER_HORMONE.index(r["Group"]),
            r["Value"],
        )
    )

    write_csv(out_dir / "Fig5I_hormone_values.csv", hormone_rows, ["Group", "Metric", "Value"])

    grouped = defaultdict(list)
    for row in hormone_rows:
        grouped[(row["Metric"], row["Group"])].append(row["Value"])

    summary = []
    for metric in PANEL_I_METRICS:
        for group in GROUP_ORDER_HORMONE:
            vals = grouped.get((metric, group), [])
            summary.append(
                {
                    "Metric": metric,
                    "Group": group,
                    "n": len(vals),
                    "mean": (sum(vals) / len(vals)) if vals else "",
                    "median": quantile(vals, 0.5),
                    "q1": quantile(vals, 0.25),
                    "q3": quantile(vals, 0.75),
                }
            )

    write_csv(
        out_dir / "Fig5I_hormone_summary.csv",
        summary,
        ["Metric", "Group", "n", "mean", "median", "q1", "q3"],
    )


def main():
    code_dir = Path(__file__).resolve().parent
    out_dir = code_dir.parent / "figure5_source_data"
    raw_dir = out_dir / "raw_support"

    expression_rows = read_csv_rows(raw_dir / "Fig5A_expression.csv")

    build_panel_a_c_d_e(expression_rows, out_dir)
    build_panel_b(raw_dir, out_dir)
    build_panel_f(raw_dir, out_dir)
    build_panel_g(raw_dir, out_dir)
    build_panel_h_i(raw_dir, out_dir)

    print("Figure 5 source-data tables generated in:", out_dir)


if __name__ == "__main__":
    main()
