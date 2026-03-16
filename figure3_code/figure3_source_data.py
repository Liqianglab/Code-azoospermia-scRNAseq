#!/usr/bin/env python3
"""Build Figure 3 source-data tables from local raw-support files.

This script only uses files under:
  ../figure3_source_data/raw_support/
and writes normalized CSV tables to:
  ../figure3_source_data/
"""

from __future__ import annotations

import csv
import io
import math
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
CELL_REF_RE = re.compile(r"^([A-Z]+)([0-9]+)$")


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


def normalize_group(value: str) -> str:
    token = str(value).strip()
    upper = token.upper()
    if upper in {"CTRL", "CONTROL"}:
        return "Ctrl"
    if upper in {"OA"}:
        return "OA"
    if "AZFC" in upper:
        return "AZFc_Del"
    if "INOA_B" in upper:
        return "iNOA_B"
    if "INOA_S" in upper or "INOS_S" in upper:
        return "iNOA_S"
    if upper == "KS":
        return "KS"
    return token


def write_csv(path: Path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_csv_rows(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def column_index(col_letters: str) -> int:
    idx = 0
    for char in col_letters:
        idx = idx * 26 + (ord(char) - ord("A") + 1)
    return idx - 1


def read_xlsx_sheet1_rows_from_bytes(xlsx_bytes: bytes):
    with zipfile.ZipFile(io.BytesIO(xlsx_bytes)) as xlsx:
        return read_xlsx_sheet1_rows(xlsx)


def read_xlsx_sheet1_rows(xlsx: zipfile.ZipFile):
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


def read_xlsx_file(path: Path):
    with zipfile.ZipFile(path) as xlsx:
        return read_xlsx_sheet1_rows(xlsx)


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


def build_fig3a_somatic_umap(raw_support: Path, out_dir: Path):
    path = raw_support / "Fig3A_cell_annotation_from_clusters.csv"
    rows = read_csv_rows(path)

    somatic_map = {
        "ECs": "Endothelial cells",
        "LCs": "Leydig cells",
        "STs": "Sertoli cells",
        "Myoid": "Myoid cells",
        "Myeloid": "Immune cells",
        "Lym": "Immune cells",
    }

    cell_rows = []
    count_map = defaultdict(int)

    for row in rows:
        original_type = row.get("cell_type", "")
        if original_type not in somatic_map:
            continue
        group = normalize_group(row.get("group", ""))
        somatic_type = somatic_map[original_type]
        cell_rows.append(
            {
                "cell_id": row.get("cell_id", ""),
                "sample_id": row.get("sample_id", ""),
                "group": group,
                "UMAP1": row.get("UMAP1", ""),
                "UMAP2": row.get("UMAP2", ""),
                "somatic_celltype": somatic_type,
                "original_celltype": original_type,
            }
        )
        count_map[(group, somatic_type)] += 1

    group_totals = defaultdict(int)
    for (group, _), count in count_map.items():
        group_totals[group] += count

    count_rows = []
    ratio_rows = []
    for (group, somatic_type), count in sorted(count_map.items()):
        total = group_totals[group]
        ratio = (count / total) if total else 0.0
        count_rows.append(
            {
                "group": group,
                "somatic_celltype": somatic_type,
                "n_cells": count,
            }
        )
        ratio_rows.append(
            {
                "group": group,
                "somatic_celltype": somatic_type,
                "ratio_in_group": f"{ratio:.8f}",
            }
        )

    write_csv(
        out_dir / "Fig3A_somatic_umap_cells.csv",
        cell_rows,
        ["cell_id", "sample_id", "group", "UMAP1", "UMAP2", "somatic_celltype", "original_celltype"],
    )
    write_csv(
        out_dir / "Fig3A_somatic_cell_counts_by_group.csv",
        count_rows,
        ["group", "somatic_celltype", "n_cells"],
    )
    write_csv(
        out_dir / "Fig3A_somatic_cell_ratio_by_group.csv",
        ratio_rows,
        ["group", "somatic_celltype", "ratio_in_group"],
    )


def build_fig3b_radar(raw_support: Path, out_dir: Path):
    xlsx_rows = read_xlsx_file(raw_support / "Fig3B_demo2_celltype_ratio_matrix.xlsx")
    if not xlsx_rows:
        return

    header = xlsx_rows[0]
    data_rows = xlsx_rows[1:]
    if len(header) < 2:
        return

    disease_col = header[0]
    celltype_cols = header[1:]

    wide_rows = []
    long_rows = []

    for row in data_rows:
        if not row:
            continue
        disease = normalize_group(row[0] if len(row) > 0 else "")
        if disease == "":
            continue

        wide = {"disease": disease}
        for idx, celltype in enumerate(celltype_cols, start=1):
            value = row[idx] if idx < len(row) else ""
            wide[celltype] = value
            long_rows.append(
                {
                    "disease": disease,
                    "celltype": celltype,
                    "value": value,
                }
            )
        wide_rows.append(wide)

    disease_order = ["AZFc_Del", "iNOA_B", "iNOA_S", "KS", "OA"]
    celltype_order = ["LCs", "STs", "ECs", "Myoids", "Immune"]

    wide_rows.sort(key=lambda row: disease_order.index(row["disease"]) if row["disease"] in disease_order else 999)
    long_rows.sort(
        key=lambda row: (
            disease_order.index(row["disease"]) if row["disease"] in disease_order else 999,
            celltype_order.index(row["celltype"]) if row["celltype"] in celltype_order else 999,
        )
    )

    transposed_rows = []
    for celltype in celltype_cols:
        trow = {"celltype": celltype}
        for wide in wide_rows:
            trow[wide["disease"]] = wide.get(celltype, "")
        transposed_rows.append(trow)

    transposed_rows.sort(key=lambda row: celltype_order.index(row["celltype"]) if row["celltype"] in celltype_order else 999)

    write_csv(
        out_dir / "Fig3B_radar_matrix_group_by_celltype.csv",
        wide_rows,
        ["disease"] + celltype_cols,
    )
    write_csv(
        out_dir / "Fig3B_radar_matrix_celltype_by_group.csv",
        transposed_rows,
        ["celltype"] + [row["disease"] for row in wide_rows],
    )
    write_csv(
        out_dir / "Fig3B_radar_matrix_long.csv",
        long_rows,
        ["disease", "celltype", "value"],
    )


def infer_group_from_venn_header(header_value: str) -> str:
    token = header_value.upper()
    if "AZFC" in token:
        return "AZFc_Del"
    if "INOA_B" in token:
        return "iNOA_B"
    if "INOA_S" in token:
        return "iNOA_S"
    if token.endswith("KS") or "-KS" in token or "_KS" in token:
        return "KS"
    if token.endswith("OA") or "-OA" in token or "_OA" in token:
        return "OA"
    return header_value


def build_fig3c_venn(raw_support: Path, out_dir: Path):
    venn_sources = {
        "LC": raw_support / "Fig3C_LC_venn_Demo_Venn_5.xlsx",
        "ST": raw_support / "Fig3C_ST_venn_Demo_Venn_5.xlsx",
        "EC": raw_support / "Fig3C_EC_venn_Demo_Venn_5.xlsx",
        "Myoid": raw_support / "Fig3C_Myoid_venn_Demo_Venn_5.xlsx",
        "Immune": raw_support / "Fig3C_Immune_venn_Demo_Venn_5.xlsx",
    }

    long_rows = []
    count_map = defaultdict(set)

    for panel_celltype, xlsx_path in venn_sources.items():
        rows = read_xlsx_file(xlsx_path)
        if not rows:
            continue
        header = rows[0]
        table_rows = rows[1:]

        max_cols = max(len(row) for row in rows) if rows else len(header)
        normalized_header = list(header) + [""] * (max_cols - len(header))

        wide_rows = []
        for row in table_rows:
            padded = list(row) + [""] * (max_cols - len(row))
            wide_rows.append({normalized_header[idx]: padded[idx] for idx in range(max_cols)})

            for col_idx, gene in enumerate(padded):
                gene_text = str(gene).strip()
                if gene_text == "":
                    continue
                col_name = normalized_header[col_idx] if col_idx < len(normalized_header) else ""
                group = infer_group_from_venn_header(col_name)
                long_rows.append(
                    {
                        "panel_celltype": panel_celltype,
                        "group": group,
                        "gene": gene_text,
                        "source_column": col_name,
                    }
                )
                count_map[(panel_celltype, group)].add(gene_text)

        write_csv(
            out_dir / f"Fig3C_venn_{panel_celltype}_wide.csv",
            wide_rows,
            normalized_header,
        )

    count_rows = []
    for (panel_celltype, group), genes in sorted(count_map.items()):
        count_rows.append(
            {
                "panel_celltype": panel_celltype,
                "group": group,
                "n_unique_genes": len(genes),
            }
        )

    long_rows.sort(key=lambda row: (row["panel_celltype"], row["group"], row["gene"]))
    write_csv(
        out_dir / "Fig3C_venn_gene_sets_long.csv",
        long_rows,
        ["panel_celltype", "group", "gene", "source_column"],
    )
    write_csv(
        out_dir / "Fig3C_venn_group_gene_counts.csv",
        count_rows,
        ["panel_celltype", "group", "n_unique_genes"],
    )


def readable_pathway(description: str) -> str:
    cleaned = description.replace("HALLMARK_", "")
    parts = [token for token in cleaned.split("_") if token]
    return " ".join(token.capitalize() for token in parts)


def build_fig3c_bubble(raw_support: Path, out_dir: Path):
    bubble_sources = {
        "LC": raw_support / "Fig3C_bubble_table_PKD1_down_LC.csv",
        "ST": raw_support / "Fig3C_bubble_table_PKD1_down_ST.csv",
        "EC": raw_support / "Fig3C_bubble_table_PKD1_down_EC.csv",
        "Myoid": raw_support / "Fig3C_bubble_table_PKD1_down_Myoid.csv",
        "Immune": raw_support / "Fig3C_bubble_table_PKD1_down_Immune.csv",
    }

    rows = []
    for panel_celltype, csv_path in bubble_sources.items():
        for raw in read_csv_rows(csv_path):
            p_adjust = to_float(raw.get("p.adjust"))
            neglog10_fdr = ""
            if p_adjust is not None:
                neglog10_fdr = min(-math.log10(max(p_adjust, 1e-300)), 20.0)

            pathway = raw.get("Description", "")
            rows.append(
                {
                    "panel_celltype": panel_celltype,
                    "group": normalize_group(raw.get("group", "")),
                    "pathway": pathway,
                    "pathway_readable": readable_pathway(pathway),
                    "setSize": raw.get("setSize", ""),
                    "NES": raw.get("NES", ""),
                    "pvalue": raw.get("pvalue", ""),
                    "p.adjust": raw.get("p.adjust", ""),
                    "qvalues": raw.get("qvalues", ""),
                    "rank": raw.get("rank", ""),
                    "leading_edge": raw.get("leading_edge", ""),
                    "neglog10_fdr": neglog10_fdr,
                }
            )

    rows.sort(key=lambda row: (row["panel_celltype"], row["group"], to_float(row["p.adjust"]) or 999.0, row["pathway"]))
    write_csv(
        out_dir / "Fig3C_bubble_pathways_long.csv",
        rows,
        [
            "panel_celltype",
            "group",
            "pathway",
            "pathway_readable",
            "setSize",
            "NES",
            "pvalue",
            "p.adjust",
            "qvalues",
            "rank",
            "leading_edge",
            "neglog10_fdr",
        ],
    )

    top_rows = []
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["panel_celltype"], row["group"])].append(row)
    for key, values in grouped.items():
        values_sorted = sorted(values, key=lambda row: to_float(row["p.adjust"]) or 999.0)
        top_rows.extend(values_sorted[:10])
    write_csv(
        out_dir / "Fig3C_bubble_top10_pathways_by_group.csv",
        top_rows,
        [
            "panel_celltype",
            "group",
            "pathway",
            "pathway_readable",
            "setSize",
            "NES",
            "pvalue",
            "p.adjust",
            "qvalues",
            "rank",
            "leading_edge",
            "neglog10_fdr",
        ],
    )


def build_violin_score_tables(raw_csv: Path, score_column: str, output_prefix: str, out_dir: Path):
    raw_rows = read_csv_rows(raw_csv)
    cell_rows = []
    grouped_scores = defaultdict(list)

    for row in raw_rows:
        score = to_float(row.get(score_column))
        if score is None:
            continue
        group = normalize_group(row.get("gname", ""))
        if group == "":
            continue
        grouped_scores[group].append(score)
        cell_rows.append(
            {
                "cell_id": row.get("cell", ""),
                "group": group,
                "score": f"{score:.10g}",
                "UMAP1": row.get("UMAP1", ""),
                "UMAP2": row.get("UMAP2", ""),
                "TSNE1": row.get("TSNE1", ""),
                "TSNE2": row.get("TSNE2", ""),
            }
        )

    group_order = ["Ctrl", "OA", "AZFc_Del", "iNOA_B", "iNOA_S", "KS"]
    cell_rows.sort(key=lambda row: (group_order.index(row["group"]) if row["group"] in group_order else 999, row["cell_id"]))

    summary_rows = []
    for group, values in grouped_scores.items():
        values_sorted = sorted(values)
        summary_rows.append(
            {
                "group": group,
                "n_cells": len(values_sorted),
                "mean": f"{sum(values_sorted) / len(values_sorted):.10g}",
                "median": f"{quantile(values_sorted, 0.5):.10g}",
                "q1": f"{quantile(values_sorted, 0.25):.10g}",
                "q3": f"{quantile(values_sorted, 0.75):.10g}",
                "min": f"{values_sorted[0]:.10g}",
                "max": f"{values_sorted[-1]:.10g}",
            }
        )

    summary_rows.sort(key=lambda row: group_order.index(row["group"]) if row["group"] in group_order else 999)

    write_csv(
        out_dir / f"{output_prefix}_cells.csv",
        cell_rows,
        ["cell_id", "group", "score", "UMAP1", "UMAP2", "TSNE1", "TSNE2"],
    )
    write_csv(
        out_dir / f"{output_prefix}_group_summary.csv",
        summary_rows,
        ["group", "n_cells", "mean", "median", "q1", "q3", "min", "max"],
    )


def build_fig3f_heatmap(raw_support: Path, out_dir: Path):
    zip_path = raw_support / "Fig3F_sasp-gene.zip"
    with zipfile.ZipFile(zip_path) as outer_zip:
        xlsx_bytes = outer_zip.read("heatmap_export/heatmap_export.xlsx")

    rows = read_xlsx_sheet1_rows_from_bytes(xlsx_bytes)
    if not rows:
        return

    header = rows[0]
    value_cols = header[1:]
    data_rows = rows[1:]

    wide_rows = []
    long_rows = []
    for row in data_rows:
        if not row:
            continue
        gene = row[0].strip() if len(row) > 0 else ""
        if gene == "":
            continue
        wide = {"gene": gene}
        for idx, group in enumerate(value_cols, start=1):
            value = row[idx] if idx < len(row) else ""
            wide[group] = value
            long_rows.append(
                {
                    "gene": gene,
                    "group": group,
                    "value": value,
                }
            )
        wide_rows.append(wide)

    write_csv(
        out_dir / "Fig3F_heatmap_matrix_wide.csv",
        wide_rows,
        ["gene"] + value_cols,
    )
    write_csv(
        out_dir / "Fig3F_heatmap_matrix_long.csv",
        long_rows,
        ["gene", "group", "value"],
    )


def main():
    script_dir = Path(__file__).resolve().parent
    source_data_dir = script_dir.parent
    out_dir = source_data_dir / "figure3_source_data"
    raw_support = out_dir / "raw_support"

    out_dir.mkdir(parents=True, exist_ok=True)

    build_fig3a_somatic_umap(raw_support, out_dir)
    build_fig3b_radar(raw_support, out_dir)
    build_fig3c_venn(raw_support, out_dir)
    build_fig3c_bubble(raw_support, out_dir)
    build_violin_score_tables(
        raw_support / "Fig3D_REACTOME_CYTOKINE_SIGNALING_expression.csv",
        "1REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM.csv",
        "Fig3D_cytokine_score",
        out_dir,
    )
    build_violin_score_tables(
        raw_support / "Fig3E_SASP_expression.csv",
        "SASP.csv",
        "Fig3E_sasp_score",
        out_dir,
    )
    build_fig3f_heatmap(raw_support, out_dir)

    print("Figure 3 source-data tables generated in:", out_dir)


if __name__ == "__main__":
    main()
