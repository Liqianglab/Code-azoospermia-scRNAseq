#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
CELL_REF_RE = re.compile(r"^([A-Z]+)([0-9]+)$")

GROUP_ORDER = ["Ctrl", "AZFc_Del", "iNOA_B", "iNOA_S", "KS"]
STAGE_ORDER = ["Stage_a", "Stage_b", "Stage_c"]
ST_CLUSTER_ORDER = ["ST1", "ST2", "ST3"]
STAGE_CANDIDATE_MAP = {"ST1": "Stage_a", "ST2": "Stage_b", "ST3": "Stage_c"}
ST_ALIAS_TO_STAGE = {"ST_a": "Stage_a", "ST_b": "Stage_b", "ST_c": "Stage_c"}
FIG4F_GROUP_ORDER = ["Ctrl", "AZFc_Del", "iNOA_B", "iNOA_S", "KS"]
FIG4F_REGULON_ORDER = [
    "JUN(23g)",
    "FOS(37g)",
    "FOSB(38g)",
    "JUND(13g)",
    "HSF1(47g)",
    "MAFB(8g)",
    "NR2F2(40g)",
    "THRB(36g)",
    "ETS2(97g)",
    "ELF1(46g)",
    "E2F6(31g)",
    "NFKB1(26g)",
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
    if token.startswith("OA"):
        return "OA"
    if token.startswith("KS"):
        return "KS"
    match = re.match(r"^(Ctrl|OA|AZFc_Del|iNOA_B|iNOA_S|KS)", token)
    if match:
        return normalize_group(match.group(1))
    return ""


def read_csv_rows(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_tsv_rows(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_tsv_values(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.reader(handle, delimiter="\t"))


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


def sort_cluster(value: str):
    return ST_CLUSTER_ORDER.index(value) if value in ST_CLUSTER_ORDER else 999


def zscore(values):
    vals = [v for v in values if v is not None]
    if not vals:
        return [None for _ in values]
    mean_val = sum(vals) / len(vals)
    var = sum((v - mean_val) ** 2 for v in vals) / len(vals)
    sd = math.sqrt(var)
    if sd == 0:
        return [0.0 if v is not None else None for v in values]
    return [((v - mean_val) / sd) if v is not None else None for v in values]


def clean_field(value):
    return str(value).strip().strip('"')


def parse_pyscenic_group_table(path: Path):
    rows = read_tsv_values(path)
    if not rows:
        return [], [], {}

    header = [clean_field(col) for col in rows[0]]
    groups = [normalize_group(col) for col in header[1:]]
    labels = []
    table = {}

    for row in rows[1:]:
        if not row:
            continue
        label = clean_field(row[0])
        if not label:
            continue
        labels.append(label)
        values = {}
        for idx, group in enumerate(groups, start=1):
            raw_value = row[idx] if idx < len(row) else ""
            values[group] = to_float(clean_field(raw_value))
        table[label] = values

    return groups, labels, table


def build_panel_a(raw_support: Path, out_dir: Path):
    pseudo_rows = read_tsv_rows(raw_support / "Fig4A_monocle_pseudotime_cells.xls")
    coord_rows = read_tsv_rows(raw_support / "Fig4A_monocle_plot_coords.xls")
    line_rows = read_tsv_rows(raw_support / "Fig4A_monocle_trajectory_lines.xls")

    coord_map = {}
    for row in coord_rows:
        cell_id = row.get("cellID", "")
        if cell_id:
            coord_map[cell_id] = {
                "traj_x": row.get("X1", ""),
                "traj_y": row.get("X2", ""),
            }

    cell_rows = []
    count_group_stage = defaultdict(int)
    count_stage_state = defaultdict(int)
    stage_pseudotime_values = defaultdict(list)

    for row in pseudo_rows:
        st_cluster = row.get("celltype", "").strip()
        if st_cluster not in STAGE_CANDIDATE_MAP:
            continue
        cell_id = row.get("barcode", "")
        group = normalize_group(row.get("group", ""))
        stage_label = STAGE_CANDIDATE_MAP.get(st_cluster, "")
        state = row.get("State", "")
        pseudotime = to_float(row.get("Pseudotime", ""))
        coords = coord_map.get(cell_id, {})

        out_row = {
            "cell_id": cell_id,
            "sample_id": row.get("sample", ""),
            "orig_ident": row.get("orig.ident", ""),
            "group": group,
            "st_cluster": st_cluster,
            "stage_label_candidate": stage_label,
            "state": state,
            "pseudotime": row.get("Pseudotime", ""),
            "traj_x": coords.get("traj_x", ""),
            "traj_y": coords.get("traj_y", ""),
            "nCount_RNA": row.get("nCount_RNA", ""),
            "nFeature_RNA": row.get("nFeature_RNA", ""),
            "Cluster": row.get("Cluster", ""),
            "Size_Factor": row.get("Size_Factor", ""),
            "num_genes_expressed": row.get("num_genes_expressed", ""),
            "UMI": row.get("UMI", ""),
        }
        cell_rows.append(out_row)

        count_group_stage[(group, st_cluster, stage_label)] += 1
        count_stage_state[(st_cluster, stage_label, state)] += 1
        if pseudotime is not None:
            stage_pseudotime_values[(st_cluster, stage_label)].append(pseudotime)

    cell_rows.sort(
        key=lambda row: (
            sort_group(row["group"]),
            sort_cluster(row["st_cluster"]),
            row["cell_id"],
        )
    )

    group_totals = defaultdict(int)
    for (group, _, _), count in count_group_stage.items():
        group_totals[group] += count

    group_stage_rows = []
    for (group, st_cluster, stage_label), count in sorted(
        count_group_stage.items(),
        key=lambda item: (sort_group(item[0][0]), sort_cluster(item[0][1])),
    ):
        total = group_totals[group]
        ratio = count / total if total else 0.0
        group_stage_rows.append(
            {
                "group": group,
                "st_cluster": st_cluster,
                "stage_label_candidate": stage_label,
                "n_cells": count,
                "ratio_in_group": f"{ratio:.8f}",
            }
        )

    stage_totals = defaultdict(int)
    for (st_cluster, stage_label, _state), count in count_stage_state.items():
        stage_totals[(st_cluster, stage_label)] += count

    stage_state_rows = []
    for (st_cluster, stage_label, state), count in sorted(
        count_stage_state.items(),
        key=lambda item: (
            sort_cluster(item[0][0]),
            item[0][2],
        ),
    ):
        total = stage_totals[(st_cluster, stage_label)]
        ratio = count / total if total else 0.0
        stage_state_rows.append(
            {
                "st_cluster": st_cluster,
                "stage_label_candidate": stage_label,
                "state": state,
                "n_cells": count,
                "ratio_in_cluster": f"{ratio:.8f}",
            }
        )

    pseudotime_rows = []
    for (st_cluster, stage_label), values in sorted(
        stage_pseudotime_values.items(), key=lambda item: sort_cluster(item[0][0])
    ):
        pseudotime_rows.append(
            {
                "st_cluster": st_cluster,
                "stage_label_candidate": stage_label,
                "n_cells": len(values),
                "pseudotime_mean": f"{(sum(values) / len(values)):.8f}" if values else "",
                "pseudotime_median": f"{quantile(values, 0.5):.8f}" if values else "",
                "pseudotime_q1": f"{quantile(values, 0.25):.8f}" if values else "",
                "pseudotime_q3": f"{quantile(values, 0.75):.8f}" if values else "",
            }
        )

    line_out_rows = []
    for row in line_rows:
        line_out_rows.append(
            {
                "line_id": row.get("ID", ""),
                "traj_x": row.get("X1", ""),
                "traj_y": row.get("X2", ""),
            }
        )

    stage_map_rows = [
        {"st_cluster": cluster, "stage_label_candidate": stage}
        for cluster, stage in STAGE_CANDIDATE_MAP.items()
    ]

    write_csv(
        out_dir / "Fig4A_ST_cells_trajectory.csv",
        cell_rows,
        [
            "cell_id",
            "sample_id",
            "orig_ident",
            "group",
            "st_cluster",
            "stage_label_candidate",
            "state",
            "pseudotime",
            "traj_x",
            "traj_y",
            "nCount_RNA",
            "nFeature_RNA",
            "Cluster",
            "Size_Factor",
            "num_genes_expressed",
            "UMI",
        ],
    )
    write_csv(
        out_dir / "Fig4A_trajectory_line_points.csv",
        line_out_rows,
        ["line_id", "traj_x", "traj_y"],
    )
    write_csv(
        out_dir / "Fig4A_stage_candidate_mapping.csv",
        stage_map_rows,
        ["st_cluster", "stage_label_candidate"],
    )
    write_csv(
        out_dir / "Fig4A_ST_counts_by_group.csv",
        group_stage_rows,
        ["group", "st_cluster", "stage_label_candidate", "n_cells", "ratio_in_group"],
    )
    write_csv(
        out_dir / "Fig4A_ST_counts_by_state.csv",
        stage_state_rows,
        ["st_cluster", "stage_label_candidate", "state", "n_cells", "ratio_in_cluster"],
    )
    write_csv(
        out_dir / "Fig4A_ST_pseudotime_summary.csv",
        pseudotime_rows,
        [
            "st_cluster",
            "stage_label_candidate",
            "n_cells",
            "pseudotime_mean",
            "pseudotime_median",
            "pseudotime_q1",
            "pseudotime_q3",
        ],
    )


def build_panel_b(raw_support: Path, out_dir: Path):
    xlsx_rows = read_xlsx_sheet1_rows(raw_support / "Fig4B_heatmap_export.xlsx")
    if not xlsx_rows:
        return
    header = xlsx_rows[0]
    if len(header) < 4:
        return

    stage_columns = [str(token).strip() for token in header[1:] if str(token).strip() != ""]
    wide_rows = []
    long_rows = []
    order_rows = []

    for idx, row in enumerate(xlsx_rows[1:], start=1):
        if not row:
            continue
        gene = str(row[0]).strip() if len(row) > 0 else ""
        if gene == "":
            continue
        value_map = {}
        for col_idx, stage_col in enumerate(stage_columns, start=1):
            raw_val = row[col_idx] if col_idx < len(row) else ""
            value_map[stage_col] = raw_val
            stage_label = STAGE_CANDIDATE_MAP.get(stage_col, "")
            long_rows.append(
                {
                    "gene": gene,
                    "st_cluster": stage_col,
                    "stage_label_candidate": stage_label,
                    "z_score": raw_val,
                }
            )

        wide_row = {"gene": gene}
        for stage_col in stage_columns:
            wide_row[stage_col] = value_map.get(stage_col, "")
        wide_row["Stage_a"] = value_map.get("ST1", "")
        wide_row["Stage_b"] = value_map.get("ST2", "")
        wide_row["Stage_c"] = value_map.get("ST3", "")
        wide_rows.append(wide_row)
        order_rows.append({"gene_order": idx, "gene": gene})

    write_csv(
        out_dir / "Fig4B_heatmap_matrix_wide.csv",
        wide_rows,
        ["gene", "ST1", "ST2", "ST3", "Stage_a", "Stage_b", "Stage_c"],
    )
    write_csv(
        out_dir / "Fig4B_heatmap_matrix_long.csv",
        long_rows,
        ["gene", "st_cluster", "stage_label_candidate", "z_score"],
    )
    write_csv(
        out_dir / "Fig4B_heatmap_gene_order.csv",
        order_rows,
        ["gene_order", "gene"],
    )


def load_stage_and_group_from_panel_a(out_dir: Path):
    rows = read_csv_rows(out_dir / "Fig4A_ST_cells_trajectory.csv")
    stage_map = {}
    group_map = {}
    for row in rows:
        cell_id = row["cell_id"]
        stage_map[cell_id] = row.get("stage_label_candidate", "")
        group_map[cell_id] = row.get("group", "")
    return stage_map, group_map


def load_expression_maps(raw_support: Path):
    expr_v1 = read_csv_rows(raw_support / "Fig4_expression_markers_v1.csv")
    expr_v2 = read_csv_rows(raw_support / "Fig4_expression_markers_v2.csv")

    value_map = defaultdict(dict)
    group_map = {}
    cluster_map = {}

    for row in expr_v1:
        cell_id = row.get("cell", "")
        if cell_id == "":
            continue
        group_map[cell_id] = normalize_group(row.get("gname", ""))
        cluster_map[cell_id] = row.get("Cluster", "")
        for key, value in row.items():
            if "normalised expression value" in key:
                val = to_float(value)
                if val is not None:
                    gene = key.replace(" normalised expression value", "").strip()
                    value_map[cell_id][gene] = val

    for row in expr_v2:
        cell_id = row.get("cell", "")
        if cell_id == "":
            continue
        if cell_id not in group_map:
            group_map[cell_id] = normalize_group(row.get("gname", ""))
        if cell_id not in cluster_map:
            cluster_map[cell_id] = row.get("Cluster", "")
        for key, value in row.items():
            if "normalised expression value" in key:
                val = to_float(value)
                if val is not None:
                    gene = key.replace(" normalised expression value", "").strip()
                    value_map[cell_id][gene] = val

    return value_map, group_map, cluster_map


def build_panel_c(raw_support: Path, out_dir: Path):
    stage_map, group_from_a = load_stage_and_group_from_panel_a(out_dir)
    value_map, group_map, cluster_map = load_expression_maps(raw_support)

    panel_genes = ["FOS", "JUN", "HOPX", "DEFB119"]
    cell_rows = []
    summary_values = defaultdict(list)

    for cell_id, gene_values in value_map.items():
        if cell_id not in stage_map:
            continue
        stage_label = stage_map[cell_id]
        group = group_map.get(cell_id, group_from_a.get(cell_id, ""))
        st_cluster = cluster_map.get(cell_id, "")
        for gene in panel_genes:
            if gene not in gene_values:
                continue
            value = gene_values[gene]
            cell_rows.append(
                {
                    "cell_id": cell_id,
                    "group": group,
                    "st_cluster": st_cluster,
                    "stage_label_candidate": stage_label,
                    "gene": gene,
                    "expression": f"{value:.8f}",
                }
            )
            summary_values[(gene, stage_label)].append(value)

    cell_rows.sort(
        key=lambda row: (
            panel_genes.index(row["gene"]) if row["gene"] in panel_genes else 999,
            sort_stage(row["stage_label_candidate"]),
            sort_group(row["group"]),
            row["cell_id"],
        )
    )

    summary_rows = []
    for gene in panel_genes:
        for stage_label in STAGE_ORDER:
            values = summary_values.get((gene, stage_label), [])
            if not values:
                continue
            summary_rows.append(
                {
                    "gene": gene,
                    "stage_label_candidate": stage_label,
                    "n_cells": len(values),
                    "mean_expression": f"{(sum(values) / len(values)):.8f}",
                    "median_expression": f"{quantile(values, 0.5):.8f}",
                    "q1_expression": f"{quantile(values, 0.25):.8f}",
                    "q3_expression": f"{quantile(values, 0.75):.8f}",
                }
            )

    write_csv(
        out_dir / "Fig4C_stage_marker_cells.csv",
        cell_rows,
        ["cell_id", "group", "st_cluster", "stage_label_candidate", "gene", "expression"],
    )
    write_csv(
        out_dir / "Fig4C_stage_marker_summary.csv",
        summary_rows,
        [
            "gene",
            "stage_label_candidate",
            "n_cells",
            "mean_expression",
            "median_expression",
            "q1_expression",
            "q3_expression",
        ],
    )


def build_panel_d(raw_support: Path, out_dir: Path):
    raw_rows = read_csv_rows(raw_support / "Fig4D_cluster_counts_by_group_raw.csv")
    count_rows = []
    group_totals = defaultdict(int)
    stage_totals = defaultdict(int)

    for row in raw_rows:
        st_cluster = row.get("Cluster", "").strip()
        if st_cluster not in STAGE_CANDIDATE_MAP:
            continue
        group = normalize_group(row.get("Sample", ""))
        n_cells = int(float(row.get("Number", "0") or 0))
        stage_label = STAGE_CANDIDATE_MAP[st_cluster]
        count_rows.append(
            {
                "group": group,
                "st_cluster": st_cluster,
                "stage_label_candidate": stage_label,
                "n_cells": n_cells,
            }
        )
        group_totals[group] += n_cells
        stage_totals[st_cluster] += n_cells

    count_rows.sort(key=lambda row: (sort_group(row["group"]), sort_cluster(row["st_cluster"])))

    ratio_rows = []
    for row in count_rows:
        total = group_totals[row["group"]]
        ratio = row["n_cells"] / total if total else 0.0
        ratio_rows.append(
            {
                **row,
                "ratio_in_group": f"{ratio:.8f}",
            }
        )

    all_group_rows = []
    all_total = sum(stage_totals.values())
    for st_cluster in ST_CLUSTER_ORDER:
        n_cells = stage_totals.get(st_cluster, 0)
        ratio = n_cells / all_total if all_total else 0.0
        all_group_rows.append(
            {
                "group": "all_group",
                "st_cluster": st_cluster,
                "stage_label_candidate": STAGE_CANDIDATE_MAP[st_cluster],
                "n_cells": n_cells,
                "ratio_in_group": f"{ratio:.8f}",
            }
        )

    selected_rows = read_csv_rows(raw_support / "Fig4D_view_stacked_bars_selected_samples.csv")
    selected_out_rows = []
    selected_group_totals = defaultdict(int)

    sample_col = ""
    cluster_col = ""
    count_col = ""
    if selected_rows:
        columns = list(selected_rows[0].keys())
        if len(columns) >= 3:
            sample_col, cluster_col, count_col = columns[0], columns[1], columns[2]

    for row in selected_rows:
        st_cluster = row.get(cluster_col, "").strip()
        if st_cluster not in STAGE_CANDIDATE_MAP:
            continue
        group = normalize_group(row.get(sample_col, ""))
        n_cells = int(float(row.get(count_col, "0") or 0))
        selected_out_rows.append(
            {
                "group": group,
                "st_cluster": st_cluster,
                "stage_label_candidate": STAGE_CANDIDATE_MAP[st_cluster],
                "n_cells": n_cells,
            }
        )
        selected_group_totals[group] += n_cells

    selected_out_rows.sort(
        key=lambda row: (sort_group(row["group"]), sort_cluster(row["st_cluster"]))
    )

    selected_ratio_rows = []
    for row in selected_out_rows:
        total = selected_group_totals[row["group"]]
        ratio = row["n_cells"] / total if total else 0.0
        selected_ratio_rows.append(
            {
                **row,
                "ratio_in_group": f"{ratio:.8f}",
            }
        )

    write_csv(
        out_dir / "Fig4D_stage_counts_by_group.csv",
        count_rows,
        ["group", "st_cluster", "stage_label_candidate", "n_cells"],
    )
    write_csv(
        out_dir / "Fig4D_stage_ratio_by_group.csv",
        ratio_rows,
        ["group", "st_cluster", "stage_label_candidate", "n_cells", "ratio_in_group"],
    )
    write_csv(
        out_dir / "Fig4D_stage_ratio_all_group.csv",
        all_group_rows,
        ["group", "st_cluster", "stage_label_candidate", "n_cells", "ratio_in_group"],
    )
    write_csv(
        out_dir / "Fig4D_stage_counts_selected_samples.csv",
        selected_out_rows,
        ["group", "st_cluster", "stage_label_candidate", "n_cells"],
    )
    write_csv(
        out_dir / "Fig4D_stage_ratio_selected_samples.csv",
        selected_ratio_rows,
        ["group", "st_cluster", "stage_label_candidate", "n_cells", "ratio_in_group"],
    )


def build_panel_e(raw_support: Path, out_dir: Path):
    stage_map, group_from_a = load_stage_and_group_from_panel_a(out_dir)
    value_map, group_map, cluster_map = load_expression_maps(raw_support)

    panel_genes = ["JUN", "EGR3", "DEFB119", "CITED1"]
    cell_rows = []
    summary_values = defaultdict(list)

    for cell_id, gene_values in value_map.items():
        if cell_id not in stage_map:
            continue
        group = group_map.get(cell_id, group_from_a.get(cell_id, ""))
        st_cluster = cluster_map.get(cell_id, "")
        stage_label = stage_map[cell_id]
        for gene in panel_genes:
            if gene not in gene_values:
                continue
            value = gene_values[gene]
            cell_rows.append(
                {
                    "cell_id": cell_id,
                    "group": group,
                    "st_cluster": st_cluster,
                    "stage_label_candidate": stage_label,
                    "gene": gene,
                    "expression": f"{value:.8f}",
                }
            )
            summary_values[(gene, group)].append(value)

    cell_rows.sort(
        key=lambda row: (
            panel_genes.index(row["gene"]) if row["gene"] in panel_genes else 999,
            sort_group(row["group"]),
            row["cell_id"],
        )
    )

    summary_rows = []
    for gene in panel_genes:
        for group in GROUP_ORDER:
            values = summary_values.get((gene, group), [])
            if not values:
                continue
            summary_rows.append(
                {
                    "gene": gene,
                    "group": group,
                    "n_cells": len(values),
                    "mean_expression": f"{(sum(values) / len(values)):.8f}",
                    "median_expression": f"{quantile(values, 0.5):.8f}",
                    "q1_expression": f"{quantile(values, 0.25):.8f}",
                    "q3_expression": f"{quantile(values, 0.75):.8f}",
                }
            )

    write_csv(
        out_dir / "Fig4E_group_marker_cells.csv",
        cell_rows,
        ["cell_id", "group", "st_cluster", "stage_label_candidate", "gene", "expression"],
    )
    write_csv(
        out_dir / "Fig4E_group_marker_summary.csv",
        summary_rows,
        [
            "gene",
            "group",
            "n_cells",
            "mean_expression",
            "median_expression",
            "q1_expression",
            "q3_expression",
        ],
    )


def load_average_matrix(path: Path):
    rows = read_csv_rows(path)
    row_map = {}
    for row in rows:
        gene = row.get("gene_name", "").strip()
        if gene:
            row_map[gene] = row
    return row_map


def build_panel_f(raw_support: Path, out_dir: Path):
    scaled_path = raw_support / "Fig4F_pyscenic_regulonActivity_byCellType_Scaled_top.xls"
    auc_path = raw_support / "Fig4F_pyscenic_regulonActivity_CellType.xls"
    annotation_path = raw_support / "Fig4F_pyscenic_heatmap_top_annotation.xls"
    aucell_path = raw_support / "Fig4F_pyscenic_aucell.xls"

    if not (scaled_path.exists() and auc_path.exists() and annotation_path.exists() and aucell_path.exists()):
        build_panel_f_proxy(raw_support, out_dir)
        return "proxy"

    _, scaled_labels, scaled_table = parse_pyscenic_group_table(scaled_path)
    _, auc_labels, auc_table = parse_pyscenic_group_table(auc_path)
    auc_set = set(auc_labels)

    selected_labels = [label for label in FIG4F_REGULON_ORDER if label in scaled_table and label in auc_set]
    if not selected_labels:
        build_panel_f_proxy(raw_support, out_dir)
        return "proxy"

    wide_rows = []
    long_rows = []
    for label in selected_labels:
        tf_gene = label.split("(", 1)[0]
        wide_row = {"regulon_label": label, "tf_gene": tf_gene}
        for group in FIG4F_GROUP_ORDER:
            auc_val = auc_table.get(label, {}).get(group)
            scaled_val = scaled_table.get(label, {}).get(group)
            wide_row[f"{group}_auc"] = "" if auc_val is None else f"{auc_val:.8f}"
            wide_row[f"{group}_activity_zscore"] = "" if scaled_val is None else f"{scaled_val:.8f}"
            long_rows.append(
                {
                    "regulon_label": label,
                    "tf_gene": tf_gene,
                    "group": group,
                    "auc": "" if auc_val is None else f"{auc_val:.8f}",
                    "activity_zscore": "" if scaled_val is None else f"{scaled_val:.8f}",
                }
            )
        wide_rows.append(wide_row)

    candidate_rows = []
    for label in scaled_labels:
        candidate_rows.append(
            {
                "regulon_label": label,
                "tf_gene": label.split("(", 1)[0],
                "in_scaled_table": "1",
                "in_auc_table": "1" if label in auc_set else "0",
                "selected_for_fig4f": "1" if label in selected_labels else "0",
            }
        )

    note_rows = [
        {
            "panel": "Fig4F",
            "status": "raw",
            "note": "Numeric regulon matrices recovered from ST细胞提取-scenic1_regulon.tar.gz and extracted into raw_support. OA group is excluded in this submission export to match Figure 4 layout.",
            "raw_support_files": ";".join(
                [
                    "raw_support/Fig4F_pyscenic_regulonActivity_CellType.xls",
                    "raw_support/Fig4F_pyscenic_regulonActivity_byCellType_Scaled_top.xls",
                    "raw_support/Fig4F_pyscenic_heatmap_top_annotation.xls",
                    "raw_support/Fig4F_pyscenic_aucell.xls",
                    "raw_support/Fig4F_pyscenic_regulons.csv",
                ]
            ),
        }
    ]

    write_csv(
        out_dir / "Fig4F_regulon_activity_matrix_wide.csv",
        wide_rows,
        [
            "regulon_label",
            "tf_gene",
            "Ctrl_auc",
            "AZFc_Del_auc",
            "iNOA_B_auc",
            "iNOA_S_auc",
            "KS_auc",
            "Ctrl_activity_zscore",
            "AZFc_Del_activity_zscore",
            "iNOA_B_activity_zscore",
            "iNOA_S_activity_zscore",
            "KS_activity_zscore",
        ],
    )
    write_csv(
        out_dir / "Fig4F_regulon_activity_matrix_long.csv",
        long_rows,
        ["regulon_label", "tf_gene", "group", "auc", "activity_zscore"],
    )
    write_csv(
        out_dir / "Fig4F_regulon_candidates_from_scaled_top.csv",
        candidate_rows,
        ["regulon_label", "tf_gene", "in_scaled_table", "in_auc_table", "selected_for_fig4f"],
    )
    write_csv(
        out_dir / "Fig4F_note.csv",
        note_rows,
        ["panel", "status", "note", "raw_support_files"],
    )

    for stale_name in [
        "Fig4F_regulon_activity_proxy_matrix_wide.csv",
        "Fig4F_regulon_activity_proxy_matrix_long.csv",
        "Fig4F_proxy_note.csv",
    ]:
        stale_path = out_dir / stale_name
        if stale_path.exists():
            stale_path.unlink()

    return "raw"


def build_panel_f_proxy(raw_support: Path, out_dir: Path):
    matrix = load_average_matrix(raw_support / "Fig4_ST_average_expression_matrix.csv")
    groups = ["Ctrl", "AZFc_Del", "iNOA_B", "iNOA_S", "KS"]

    regulon_labels = [
        ("JUN", "JUN(239)"),
        ("FOS", "FOS(37)"),
        ("FOSB", "FOSB(39)"),
        ("JUND", "JUND(13)"),
        ("HSF1", "HSF1(47)"),
        ("MAFB", "MAFB(86)"),
        ("NR2F2", "NR2F2(40)"),
        ("THRB", "THRB(36)"),
        ("ETS2", "ETS2(97)"),
        ("ELF1", "ELF1(46)"),
        ("E2F3", "E2F3(19)"),
        ("NFKB1", "NFKB1(26)"),
    ]

    wide_rows = []
    long_rows = []

    for gene, label in regulon_labels:
        row = matrix.get(gene, {})
        values = [to_float(row.get(group, "")) for group in groups]
        z_vals = zscore(values)
        wide_row = {
            "regulon_label": label,
            "tf_gene": gene,
        }
        for group, val, z_val in zip(groups, values, z_vals):
            wide_row[f"{group}_avg_expr"] = "" if val is None else f"{val:.8f}"
            wide_row[f"{group}_zscore_proxy"] = "" if z_val is None else f"{z_val:.8f}"
            long_rows.append(
                {
                    "regulon_label": label,
                    "tf_gene": gene,
                    "group": group,
                    "avg_expr": "" if val is None else f"{val:.8f}",
                    "activity_zscore_proxy": "" if z_val is None else f"{z_val:.8f}",
                }
            )
        wide_rows.append(wide_row)

    long_rows.sort(
        key=lambda row: (
            [label for _, label in regulon_labels].index(row["regulon_label"]),
            groups.index(row["group"]) if row["group"] in groups else 999,
        )
    )

    note_rows = [
        {
            "panel": "Fig4F",
            "status": "proxy",
            "note": "Original numeric regulon AUC matrix file was not found in this workspace. Proxy matrix is computed from group-average TF expression in Fig4_ST_average_expression_matrix.csv and row-wise z-scored across groups. OA group is excluded in this submission export to match Figure 4 layout.",
            "raw_support_file": "raw_support/Fig4_ST_average_expression_matrix.csv",
        }
    ]

    write_csv(
        out_dir / "Fig4F_regulon_activity_proxy_matrix_wide.csv",
        wide_rows,
        [
            "regulon_label",
            "tf_gene",
            "Ctrl_avg_expr",
            "AZFc_Del_avg_expr",
            "iNOA_B_avg_expr",
            "iNOA_S_avg_expr",
            "KS_avg_expr",
            "Ctrl_zscore_proxy",
            "AZFc_Del_zscore_proxy",
            "iNOA_B_zscore_proxy",
            "iNOA_S_zscore_proxy",
            "KS_zscore_proxy",
        ],
    )
    write_csv(
        out_dir / "Fig4F_regulon_activity_proxy_matrix_long.csv",
        long_rows,
        ["regulon_label", "tf_gene", "group", "avg_expr", "activity_zscore_proxy"],
    )
    write_csv(
        out_dir / "Fig4F_proxy_note.csv",
        note_rows,
        ["panel", "status", "note", "raw_support_file"],
    )


def build_panel_g(raw_support: Path, out_dir: Path):
    raw_path = raw_support / "Fig4G_BTB_Integrity_Score_expression.csv"
    if not raw_path.exists():
        build_panel_g_proxy(raw_support, out_dir)
        return "proxy"

    rows = read_csv_rows(raw_path)
    if not rows:
        build_panel_g_proxy(raw_support, out_dir)
        return "proxy"

    columns = list(rows[0].keys())
    score_col = ""
    for column in columns:
        upper = column.upper()
        if "BTB" in upper and "INTEGRITY" in upper:
            score_col = column
            break
    if not score_col:
        for column in columns:
            upper = column.upper()
            if "BTB" in upper and "SCORE" in upper:
                score_col = column
                break
    if not score_col:
        build_panel_g_proxy(raw_support, out_dir)
        return "proxy"

    cell_rows = []
    group_values = defaultdict(list)
    group_stage_values = defaultdict(list)

    for row in rows:
        cell_id = row.get("cell", "")
        group = infer_group_from_cell_id(cell_id)
        if group not in GROUP_ORDER:
            continue
        major_cell_type = str(row.get("Major cell types", "")).strip()
        stage_label = ST_ALIAS_TO_STAGE.get(major_cell_type, "")
        btb_score = to_float(row.get(score_col, ""))

        out_row = {
            "cell_id": cell_id,
            "group": group,
            "major_cell_type": major_cell_type,
            "stage_label_candidate": stage_label,
            "btb_score": "" if btb_score is None else f"{btb_score:.8f}",
            "UMAP1": row.get("UMAP1", ""),
            "UMAP2": row.get("UMAP2", ""),
            "TSNE1": row.get("TSNE1", ""),
            "TSNE2": row.get("TSNE2", ""),
        }
        cell_rows.append(out_row)

        if btb_score is not None and group in GROUP_ORDER:
            group_values[group].append(btb_score)
            group_stage_values[(group, major_cell_type, stage_label)].append(btb_score)

    cell_rows.sort(key=lambda row: (sort_group(row["group"]), row["cell_id"]))

    summary_rows = []
    for group in GROUP_ORDER:
        values = group_values.get(group, [])
        if not values:
            continue
        summary_rows.append(
            {
                "group": group,
                "n_cells": len(values),
                "mean_btb_score": f"{(sum(values) / len(values)):.8f}",
                "median_btb_score": f"{quantile(values, 0.5):.8f}",
                "q1_btb_score": f"{quantile(values, 0.25):.8f}",
                "q3_btb_score": f"{quantile(values, 0.75):.8f}",
            }
        )

    stage_rows = []
    for (group, major_cell_type, stage_label), values in sorted(
        group_stage_values.items(),
        key=lambda item: (
            sort_group(item[0][0]),
            sort_stage(item[0][2]),
            item[0][1],
        ),
    ):
        stage_rows.append(
            {
                "group": group,
                "major_cell_type": major_cell_type,
                "stage_label_candidate": stage_label,
                "n_cells": len(values),
                "mean_btb_score": f"{(sum(values) / len(values)):.8f}",
                "median_btb_score": f"{quantile(values, 0.5):.8f}",
                "q1_btb_score": f"{quantile(values, 0.25):.8f}",
                "q3_btb_score": f"{quantile(values, 0.75):.8f}",
            }
        )

    note_rows = [
        {
            "panel": "Fig4G",
            "status": "raw",
            "note": "Cell-level BTB integrity scores recovered from ST细胞提取 - BTB_Integrity_Score_expression.csv. OA group is excluded in this submission export to match Figure 4 layout.",
            "raw_support_file": "raw_support/Fig4G_BTB_Integrity_Score_expression.csv",
            "score_column_used": score_col,
        }
    ]

    write_csv(
        out_dir / "Fig4G_BTB_score_cells.csv",
        cell_rows,
        [
            "cell_id",
            "group",
            "major_cell_type",
            "stage_label_candidate",
            "btb_score",
            "UMAP1",
            "UMAP2",
            "TSNE1",
            "TSNE2",
        ],
    )
    write_csv(
        out_dir / "Fig4G_BTB_group_score_summary.csv",
        summary_rows,
        ["group", "n_cells", "mean_btb_score", "median_btb_score", "q1_btb_score", "q3_btb_score"],
    )
    write_csv(
        out_dir / "Fig4G_BTB_stage_score_summary.csv",
        stage_rows,
        [
            "group",
            "major_cell_type",
            "stage_label_candidate",
            "n_cells",
            "mean_btb_score",
            "median_btb_score",
            "q1_btb_score",
            "q3_btb_score",
        ],
    )
    write_csv(
        out_dir / "Fig4G_note.csv",
        note_rows,
        ["panel", "status", "note", "raw_support_file", "score_column_used"],
    )

    for stale_name in [
        "Fig4G_BTB_signature_genes_used.csv",
        "Fig4G_BTB_gene_expression_by_group.csv",
        "Fig4G_BTB_group_score_proxy.csv",
        "Fig4G_proxy_note.csv",
    ]:
        stale_path = out_dir / stale_name
        if stale_path.exists():
            stale_path.unlink()

    return "raw"


def load_gene_list(path: Path):
    with path.open("r", encoding="utf-8-sig") as handle:
        lines = [line.strip() for line in handle.readlines()]
    return [line for line in lines if line and not line.startswith("#")]


def load_csv_gene_list(path: Path):
    rows = read_csv_rows(path)
    if not rows:
        return []
    first_col = list(rows[0].keys())[0]
    genes = []
    for row in rows:
        gene = str(row.get(first_col, "")).strip()
        if gene and gene != first_col:
            genes.append(gene)
    return genes


def build_panel_g_proxy(raw_support: Path, out_dir: Path):
    matrix = load_average_matrix(raw_support / "Fig4_ST_average_expression_matrix.csv")
    groups = ["Ctrl", "AZFc_Del", "iNOA_B", "iNOA_S", "KS"]

    signature_genes = load_csv_gene_list(raw_support / "Fig4G_BTB_signature_genes.csv")
    if not signature_genes:
        signature_genes = load_gene_list(raw_support / "Fig4G_BTB_core_for_heatmap.txt")

    gene_used_rows = []
    present_genes = []
    for gene in signature_genes:
        in_matrix = gene in matrix
        gene_used_rows.append(
            {
                "gene": gene,
                "in_ST_average_matrix": "1" if in_matrix else "0",
            }
        )
        if in_matrix:
            present_genes.append(gene)

    gene_group_rows = []
    for gene in present_genes:
        row = matrix[gene]
        for group in groups:
            val = to_float(row.get(group, ""))
            gene_group_rows.append(
                {
                    "gene": gene,
                    "group": group,
                    "avg_expr": "" if val is None else f"{val:.8f}",
                }
            )

    group_scores = []
    for group in groups:
        vals = []
        for gene in present_genes:
            val = to_float(matrix[gene].get(group, ""))
            if val is not None:
                vals.append(val)
        score = (sum(vals) / len(vals)) if vals else None
        group_scores.append(score)

    z_scores = zscore(group_scores)
    score_rows = []
    for group, score, z_val in zip(groups, group_scores, z_scores):
        score_rows.append(
            {
                "group": group,
                "n_signature_genes_total": len(signature_genes),
                "n_signature_genes_in_matrix": len(present_genes),
                "score_mean_expr_proxy": "" if score is None else f"{score:.8f}",
                "score_z_proxy": "" if z_val is None else f"{z_val:.8f}",
            }
        )

    note_rows = [
        {
            "panel": "Fig4G",
            "status": "proxy",
            "note": "Cell-level BTB scores used for violin were not found as numeric table in this workspace. Proxy score is mean expression of BTB signature genes computed from group-average ST matrix. OA group is excluded in this submission export to match Figure 4 layout.",
            "raw_support_files": "raw_support/Fig4G_BTB_signature_genes.csv;raw_support/Fig4_ST_average_expression_matrix.csv",
        }
    ]

    write_csv(
        out_dir / "Fig4G_BTB_signature_genes_used.csv",
        gene_used_rows,
        ["gene", "in_ST_average_matrix"],
    )
    write_csv(
        out_dir / "Fig4G_BTB_gene_expression_by_group.csv",
        gene_group_rows,
        ["gene", "group", "avg_expr"],
    )
    write_csv(
        out_dir / "Fig4G_BTB_group_score_proxy.csv",
        score_rows,
        [
            "group",
            "n_signature_genes_total",
            "n_signature_genes_in_matrix",
            "score_mean_expr_proxy",
            "score_z_proxy",
        ],
    )
    write_csv(
        out_dir / "Fig4G_proxy_note.csv",
        note_rows,
        ["panel", "status", "note", "raw_support_files"],
    )


def build_panel_hi_note(project_root: Path, out_dir: Path):
    note_rows = [
        {
            "panel": "Fig4H",
            "data_type": "representative_IF_images",
            "numeric_source_data": "not_applicable",
            "figure_file_reference": str(project_root / "Figure 4_20251210.ai"),
        },
        {
            "panel": "Fig4I",
            "data_type": "representative_IF_images",
            "numeric_source_data": "not_applicable",
            "figure_file_reference": str(project_root / "Figure 4_20251210.ai"),
        },
    ]
    write_csv(
        out_dir / "Fig4HI_image_panels_note.csv",
        note_rows,
        ["panel", "data_type", "numeric_source_data", "figure_file_reference"],
    )


def write_mapping_and_readme(project_root: Path, out_dir: Path, panel_f_mode: str, panel_g_mode: str):
    panel_f_is_raw = panel_f_mode == "raw"
    panel_f_generated = (
        "Fig4F_regulon_activity_matrix_wide.csv;Fig4F_regulon_activity_matrix_long.csv;Fig4F_regulon_candidates_from_scaled_top.csv;Fig4F_note.csv"
        if panel_f_is_raw
        else "Fig4F_regulon_activity_proxy_matrix_wide.csv;Fig4F_regulon_activity_proxy_matrix_long.csv;Fig4F_proxy_note.csv"
    )
    panel_f_raw_input = (
        "raw_support/Fig4F_pyscenic_regulonActivity_CellType.xls;raw_support/Fig4F_pyscenic_regulonActivity_byCellType_Scaled_top.xls;raw_support/Fig4F_pyscenic_heatmap_top_annotation.xls;raw_support/Fig4F_pyscenic_aucell.xls;raw_support/Fig4F_pyscenic_regulons.csv"
        if panel_f_is_raw
        else "raw_support/Fig4F_regulon_group_heatmap_source.pdf;raw_support/Fig4F_regulon_auc_group_source.pdf;raw_support/Fig4_ST_average_expression_matrix.csv"
    )
    panel_f_status = "available" if panel_f_is_raw else "proxy_only"
    panel_f_availability = "full" if panel_f_is_raw else "partial_proxy"
    panel_f_availability_note = (
        "Original numeric regulon AUC/Scaled matrices recovered from pyscenic outputs; OA group excluded to match Figure 4."
        if panel_f_is_raw
        else "Original numeric regulon AUC matrix not found; proxy generated from group-average TF expression, with OA group excluded to match Figure 4."
    )

    panel_g_is_raw = panel_g_mode == "raw"
    panel_g_generated = (
        "Fig4G_BTB_score_cells.csv;Fig4G_BTB_group_score_summary.csv;Fig4G_BTB_stage_score_summary.csv;Fig4G_note.csv"
        if panel_g_is_raw
        else "Fig4G_BTB_signature_genes_used.csv;Fig4G_BTB_gene_expression_by_group.csv;Fig4G_BTB_group_score_proxy.csv;Fig4G_proxy_note.csv"
    )
    panel_g_raw_input = (
        "raw_support/Fig4G_BTB_Integrity_Score_expression.csv"
        if panel_g_is_raw
        else "raw_support/Fig4G_BTB_signature_genes.csv;raw_support/Fig4G_BTB_core_for_heatmap.txt;raw_support/Fig4_ST_average_expression_matrix.csv"
    )
    panel_g_status = "available" if panel_g_is_raw else "proxy_only"
    panel_g_availability = "full" if panel_g_is_raw else "partial_proxy"
    panel_g_availability_note = (
        "Cell-level BTB integrity score table recovered from ST extraction; OA group excluded to match Figure 4."
        if panel_g_is_raw
        else "Original cell-level BTB score table not found; proxy generated from BTB signature + group-average ST expression, with OA group excluded to match Figure 4."
    )

    mapping_rows = [
        {
            "panel": "A",
            "subpanel": "ST UMAP / trajectory",
            "generated_source_data": "Fig4A_ST_cells_trajectory.csv;Fig4A_trajectory_line_points.csv;Fig4A_ST_counts_by_group.csv;Fig4A_ST_counts_by_state.csv;Fig4A_ST_pseudotime_summary.csv;Fig4A_stage_candidate_mapping.csv",
            "raw_input": "raw_support/Fig4A_monocle_pseudotime_cells.xls;raw_support/Fig4A_monocle_plot_coords.xls;raw_support/Fig4A_monocle_trajectory_lines.xls",
            "code": "../figure4_code/figure4_source_data.py",
            "status": "available",
        },
        {
            "panel": "B",
            "subpanel": "Stage heatmap matrix",
            "generated_source_data": "Fig4B_heatmap_matrix_wide.csv;Fig4B_heatmap_matrix_long.csv;Fig4B_heatmap_gene_order.csv",
            "raw_input": "raw_support/Fig4B_heatmap_export.xlsx;raw_support/Fig4B_ST_heatmap_export_bundle.zip",
            "code": "../figure4_code/figure4_source_data.py;../figure4_code/raw_original_scripts/Fig4_ST_heatmap_clustergvis_original.R",
            "status": "available",
        },
        {
            "panel": "C",
            "subpanel": "Stage marker violins",
            "generated_source_data": "Fig4C_stage_marker_cells.csv;Fig4C_stage_marker_summary.csv",
            "raw_input": "raw_support/Fig4_expression_markers_v1.csv;raw_support/Fig4_expression_markers_v2.csv;raw_support/Fig4A_monocle_pseudotime_cells.xls",
            "code": "../figure4_code/figure4_source_data.py",
            "status": "available",
        },
        {
            "panel": "D",
            "subpanel": "Stage composition by group",
            "generated_source_data": "Fig4D_stage_counts_by_group.csv;Fig4D_stage_ratio_by_group.csv;Fig4D_stage_ratio_all_group.csv;Fig4D_stage_counts_selected_samples.csv;Fig4D_stage_ratio_selected_samples.csv",
            "raw_input": "raw_support/Fig4D_cluster_counts_by_group_raw.csv;raw_support/Fig4D_view_stacked_bars_selected_samples.csv",
            "code": "../figure4_code/figure4_source_data.py;../figure4_code/raw_original_scripts/Fig4_cluster_ratio_original.R",
            "status": "available",
        },
        {
            "panel": "E",
            "subpanel": "Group marker violins",
            "generated_source_data": "Fig4E_group_marker_cells.csv;Fig4E_group_marker_summary.csv",
            "raw_input": "raw_support/Fig4_expression_markers_v1.csv;raw_support/Fig4_expression_markers_v2.csv;raw_support/Fig4A_monocle_pseudotime_cells.xls",
            "code": "../figure4_code/figure4_source_data.py",
            "status": "available",
        },
        {
            "panel": "F",
            "subpanel": "ST regulon activity heatmap",
            "generated_source_data": panel_f_generated,
            "raw_input": panel_f_raw_input,
            "code": "../figure4_code/figure4_source_data.py",
            "status": panel_f_status,
        },
        {
            "panel": "G",
            "subpanel": "BTB integrity score",
            "generated_source_data": panel_g_generated,
            "raw_input": panel_g_raw_input,
            "code": "../figure4_code/figure4_source_data.py",
            "status": panel_g_status,
        },
        {
            "panel": "H-I",
            "subpanel": "IF images",
            "generated_source_data": "Fig4HI_image_panels_note.csv",
            "raw_input": "Figure 4_20251210.ai;Figure 4_20251210.pdf",
            "code": "NA",
            "status": "image_only",
        },
    ]

    write_csv(
        out_dir / "Fig4_file_mapping.csv",
        mapping_rows,
        ["panel", "subpanel", "generated_source_data", "raw_input", "code", "status"],
    )

    availability_rows = [
        {
            "panel": "A",
            "availability": "full",
            "note": "Cell-level and trajectory line data available.",
        },
        {
            "panel": "B",
            "availability": "full",
            "note": "Heatmap matrix file available in ST heatmap export bundle.",
        },
        {
            "panel": "C",
            "availability": "full",
            "note": "Per-cell marker expression available from expression marker tables.",
        },
        {
            "panel": "D",
            "availability": "full",
            "note": "Stage composition counts available.",
        },
        {
            "panel": "E",
            "availability": "full",
            "note": "Per-cell marker expression by group available.",
        },
        {
            "panel": "F",
            "availability": panel_f_availability,
            "note": panel_f_availability_note,
        },
        {
            "panel": "G",
            "availability": panel_g_availability,
            "note": panel_g_availability_note,
        },
        {
            "panel": "H-I",
            "availability": "image_only",
            "note": "Representative IF image panels; numeric source data not applicable.",
        },
    ]
    write_csv(
        out_dir / "Fig4_panels_data_availability.csv",
        availability_rows,
        ["panel", "availability", "note"],
    )

    panel_f_readme_tables = (
        [
            "- Fig4F_regulon_activity_matrix_wide.csv",
            "- Fig4F_regulon_activity_matrix_long.csv",
            "- Fig4F_regulon_candidates_from_scaled_top.csv",
            "- Fig4F_note.csv",
        ]
        if panel_f_is_raw
        else [
            "- Fig4F_regulon_activity_proxy_matrix_wide.csv",
            "- Fig4F_regulon_activity_proxy_matrix_long.csv",
            "- Fig4F_proxy_note.csv",
        ]
    )
    panel_f_readme_note = (
        "- Panel F uses raw pyscenic regulon matrices extracted from ST细胞提取-scenic1_regulon.tar.gz, and excludes OA to match Figure 4."
        if panel_f_is_raw
        else "- Panel F is proxy because original numeric matrix files were not found in the current workspace; OA is excluded to match Figure 4."
    )
    panel_g_readme_tables = (
        [
            "- Fig4G_BTB_score_cells.csv",
            "- Fig4G_BTB_group_score_summary.csv",
            "- Fig4G_BTB_stage_score_summary.csv",
            "- Fig4G_note.csv",
        ]
        if panel_g_is_raw
        else [
            "- Fig4G_BTB_signature_genes_used.csv",
            "- Fig4G_BTB_gene_expression_by_group.csv",
            "- Fig4G_BTB_group_score_proxy.csv",
            "- Fig4G_proxy_note.csv",
        ]
    )
    panel_g_readme_note = (
        "- Panel G uses raw cell-level BTB integrity score table from ST细胞提取 - BTB_Integrity_Score_expression.csv, and excludes OA to match Figure 4."
        if panel_g_is_raw
        else "- Panel G is proxy because original cell-level BTB score table was not found; OA is excluded to match Figure 4."
    )

    readme_lines = [
        "Figure 4 source data (Panels A-I)",
        "",
        "Generated tables:",
        "- Fig4A_ST_cells_trajectory.csv",
        "- Fig4A_trajectory_line_points.csv",
        "- Fig4A_stage_candidate_mapping.csv",
        "- Fig4A_ST_counts_by_group.csv",
        "- Fig4A_ST_counts_by_state.csv",
        "- Fig4A_ST_pseudotime_summary.csv",
        "- Fig4B_heatmap_matrix_wide.csv",
        "- Fig4B_heatmap_matrix_long.csv",
        "- Fig4B_heatmap_gene_order.csv",
        "- Fig4C_stage_marker_cells.csv",
        "- Fig4C_stage_marker_summary.csv",
        "- Fig4D_stage_counts_by_group.csv",
        "- Fig4D_stage_ratio_by_group.csv",
        "- Fig4D_stage_ratio_all_group.csv",
        "- Fig4D_stage_counts_selected_samples.csv",
        "- Fig4D_stage_ratio_selected_samples.csv",
        "- Fig4E_group_marker_cells.csv",
        "- Fig4E_group_marker_summary.csv",
    ]
    readme_lines.extend(panel_f_readme_tables)
    readme_lines.extend(panel_g_readme_tables)
    readme_lines.extend(
        [
        "- Fig4HI_image_panels_note.csv",
        "- Fig4_file_mapping.csv",
        "- Fig4_panels_data_availability.csv",
        "",
        "Notes:",
        "- ST stage labels in this export use candidate mapping ST1->Stage_a, ST2->Stage_b, ST3->Stage_c (see Fig4A_stage_candidate_mapping.csv).",
        panel_f_readme_note,
        panel_g_readme_note,
        "",
        "Build script:",
        "- ../figure4_code/figure4_source_data.py",
        ]
    )
    (out_dir / "README.txt").write_text("\n".join(readme_lines), encoding="utf-8")

    code_readme_lines = [
        "Figure 4 code files",
        "",
        "1) figure4_source_data.py",
        "- Build Figure 4 source-data CSV tables from local files under:",
        "  ../figure4_source_data/raw_support/",
        "- No external Python dependencies are required.",
        "",
        "2) raw_original_scripts/",
        "- Original scripts recovered from workspace:",
        "  - Fig4_cluster_ratio_original.R",
        "  - Fig4_ST_heatmap_clustergvis_original.R",
        "  - Fig4_ST_heatmap_clustergvis_alt_original.R",
        "",
        "Run:",
        "- python3 figure4_source_data.py",
    ]
    (project_root / "source data" / "figure4_code" / "README.txt").write_text(
        "\n".join(code_readme_lines),
        encoding="utf-8",
    )


def main():
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent
    out_dir = project_root / "source data" / "figure4_source_data"
    raw_support = out_dir / "raw_support"

    build_panel_a(raw_support, out_dir)
    build_panel_b(raw_support, out_dir)
    build_panel_c(raw_support, out_dir)
    build_panel_d(raw_support, out_dir)
    build_panel_e(raw_support, out_dir)
    panel_f_mode = build_panel_f(raw_support, out_dir)
    panel_g_mode = build_panel_g(raw_support, out_dir)
    build_panel_hi_note(project_root, out_dir)
    write_mapping_and_readme(project_root, out_dir, panel_f_mode, panel_g_mode)

    print(f"Wrote Figure 4 source data to: {out_dir}")


if __name__ == "__main__":
    main()
