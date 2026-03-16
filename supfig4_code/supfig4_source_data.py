#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import re
import shutil
import statistics
import tarfile
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
CELL_REF_RE = re.compile(r"^([A-Z]+)([0-9]+)$")

STAGE_MAP = {"ST1": "Stage_a", "ST2": "Stage_b", "ST3": "Stage_c"}
GROUP_ORDER = {"AZFc_Del": 0, "iNOA_B": 1, "iNOA_S": 2, "KS": 3}
TARGET_REGULONS = [
    "NFE2L1",
    "JUNB",
    "ATF6",
    "CEBPD",
    "FOS",
    "FOSB",
    "ATF3",
    "SETDB1",
    "RFX3",
    "ELF1",
    "ETV5",
    "SMAD1",
    "ELK4",
]

STALE_OUTPUT_FILES = [
    "SupFig4D_ST1vsST3_gsea_all.csv",
    "SupFig4D_ST1vsST3_top20.csv",
    "SupFig4F_ST3vsST1_gsea_all.csv",
    "SupFig4F_ST3vsST1_top20.csv",
]


def to_float(value):
    if value is None:
        return None
    text = str(value).strip().strip('"')
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
    if upper == "OA":
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


def read_tsv_rows(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv_rows(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


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


def remove_stale_outputs(out_dir: Path):
    for filename in STALE_OUTPUT_FILES:
        stale_path = out_dir / filename
        if stale_path.exists():
            stale_path.unlink()


def write_source_data_workbook(out_dir: Path):
    try:
        from openpyxl import Workbook
    except Exception:
        return ""

    workbook_spec = [
        ("A_cells", "SupFig4A_ST_cells_pseudotime.csv"),
        ("A_markers", "SupFig4A_stage_markers_from_diffgene.csv"),
        ("B_target", "SupFig4B_regulon_stage_matrix_target.csv"),
        ("B_detected", "SupFig4B_regulon_stage_matrix_detected.csv"),
        ("C_cells", "SupFig4C_glycolysis_score_cells.csv"),
        ("C_summary", "SupFig4C_glycolysis_score_summary.csv"),
        ("C_gene_set", "SupFig4C_gene_set_candidates.csv"),
        ("D_all", "SupFig4D_ST1_disease_vs_ctrl_gsea_all.csv"),
        ("D_top10", "SupFig4D_ST1_disease_vs_ctrl_top10.csv"),
        ("E_cells", "SupFig4E_inflammatory_score_cells.csv"),
        ("E_summary", "SupFig4E_inflammatory_score_summary.csv"),
        ("E_gene_set", "SupFig4E_gene_set_candidates.csv"),
        ("F_all", "SupFig4F_ST3_disease_vs_ctrl_gsea_all.csv"),
        ("F_top10", "SupFig4F_ST3_disease_vs_ctrl_top10.csv"),
        ("missing_note", "SupFig4C_E_missing_cell_score_raw.csv"),
        ("file_mapping", "SupFig4_file_mapping.csv"),
        ("availability", "SupFig4_panels_data_availability.csv"),
    ]

    wb = Workbook()
    wb.remove(wb.active)

    for sheet_name, csv_name in workbook_spec:
        path = out_dir / csv_name
        if not path.exists():
            continue
        ws = wb.create_sheet(title=sheet_name)
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.reader(handle):
                ws.append(row)

    out_path = out_dir / "SourceData_SupFig4.xlsx"
    wb.save(out_path)
    return out_path.name


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

        sheet_root = ET.fromstring(xlsx.read("xl/worksheets/sheet1.xml"))
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
            rows.append([row_map.get(idx, "") for idx in range(max_col + 1)] if max_col >= 0 else [])
        return rows


def copy_panel_raw_to_support(panel_raw_dir: Path, raw_support_dir: Path):
    raw_support_dir.mkdir(parents=True, exist_ok=True)
    copied = []
    for file_path in sorted(panel_raw_dir.glob("panel_*/raw_data/*")):
        if file_path.is_file():
            target = raw_support_dir / file_path.name
            shutil.copy2(file_path, target)
            copied.append(target.name)
    return copied


def copy_external_st_gsea_raw(project_root: Path, raw_support_dir: Path):
    external_dir = project_root.parent / "睾丸文章数据" / "ST" / "疾病的GSEA"
    file_map = {
        "ST1-diseasevsctrl.csv": "SupFig4D_ST1-diseasevsctrl.csv",
        "ST3-diseasevsctrl.csv": "SupFig4F_ST3-diseasevsctrl.csv",
        "ST1_疾病vsctrl-NES-10.pdf": "SupFig4D_ST1_疾病vsctrl-NES-10.pdf",
        "ST3_疾病vsctrl-NES-10.pdf": "SupFig4F_ST3_疾病vsctrl-NES-10.pdf",
        "气泡图统一标尺.R": "SupFig4DF_气泡图统一标尺.R",
    }
    copied = []
    missing = []
    for source_name, target_name in file_map.items():
        source_path = external_dir / source_name
        target_path = raw_support_dir / target_name
        if source_path.is_file():
            shutil.copy2(source_path, target_path)
            copied.append(target_name)
        else:
            missing.append(str(source_path))
    if missing:
        raise FileNotFoundError(
            "Missing expected ST disease GSEA raw files:\n- " + "\n- ".join(missing)
        )
    return copied


def build_panel_a(out_dir: Path, raw_support_dir: Path):
    pseudo_path = raw_support_dir / "SupFig4A_out.monocle_Pseudotime.xls"
    coord_path = raw_support_dir / "SupFig4A_out.plotpseudotimedata.xls"

    pseudo_rows = read_tsv_rows(pseudo_path)
    coord_rows = read_tsv_rows(coord_path)
    coord_map = {}
    for row in coord_rows:
        cell_id = row.get("cellID", "").strip()
        if cell_id:
            coord_map[cell_id] = (row.get("X1", ""), row.get("X2", ""))

    cell_rows = []
    for row in pseudo_rows:
        cell_id = row.get("barcode", "").strip()
        cluster = row.get("celltype", "").strip()
        stage = STAGE_MAP.get(cluster, "")
        group = normalize_group(row.get("group", ""))
        traj_x, traj_y = coord_map.get(cell_id, ("", ""))
        cell_rows.append(
            {
                "cell_id": cell_id,
                "sample_id": row.get("sample", ""),
                "group": group,
                "st_cluster": cluster,
                "stage_label": stage,
                "state": row.get("State", ""),
                "pseudotime": row.get("Pseudotime", ""),
                "traj_x": traj_x,
                "traj_y": traj_y,
                "nCount_RNA": row.get("nCount_RNA", ""),
                "nFeature_RNA": row.get("nFeature_RNA", ""),
            }
        )

    write_csv(
        out_dir / "SupFig4A_ST_cells_pseudotime.csv",
        cell_rows,
        [
            "cell_id",
            "sample_id",
            "group",
            "st_cluster",
            "stage_label",
            "state",
            "pseudotime",
            "traj_x",
            "traj_y",
            "nCount_RNA",
            "nFeature_RNA",
        ],
    )

    diff_files = [
        ("SupFig4A_clusterImmatureSTs_1_diffgenes.xls", "ST1"),
        ("SupFig4A_clusterImmatureSTs_2_diffgenes.xls", "ST2"),
        ("SupFig4A_clusterMatureSTs_diffgenes.xls", "ST3"),
    ]
    marker_rows = []
    for filename, cluster in diff_files:
        stage = STAGE_MAP[cluster]
        rows = read_tsv_rows(raw_support_dir / filename)
        for idx, row in enumerate(rows, start=1):
            marker_rows.append(
                {
                    "st_cluster": cluster,
                    "stage_label": stage,
                    "rank_in_file": idx,
                    "gene": row.get("names", ""),
                    "scores": row.get("scores", ""),
                    "logfoldchanges": row.get("logfoldchanges", ""),
                    "pvals": row.get("pvals", ""),
                    "pvals_adj": row.get("pvals_adj", ""),
                    "pct_nz_group": row.get("pct_nz_group", ""),
                    "pct_nz_reference": row.get("pct_nz_reference", ""),
                    "source_file": filename,
                }
            )

    write_csv(
        out_dir / "SupFig4A_stage_markers_from_diffgene.csv",
        marker_rows,
        [
            "st_cluster",
            "stage_label",
            "rank_in_file",
            "gene",
            "scores",
            "logfoldchanges",
            "pvals",
            "pvals_adj",
            "pct_nz_group",
            "pct_nz_reference",
            "source_file",
        ],
    )


def stage_map_from_panel_a(path: Path):
    rows = read_csv_rows(path)
    result = {}
    for row in rows:
        cell_id = row.get("cell_id", "")
        stage = row.get("stage_label", "")
        if cell_id and stage:
            result[cell_id] = stage
    return result


def extract_scenic_from_tar(raw_support_dir: Path):
    tar_path = raw_support_dir / "SupFig4B_ST细胞提取-scenic1_regulon.tar.gz"
    targets = {
        "pyscenic_aucell.xls",
        "pyscenic_heatmap_top.xls",
        "pyscenic_heatmap_top_annotation.xls",
        "pyscenic_regulonActivity_CellType.xls",
        "pyscenic_regulonActivity_byCellType_Scaled_top.xls",
        "pyscenic_regulons.csv",
    }
    extracted = {}
    with tarfile.open(tar_path, "r:gz") as tar:
        for member in tar.getmembers():
            name = Path(member.name).name
            if name in targets:
                out_path = raw_support_dir / f"SupFig4B_{name}"
                with tar.extractfile(member) as src, out_path.open("wb") as dst:
                    if src is not None:
                        shutil.copyfileobj(src, dst)
                extracted[name] = out_path
    return extracted


def parse_aucell_stage_means(aucell_path: Path, stage_map: dict[str, str]):
    stage_sums = defaultdict(lambda: {"Stage_a": 0.0, "Stage_b": 0.0, "Stage_c": 0.0})
    stage_counts = defaultdict(lambda: {"Stage_a": 0, "Stage_b": 0, "Stage_c": 0})
    labels_by_gene = defaultdict(list)

    with aucell_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            return [], [], set()
        cell_ids = [cell.strip() for cell in header[1:]]
        cell_stages = [stage_map.get(cell, "") for cell in cell_ids]

        for row in reader:
            if not row:
                continue
            label = str(row[0]).strip().strip('"')
            if not label:
                continue
            gene = label.split("(", 1)[0].strip()
            labels_by_gene[gene].append(label)
            values = row[1:]
            for idx, raw_value in enumerate(values):
                if idx >= len(cell_stages):
                    break
                stage = cell_stages[idx]
                if stage not in {"Stage_a", "Stage_b", "Stage_c"}:
                    continue
                value = to_float(raw_value)
                if value is None:
                    continue
                stage_sums[label][stage] += value
                stage_counts[label][stage] += 1

    detected_rows = []
    target_rows = []
    available_genes = set()

    for label in sorted(stage_sums.keys()):
        gene = label.split("(", 1)[0].strip()
        means = {}
        for stage in ["Stage_a", "Stage_b", "Stage_c"]:
            count = stage_counts[label][stage]
            means[stage] = stage_sums[label][stage] / count if count else None
        detected_rows.append(
            {
                "regulon_label": label,
                "gene": gene,
                "Stage_a": "" if means["Stage_a"] is None else f"{means['Stage_a']:.8f}",
                "Stage_b": "" if means["Stage_b"] is None else f"{means['Stage_b']:.8f}",
                "Stage_c": "" if means["Stage_c"] is None else f"{means['Stage_c']:.8f}",
                "n_stage_a": stage_counts[label]["Stage_a"],
                "n_stage_b": stage_counts[label]["Stage_b"],
                "n_stage_c": stage_counts[label]["Stage_c"],
            }
        )

    label_priority = {}
    for gene, labels in labels_by_gene.items():
        label_priority[gene] = sorted(labels)[0]

    for gene in TARGET_REGULONS:
        label = label_priority.get(gene, "")
        if label and label in stage_sums:
            available_genes.add(gene)
            means = {}
            for stage in ["Stage_a", "Stage_b", "Stage_c"]:
                count = stage_counts[label][stage]
                means[stage] = stage_sums[label][stage] / count if count else None
            target_rows.append(
                {
                    "target_gene": gene,
                    "regulon_label": label,
                    "available_in_aucell": "1",
                    "Stage_a": "" if means["Stage_a"] is None else f"{means['Stage_a']:.8f}",
                    "Stage_b": "" if means["Stage_b"] is None else f"{means['Stage_b']:.8f}",
                    "Stage_c": "" if means["Stage_c"] is None else f"{means['Stage_c']:.8f}",
                    "n_stage_a": stage_counts[label]["Stage_a"],
                    "n_stage_b": stage_counts[label]["Stage_b"],
                    "n_stage_c": stage_counts[label]["Stage_c"],
                }
            )
        else:
            target_rows.append(
                {
                    "target_gene": gene,
                    "regulon_label": "",
                    "available_in_aucell": "0",
                    "Stage_a": "",
                    "Stage_b": "",
                    "Stage_c": "",
                    "n_stage_a": "",
                    "n_stage_b": "",
                    "n_stage_c": "",
                }
            )

    return detected_rows, target_rows, available_genes


def build_panel_b(out_dir: Path, raw_support_dir: Path):
    extracted = extract_scenic_from_tar(raw_support_dir)
    stage_map = stage_map_from_panel_a(out_dir / "SupFig4A_ST_cells_pseudotime.csv")

    aucell_path = extracted.get("pyscenic_aucell.xls")
    if aucell_path is None:
        write_csv(
            out_dir / "SupFig4B_regulon_stage_matrix_target.csv",
            [],
            ["target_gene", "regulon_label", "available_in_aucell", "Stage_a", "Stage_b", "Stage_c", "n_stage_a", "n_stage_b", "n_stage_c"],
        )
        write_csv(
            out_dir / "SupFig4B_regulon_stage_matrix_detected.csv",
            [],
            ["regulon_label", "gene", "Stage_a", "Stage_b", "Stage_c", "n_stage_a", "n_stage_b", "n_stage_c"],
        )
        return set()

    detected_rows, target_rows, available_genes = parse_aucell_stage_means(aucell_path, stage_map)
    write_csv(
        out_dir / "SupFig4B_regulon_stage_matrix_detected.csv",
        detected_rows,
        ["regulon_label", "gene", "Stage_a", "Stage_b", "Stage_c", "n_stage_a", "n_stage_b", "n_stage_c"],
    )
    write_csv(
        out_dir / "SupFig4B_regulon_stage_matrix_target.csv",
        target_rows,
        ["target_gene", "regulon_label", "available_in_aucell", "Stage_a", "Stage_b", "Stage_c", "n_stage_a", "n_stage_b", "n_stage_c"],
    )
    return available_genes


def read_single_column_genes(path: Path):
    genes = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        first = True
        for row in reader:
            if not row:
                continue
            value = str(row[0]).strip()
            if value == "":
                continue
            if first:
                first = False
                if value.lower() in {"gene", "genes", "sasp", "reactome_cytokine_signaling_in_immune_system", "canonical glycolysis"}:
                    continue
            genes.append(value)
    return genes


def build_panel_c_e(out_dir: Path, raw_support_dir: Path):
    glycolysis_genes = read_single_column_genes(raw_support_dir / "SupFig4C_canonical_glycolysis.csv")
    immune_genes = read_single_column_genes(raw_support_dir / "SupFig4C_MMUNE_SYSTEM_candidate.csv")
    immune_genes_e = read_single_column_genes(raw_support_dir / "SupFig4E_MMUNE_SYSTEM.csv")
    sasp_genes = read_single_column_genes(raw_support_dir / "SupFig4E_SASP.csv")
    glycolysis_score_path = raw_support_dir / "SupFig4C_ST_Gene_set_enrichment1_expression.csv"
    inflammatory_score_path = raw_support_dir / "SupFig4E_ST_Gene_set_enrichment7_expression.csv"

    panel_c_rows = []
    for gene in glycolysis_genes:
        panel_c_rows.append({"panel": "C", "gene_set": "canonical_glycolysis", "gene": gene, "source_file": "SupFig4C_canonical_glycolysis.csv"})
    for gene in immune_genes:
        panel_c_rows.append({"panel": "C", "gene_set": "MMUNE_SYSTEM_candidate", "gene": gene, "source_file": "SupFig4C_MMUNE_SYSTEM_candidate.csv"})

    panel_e_rows = []
    for gene in sasp_genes:
        panel_e_rows.append({"panel": "E", "gene_set": "SASP", "gene": gene, "source_file": "SupFig4E_SASP.csv"})
    for gene in immune_genes_e:
        panel_e_rows.append({"panel": "E", "gene_set": "MMUNE_SYSTEM", "gene": gene, "source_file": "SupFig4E_MMUNE_SYSTEM.csv"})

    write_csv(
        out_dir / "SupFig4C_gene_set_candidates.csv",
        panel_c_rows,
        ["panel", "gene_set", "gene", "source_file"],
    )
    write_csv(
        out_dir / "SupFig4E_gene_set_candidates.csv",
        panel_e_rows,
        ["panel", "gene_set", "gene", "source_file"],
    )

    missing_rows = []
    if glycolysis_score_path.exists():
        score_rows = []
        score_by_stage = defaultdict(list)
        with glycolysis_score_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                stage_raw = str(row.get("Major cell types", "")).strip()
                stage_label = stage_raw
                if stage_raw == "ST_a":
                    stage_label = "Stage_a"
                elif stage_raw == "ST_b":
                    stage_label = "Stage_b"
                elif stage_raw == "ST_c":
                    stage_label = "Stage_c"
                score = to_float(row.get("Lactate_Glycolysis_Signature.csv", ""))
                if score is None:
                    continue
                score_by_stage[stage_label].append(score)
                score_rows.append(
                    {
                        "cell_id": row.get("cell", ""),
                        "glycolysis_score": f"{score:.8f}",
                        "st_stage_raw": stage_raw,
                        "stage_label": stage_label,
                        "umap1": row.get("UMAP1", ""),
                        "umap2": row.get("UMAP2", ""),
                        "tsne1": row.get("TSNE1", ""),
                        "tsne2": row.get("TSNE2", ""),
                    }
                )
        write_csv(
            out_dir / "SupFig4C_glycolysis_score_cells.csv",
            score_rows,
            ["cell_id", "glycolysis_score", "st_stage_raw", "stage_label", "umap1", "umap2", "tsne1", "tsne2"],
        )

        summary_rows = []
        ordered_stages = ["Stage_a", "Stage_b", "Stage_c"]
        for stage in ordered_stages:
            values = score_by_stage.get(stage, [])
            if not values:
                continue
            summary_rows.append(
                {
                    "stage_label": stage,
                    "n_cells": len(values),
                    "mean_score": f"{statistics.mean(values):.8f}",
                    "median_score": f"{statistics.median(values):.8f}",
                    "min_score": f"{min(values):.8f}",
                    "max_score": f"{max(values):.8f}",
                }
            )
        write_csv(
            out_dir / "SupFig4C_glycolysis_score_summary.csv",
            summary_rows,
            ["stage_label", "n_cells", "mean_score", "median_score", "min_score", "max_score"],
        )
        missing_rows.append(
            {
                "panel": "C",
                "status": "available",
                "note": "已定位并导出按ST_a/ST_b/ST_c分组的cell-level glycolysis score原始表。",
            }
        )
    else:
        missing_rows.append(
            {
                "panel": "C",
                "status": "missing_cell_level_score_table",
                "note": "未定位到按Stage_a/Stage_b/Stage_c分组的glycolytic score细胞级原始表。",
            }
        )

    if inflammatory_score_path.exists():
        score_rows = []
        score_by_stage = defaultdict(list)
        with inflammatory_score_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            meta_columns = {"cell", "Major cell types", "UMAP1", "UMAP2", "TSNE1", "TSNE2"}
            score_columns = [column for column in fieldnames if column not in meta_columns]
            score_column = score_columns[0] if score_columns else ""

            for row in reader:
                stage_raw = str(row.get("Major cell types", "")).strip()
                stage_label = stage_raw
                if stage_raw == "ST_a":
                    stage_label = "Stage_a"
                elif stage_raw == "ST_b":
                    stage_label = "Stage_b"
                elif stage_raw == "ST_c":
                    stage_label = "Stage_c"
                score = to_float(row.get(score_column, ""))
                if score is None:
                    continue
                score_by_stage[stage_label].append(score)
                score_rows.append(
                    {
                        "cell_id": row.get("cell", ""),
                        "inflammatory_score": f"{score:.8f}",
                        "score_column": score_column,
                        "st_stage_raw": stage_raw,
                        "stage_label": stage_label,
                        "umap1": row.get("UMAP1", ""),
                        "umap2": row.get("UMAP2", ""),
                        "tsne1": row.get("TSNE1", ""),
                        "tsne2": row.get("TSNE2", ""),
                    }
                )

        write_csv(
            out_dir / "SupFig4E_inflammatory_score_cells.csv",
            score_rows,
            ["cell_id", "inflammatory_score", "score_column", "st_stage_raw", "stage_label", "umap1", "umap2", "tsne1", "tsne2"],
        )

        summary_rows = []
        ordered_stages = ["Stage_a", "Stage_b", "Stage_c"]
        score_column_value = score_rows[0]["score_column"] if score_rows else ""
        for stage in ordered_stages:
            values = score_by_stage.get(stage, [])
            if not values:
                continue
            summary_rows.append(
                {
                    "stage_label": stage,
                    "score_column": score_column_value,
                    "n_cells": len(values),
                    "mean_score": f"{statistics.mean(values):.8f}",
                    "median_score": f"{statistics.median(values):.8f}",
                    "min_score": f"{min(values):.8f}",
                    "max_score": f"{max(values):.8f}",
                }
            )
        write_csv(
            out_dir / "SupFig4E_inflammatory_score_summary.csv",
            summary_rows,
            ["stage_label", "score_column", "n_cells", "mean_score", "median_score", "min_score", "max_score"],
        )
        missing_rows.append(
            {
                "panel": "E",
                "status": "available",
                "note": "已定位并导出按ST_a/ST_b/ST_c分组的cell-level inflammatory/immune score原始表。",
            }
        )
    else:
        missing_rows.append(
            {
                "panel": "E",
                "status": "missing_cell_level_score_table",
                "note": "未定位到按Stage_a/Stage_b/Stage_c分组的inflammatory response score细胞级原始表。",
            }
        )
    write_csv(out_dir / "SupFig4C_E_missing_cell_score_raw.csv", missing_rows, ["panel", "status", "note"])


def humanize_pathway(description: str):
    text = str(description).strip()
    if text.startswith("KEGG_"):
        text = text[5:]
    return text.replace("_", " ").title()


def parse_disease_gsea_csv(path: Path, panel_stage: str):
    rows = read_csv_rows(path)
    out_rows = []
    for row in rows:
        desc = str(row.get("Description", "")).strip()
        if desc == "":
            continue
        group = normalize_group(row.get("group", ""))
        nes = to_float(row.get("NES"))
        pvalue = to_float(row.get("pvalue"))
        padj = to_float(row.get("p.adjust"))
        out_rows.append(
            {
                "panel_stage": panel_stage,
                "group": group,
                "Description": desc,
                "pathway_readable": humanize_pathway(desc),
                "setSize": row.get("setSize", ""),
                "enrichmentScore": row.get("enrichmentScore", ""),
                "NES": row.get("NES", ""),
                "pvalue": row.get("pvalue", ""),
                "p.adjust": row.get("p.adjust", ""),
                "qvalues": row.get("qvalues", ""),
                "rank": row.get("rank", ""),
                "leading_edge": row.get("leading_edge", ""),
                "neglog10_fdr": "" if padj is None else f"{min(-math.log10(max(padj, 1e-300)), 20.0):.8f}",
                "_pvalue": pvalue,
                "_abs_nes": abs(nes) if nes is not None else None,
            }
        )
    return out_rows


def pick_top_by_group(rows, top_n: int):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["group"]].append(row)
    out = []
    for group in sorted(grouped, key=lambda value: GROUP_ORDER.get(value, 99)):
        ordered = sorted(
            grouped[group],
            key=lambda row: (
                row["_pvalue"] if row["_pvalue"] is not None else 999.0,
                -(row["_abs_nes"] if row["_abs_nes"] is not None else 0.0),
                row["pathway_readable"],
            ),
        )
        for idx, row in enumerate(ordered[:top_n], start=1):
            out.append(
                {
                    "rank_in_group": idx,
                    "panel_stage": row["panel_stage"],
                    "group": row["group"],
                    "Description": row["Description"],
                    "pathway_readable": row["pathway_readable"],
                    "NES": row["NES"],
                    "pvalue": row["pvalue"],
                    "p.adjust": row["p.adjust"],
                    "neglog10_fdr": row["neglog10_fdr"],
                    "setSize": row["setSize"],
                    "enrichmentScore": row["enrichmentScore"],
                    "rank": row["rank"],
                }
            )
    return out


def sort_disease_rows(rows):
    return sorted(
        rows,
        key=lambda row: (
            GROUP_ORDER.get(row["group"], 99),
            row["_pvalue"] if row["_pvalue"] is not None else 999.0,
            -(row["_abs_nes"] if row["_abs_nes"] is not None else 0.0),
            row["pathway_readable"],
        ),
    )


def build_panel_d_f(out_dir: Path, raw_support_dir: Path):
    d_rows = parse_disease_gsea_csv(raw_support_dir / "SupFig4D_ST1-diseasevsctrl.csv", "ST1")
    f_rows = parse_disease_gsea_csv(raw_support_dir / "SupFig4F_ST3-diseasevsctrl.csv", "ST3")
    d_rows_sorted = sort_disease_rows(d_rows)
    f_rows_sorted = sort_disease_rows(f_rows)
    d_top = pick_top_by_group(d_rows_sorted, 10)
    f_top = pick_top_by_group(f_rows_sorted, 10)

    all_cols = [
        "panel_stage",
        "group",
        "Description",
        "pathway_readable",
        "setSize",
        "enrichmentScore",
        "NES",
        "pvalue",
        "p.adjust",
        "qvalues",
        "rank",
        "leading_edge",
        "neglog10_fdr",
    ]
    top_cols = [
        "rank_in_group",
        "panel_stage",
        "group",
        "Description",
        "pathway_readable",
        "NES",
        "pvalue",
        "p.adjust",
        "neglog10_fdr",
        "setSize",
        "enrichmentScore",
        "rank",
    ]

    write_csv(out_dir / "SupFig4D_ST1_disease_vs_ctrl_gsea_all.csv", d_rows_sorted, all_cols)
    write_csv(out_dir / "SupFig4F_ST3_disease_vs_ctrl_gsea_all.csv", f_rows_sorted, all_cols)
    write_csv(out_dir / "SupFig4D_ST1_disease_vs_ctrl_top10.csv", d_top, top_cols)
    write_csv(out_dir / "SupFig4F_ST3_disease_vs_ctrl_top10.csv", f_top, top_cols)


def write_mapping_and_readme(
    out_dir: Path,
    available_regulon_genes: set[str],
    copied_raw_files: list[str],
    workbook_name: str,
):
    mapping_rows = [
        {
            "panel": "A",
            "subpanel": "Pseudotime heatmap input",
            "generated_source_data": "SupFig4A_ST_cells_pseudotime.csv;SupFig4A_stage_markers_from_diffgene.csv",
            "raw_input": "raw_support/SupFig4A_out.monocle_Pseudotime.xls;raw_support/SupFig4A_out.plotpseudotimedata.xls;raw_support/SupFig4A_clusterImmatureSTs_1_diffgenes.xls;raw_support/SupFig4A_clusterImmatureSTs_2_diffgenes.xls;raw_support/SupFig4A_clusterMatureSTs_diffgenes.xls",
            "code": "../supfig4_code/supfig4_source_data.py",
            "status": "available",
        },
        {
            "panel": "B",
            "subpanel": "Regulon x Stage heatmap",
            "generated_source_data": "SupFig4B_regulon_stage_matrix_target.csv;SupFig4B_regulon_stage_matrix_detected.csv",
            "raw_input": "raw_support/SupFig4B_ST细胞提取-scenic1_regulon.tar.gz;raw_support/SupFig4B_out.monocle_Pseudotime.xls;raw_support/SupFig4B_pyscenic_aucell.xls",
            "code": "../supfig4_code/supfig4_source_data.py;../supfig4_code/supfig4_panelB_heatmap_plot.py",
            "status": "partial",
        },
        {
            "panel": "C",
            "subpanel": "Glycolytic activity score violin",
            "generated_source_data": "SupFig4C_glycolysis_score_cells.csv;SupFig4C_glycolysis_score_summary.csv;SupFig4C_gene_set_candidates.csv;SupFig4C_E_missing_cell_score_raw.csv",
            "raw_input": "raw_support/SupFig4C_ST_Gene_set_enrichment1_expression.csv;raw_support/SupFig4C_canonical_glycolysis.csv;raw_support/SupFig4C_MMUNE_SYSTEM_candidate.csv",
            "code": "../supfig4_code/supfig4_source_data.py",
            "status": "available",
        },
        {
            "panel": "D",
            "subpanel": "ST1 enriched pathways bubble",
            "generated_source_data": "SupFig4D_ST1_disease_vs_ctrl_gsea_all.csv;SupFig4D_ST1_disease_vs_ctrl_top10.csv",
            "raw_input": "raw_support/SupFig4D_ST1-diseasevsctrl.csv;raw_support/SupFig4D_ST1_疾病vsctrl-NES-10.pdf;raw_support/SupFig4DF_气泡图统一标尺.R",
            "code": "../supfig4_code/supfig4_source_data.py;../supfig4_code/supfig4_panelDF_bubble_proxy_plot.py;../supfig4_code/supfig4_panelDF_bubble_original.R",
            "status": "available",
        },
        {
            "panel": "E",
            "subpanel": "Inflammatory response score violin",
            "generated_source_data": "SupFig4E_inflammatory_score_cells.csv;SupFig4E_inflammatory_score_summary.csv;SupFig4E_gene_set_candidates.csv;SupFig4C_E_missing_cell_score_raw.csv",
            "raw_input": "raw_support/SupFig4E_ST_Gene_set_enrichment7_expression.csv;raw_support/SupFig4E_SASP.csv;raw_support/SupFig4E_MMUNE_SYSTEM.csv",
            "code": "../supfig4_code/supfig4_source_data.py",
            "status": "available",
        },
        {
            "panel": "F",
            "subpanel": "ST3 enriched pathways bubble",
            "generated_source_data": "SupFig4F_ST3_disease_vs_ctrl_gsea_all.csv;SupFig4F_ST3_disease_vs_ctrl_top10.csv",
            "raw_input": "raw_support/SupFig4F_ST3-diseasevsctrl.csv;raw_support/SupFig4F_ST3_疾病vsctrl-NES-10.pdf;raw_support/SupFig4DF_气泡图统一标尺.R",
            "code": "../supfig4_code/supfig4_source_data.py;../supfig4_code/supfig4_panelDF_bubble_proxy_plot.py;../supfig4_code/supfig4_panelDF_bubble_original.R",
            "status": "available",
        },
    ]
    write_csv(
        out_dir / "SupFig4_file_mapping.csv",
        mapping_rows,
        ["panel", "subpanel", "generated_source_data", "raw_input", "code", "status"],
    )

    availability_rows = [
        {"panel": "A", "status": "available", "note": "可导出细胞拟时序与阶段差异基因列表。"},
        {
            "panel": "B",
            "status": "partial",
            "note": f"目标regulon中在aucell检测到 {len(available_regulon_genes)}/{len(TARGET_REGULONS)} 个：{';'.join(sorted(available_regulon_genes)) if available_regulon_genes else 'none'}。",
        },
        {"panel": "C", "status": "available", "note": "已定位并导出cell-level glycolysis score表（ST_a/ST_b/ST_c）。"},
        {"panel": "D", "status": "available", "note": "已定位并导出ST1-diseasevsctrl.csv（疾病vs对照，按group通路富集气泡图原始表）。"},
        {"panel": "E", "status": "available", "note": "已定位并导出cell-level inflammatory/immune score表（ST_a/ST_b/ST_c）。"},
        {"panel": "F", "status": "available", "note": "已定位并导出ST3-diseasevsctrl.csv（疾病vs对照，按group通路富集气泡图原始表）。"},
    ]
    write_csv(out_dir / "SupFig4_panels_data_availability.csv", availability_rows, ["panel", "status", "note"])

    readme_lines = [
        "Supplementary Figure 4 source data (Panels A-F)",
        "",
        "Generated files:",
        "- SupFig4A_ST_cells_pseudotime.csv",
        "- SupFig4A_stage_markers_from_diffgene.csv",
        "- SupFig4B_regulon_stage_matrix_target.csv",
        "- SupFig4B_regulon_stage_matrix_detected.csv",
        "- SupFig4C_gene_set_candidates.csv",
        "- SupFig4C_glycolysis_score_cells.csv",
        "- SupFig4C_glycolysis_score_summary.csv",
        "- SupFig4D_ST1_disease_vs_ctrl_gsea_all.csv",
        "- SupFig4D_ST1_disease_vs_ctrl_top10.csv",
        "- SupFig4E_gene_set_candidates.csv",
        "- SupFig4E_inflammatory_score_cells.csv",
        "- SupFig4E_inflammatory_score_summary.csv",
        "- SupFig4F_ST3_disease_vs_ctrl_gsea_all.csv",
        "- SupFig4F_ST3_disease_vs_ctrl_top10.csv",
        "- SupFig4C_E_missing_cell_score_raw.csv",
        "- SupFig4_file_mapping.csv",
        "- SupFig4_panels_data_availability.csv",
    ]
    if workbook_name:
        readme_lines.append(f"- {workbook_name}")
    readme_lines.extend(
        [
            "",
            "Raw support copied from panel_raw:",
            f"- {len(copied_raw_files)} files in ./raw_support/",
            "",
            "Notes:",
            "- Panel B target regulons follow Supplementary Figure 4 label set and are joined to ST stage via monocle cell IDs.",
            "- Panel C glycolysis score table has been located and exported at cell level.",
            "- Panel E inflammatory/immune score table has been located and exported at cell level.",
            "- Panel D/F now use original ST disease-vs-control KEGG GSEA tables (ST1-diseasevsctrl.csv / ST3-diseasevsctrl.csv).",
            "- Legacy ST1vsST3/ ST3vsST1 summary CSVs are removed to avoid version confusion.",
            "- Plotting scripts for reconstructed/proxy views are in ../supfig4_code/.",
        ]
    )
    (out_dir / "README.txt").write_text("\n".join(readme_lines), encoding="utf-8")


def main():
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent
    out_dir = project_root / "source data" / "supfig4_source_data"
    raw_support_dir = out_dir / "raw_support"
    panel_raw_dir = project_root / "source data" / "supfig4_panel_raw"

    copied_from_panel_raw = copy_panel_raw_to_support(panel_raw_dir, raw_support_dir)
    copied_from_external = copy_external_st_gsea_raw(project_root, raw_support_dir)
    copied_raw_files = sorted(set(copied_from_panel_raw + copied_from_external))
    remove_stale_outputs(out_dir)
    build_panel_a(out_dir, raw_support_dir)
    available_regulon_genes = build_panel_b(out_dir, raw_support_dir)
    build_panel_c_e(out_dir, raw_support_dir)
    build_panel_d_f(out_dir, raw_support_dir)
    workbook_name = write_source_data_workbook(out_dir)
    write_mapping_and_readme(out_dir, available_regulon_genes, copied_raw_files, workbook_name)
    print(f"Wrote source data to: {out_dir}")


if __name__ == "__main__":
    main()
