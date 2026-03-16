#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path

GROUP_ORDER = ["Ctrl", "AZFc_Del", "iNOA_B", "iNOA_S", "KS"]
GROUP_MAP_FROM_PREFIX = {
    "Ctrl": "Ctrl",
    "AZFc_Del": "AZFc_Del",
    "iNOA_B": "iNOA_B",
    "iNOA_S": "iNOA_S",
    "KS": "KS",
    "OA": "OA",
}
STAGE_MAP = {"LC1": "LC_a", "LC2": "LC_b", "LC3": "LC_c"}
STAGE_ORDER = ["LC_a", "LC_b", "LC_c"]
PANEL_A_SCORE_SPECS = [
    {
        "score_id": "Luteinizing_hormone",
        "score_label": "Luteinizing hormone Score",
        "file_name": "B3_LCs_P20073103.diff_PRO.h5ad - Gene set enrichment17_expression.csv",
    },
    {
        "score_id": "Steroid_metabolism",
        "score_label": "Steroid Metabolism Score",
        "file_name": "B3_LCs_P20073103.diff_PRO.h5ad - Gene set enrichment10 - Gene set enrichment14_expression.csv",
    },
    {
        "score_id": "Glucocorticoid_receptor_pathway",
        "score_label": "Glucocorticoid receptor pathway score",
        "file_name": "B3_LCs_P20073103.diff_PRO.h5ad - Gene set enrichment7_expression.csv",
    },
]
PANEL_A_META_COLUMNS = {"cell", "Major cell types", "UMAP1", "UMAP2", "TSNE1", "TSNE2"}


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


def infer_group_from_cell_id(cell_id: str) -> str:
    token = str(cell_id).strip()
    for prefix, group in GROUP_MAP_FROM_PREFIX.items():
        if token.startswith(prefix):
            return group
    match = re.match(r"^(Ctrl|OA|AZFc_Del|iNOA_B|iNOA_S|KS)", token)
    if match:
        return GROUP_MAP_FROM_PREFIX.get(match.group(1), match.group(1))
    return ""


def sort_group(value: str):
    return GROUP_ORDER.index(value) if value in GROUP_ORDER else 999


def sort_stage(value: str):
    return STAGE_ORDER.index(value) if value in STAGE_ORDER else 999


def read_gene_list(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        rows = list(reader)
    if not rows:
        return "", []
    title = rows[0][0].strip() if rows[0] else ""
    genes = []
    seen = set()
    for row in rows[1:]:
        if not row:
            continue
        gene = row[0].strip()
        if gene and gene not in seen:
            genes.append(gene)
            seen.add(gene)
    return title, genes


def panel_a_score_order(score_id: str):
    ids = [item["score_id"] for item in PANEL_A_SCORE_SPECS]
    return ids.index(score_id) if score_id in ids else 999


def pick_panel_a_score_column(columns):
    candidates = [col for col in columns if col not in PANEL_A_META_COLUMNS]
    if not candidates:
        return ""
    if len(candidates) == 1:
        return candidates[0]
    preferred = ["Luteinizing hormone.csv", "steroid metabolism.csv", "Glucocorticoid receptor pathway"]
    for item in preferred:
        if item in candidates:
            return item
    return candidates[0]


def panel_a_summary_rows(cell_rows, value_column: str):
    summary_rows = []
    grouped = defaultdict(list)
    labels = {item["score_id"]: item["score_label"] for item in PANEL_A_SCORE_SPECS}

    for row in cell_rows:
        score = str(row.get("score_id", "")).strip()
        stage = str(row.get("stage", "")).strip()
        value = to_float(row.get(value_column))
        if score and stage and value is not None:
            grouped[(score, stage)].append(value)

    for score in [item["score_id"] for item in PANEL_A_SCORE_SPECS]:
        for stage in STAGE_ORDER:
            values = grouped.get((score, stage), [])
            summary_rows.append(
                {
                    "score_id": score,
                    "score_label": labels.get(score, ""),
                    "stage": stage,
                    "n_cells": len(values),
                    "mean": (sum(values) / len(values)) if values else "",
                    "median": quantile(values, 0.5),
                    "q1": quantile(values, 0.25),
                    "q3": quantile(values, 0.75),
                    "min": min(values) if values else "",
                    "max": max(values) if values else "",
                }
            )
    return summary_rows


def build_panel_a_from_raw_scores(raw_dir: Path, out_dir: Path):
    score_files = [raw_dir / item["file_name"] for item in PANEL_A_SCORE_SPECS]
    if not all(path.exists() for path in score_files):
        return False

    cell_rows = []
    for score_spec in PANEL_A_SCORE_SPECS:
        path = raw_dir / score_spec["file_name"]
        rows = read_csv_rows(path)
        if not rows:
            return False

        score_column = pick_panel_a_score_column(list(rows[0].keys()))
        if score_column == "":
            return False

        for row in rows:
            cluster = str(row.get("Major cell types", "")).strip()
            stage = STAGE_MAP.get(cluster, "")
            if stage == "":
                continue

            cell_id = str(row.get("cell", "")).strip()
            group = infer_group_from_cell_id(cell_id)
            score_value = to_float(row.get(score_column))
            if score_value is None:
                continue

            cell_rows.append(
                {
                    "cell": cell_id,
                    "group": group,
                    "cluster": cluster,
                    "stage": stage,
                    "score_id": score_spec["score_id"],
                    "score_label": score_spec["score_label"],
                    "score_value": score_value,
                    "score_column": score_column,
                    "score_source_file": path.name,
                    "UMAP1": row.get("UMAP1", ""),
                    "UMAP2": row.get("UMAP2", ""),
                    "TSNE1": row.get("TSNE1", ""),
                    "TSNE2": row.get("TSNE2", ""),
                }
            )

    cell_rows.sort(
        key=lambda r: (
            panel_a_score_order(r["score_id"]),
            sort_stage(r["stage"]),
            sort_group(r["group"]),
            r["cell"],
        )
    )

    write_csv(
        out_dir / "SupFig5A_score_cells.csv",
        cell_rows,
        [
            "cell",
            "group",
            "cluster",
            "stage",
            "score_id",
            "score_label",
            "score_value",
            "score_column",
            "score_source_file",
            "UMAP1",
            "UMAP2",
            "TSNE1",
            "TSNE2",
        ],
    )

    summary_rows = panel_a_summary_rows(cell_rows, "score_value")
    write_csv(
        out_dir / "SupFig5A_score_summary.csv",
        summary_rows,
        [
            "score_id",
            "score_label",
            "stage",
            "n_cells",
            "mean",
            "median",
            "q1",
            "q3",
            "min",
            "max",
        ],
    )

    write_csv(
        out_dir / "SupFig5A_score_note.csv",
        [
            {
                "note": "Panel A raw score tables were found and exported directly from three Gene set enrichment files.",
                "raw_score_tables": ";".join([f"raw_support/{item['file_name']}" for item in PANEL_A_SCORE_SPECS]),
                "stage_mapping": "LC1->LC_a;LC2->LC_b;LC3->LC_c",
            }
        ],
        ["note", "raw_score_tables", "stage_mapping"],
    )
    return True


def build_panel_a(raw_dir: Path, out_dir: Path):
    if build_panel_a_from_raw_scores(raw_dir, out_dir):
        return

    expression_rows = read_csv_rows(raw_dir / "SupFig5A_expression.csv")
    if not expression_rows:
        return

    expr_columns = list(expression_rows[0].keys())
    expr_gene_to_col = {}
    suffix = " normalised expression value"
    for col in expr_columns:
        if col.endswith(suffix):
            gene = col[: -len(suffix)]
            expr_gene_to_col[gene] = col

    overlap_rows = []
    cell_rows = []
    legacy_map = {
        "Luteinizing_hormone": raw_dir / "SupFig5A_Luteinizing_hormone_genes.csv",
        "Steroid_metabolism": raw_dir / "SupFig5A_steroid_metabolism_genes.csv",
        "Glucocorticoid_receptor_pathway": raw_dir / "SupFig5A_glucocorticoid_receptor_pathway_genes.csv",
    }

    for score in PANEL_A_SCORE_SPECS:
        gene_file = legacy_map[score["score_id"]]
        title, genes = read_gene_list(gene_file)

        overlap_genes = [gene for gene in genes if gene in expr_gene_to_col]
        missing_genes = [gene for gene in genes if gene not in expr_gene_to_col]

        overlap_rows.append(
            {
                "score_id": score["score_id"],
                "score_label": score["score_label"],
                "gene_set_title": title,
                "gene_set_size": len(genes),
                "overlap_gene_count": len(overlap_genes),
                "overlap_genes": ";".join(overlap_genes),
                "missing_gene_count": len(missing_genes),
                "missing_genes_preview": ";".join(missing_genes[:30]),
                "source_gene_file": str(gene_file.name),
            }
        )

        for row in expression_rows:
            cluster = str(row.get("Major cell types", "")).strip()
            stage = STAGE_MAP.get(cluster, "")
            if stage == "":
                continue

            cell_id = str(row.get("cell", "")).strip()
            group = infer_group_from_cell_id(cell_id)

            values = []
            for gene in overlap_genes:
                col = expr_gene_to_col.get(gene)
                value = to_float(row.get(col)) if col else None
                if value is not None:
                    values.append(value)

            score_value = (sum(values) / len(values)) if values else ""
            cell_rows.append(
                {
                    "cell": cell_id,
                    "group": group,
                    "cluster": cluster,
                    "stage": stage,
                    "score_id": score["score_id"],
                    "score_label": score["score_label"],
                    "score_proxy": score_value,
                    "overlap_gene_count": len(overlap_genes),
                    "overlap_genes": ";".join(overlap_genes),
                }
            )

    cell_rows.sort(
        key=lambda r: (
            panel_a_score_order(r["score_id"]),
            sort_stage(r["stage"]),
            sort_group(r["group"]),
            r["cell"],
        )
    )

    write_csv(
        out_dir / "SupFig5A_score_proxy_cells.csv",
        cell_rows,
        [
            "cell",
            "group",
            "cluster",
            "stage",
            "score_id",
            "score_label",
            "score_proxy",
            "overlap_gene_count",
            "overlap_genes",
        ],
    )

    summary_rows = panel_a_summary_rows(cell_rows, "score_proxy")
    write_csv(
        out_dir / "SupFig5A_score_proxy_summary.csv",
        summary_rows,
        [
            "score_id",
            "score_label",
            "stage",
            "n_cells",
            "mean",
            "median",
            "q1",
            "q3",
            "min",
            "max",
        ],
    )

    write_csv(
        out_dir / "SupFig5A_gene_set_overlap_details.csv",
        overlap_rows,
        [
            "score_id",
            "score_label",
            "gene_set_title",
            "gene_set_size",
            "overlap_gene_count",
            "overlap_genes",
            "missing_gene_count",
            "missing_genes_preview",
            "source_gene_file",
        ],
    )

    write_csv(
        out_dir / "SupFig5A_proxy_note.csv",
        [
            {
                "note": "Cell-level score tables for panel A were not found in this folder. Proxy scores were reconstructed as mean normalized expression of overlap genes between each gene set and SupFig5A_expression.csv marker columns.",
                "raw_violin_panels": "raw_support/SupFig5A_Luteinizing_hormone_violin.pdf;raw_support/SupFig5A_steroid_metabolism_score_violin.pdf;raw_support/SupFig5A_Glucocorticoid_receptor_pathway_score_violin.pdf",
                "raw_expression": "raw_support/SupFig5A_expression.csv",
            }
        ],
        ["note", "raw_violin_panels", "raw_expression"],
    )


def parse_go_rows(path: Path, direction: str):
    rows = read_csv_rows(path)
    parsed = []
    for row in rows:
        group = str(row.get("ONTOLOGY", "")).strip()
        pvalue = to_float(row.get("pvalue"))
        padj = to_float(row.get("p.adjust"))
        count = to_float(row.get("Count"))

        parsed.append(
            {
                "direction": direction,
                "group": group,
                "ID": row.get("ID", ""),
                "Description": row.get("Description", ""),
                "GeneRatio": row.get("GeneRatio", ""),
                "BgRatio": row.get("BgRatio", ""),
                "pvalue": pvalue,
                "p.adjust": padj,
                "qvalue": to_float(row.get("qvalue")),
                "geneID": row.get("geneID", ""),
                "Count": int(count) if count is not None else "",
                "minus_log10_pvalue": (-math.log10(pvalue)) if (pvalue and pvalue > 0) else "",
                "minus_log10_p_adjust": (-math.log10(padj)) if (padj and padj > 0) else "",
            }
        )
    return parsed


def build_panel_b(raw_dir: Path, out_dir: Path):
    down_rows = parse_go_rows(raw_dir / "SupFig5B_GO_enrichment_down.csv", "Down")
    up_rows = parse_go_rows(raw_dir / "SupFig5B_GO_enrichment_up.csv", "Up")

    up_rows.sort(key=lambda r: (sort_group(r["group"]), r["pvalue"] if r["pvalue"] is not None else 999))
    down_rows.sort(key=lambda r: (sort_group(r["group"]), r["pvalue"] if r["pvalue"] is not None else 999))

    write_csv(
        out_dir / "SupFig5B_GO_up_long.csv",
        up_rows,
        [
            "direction",
            "group",
            "Description",
            "ID",
            "Count",
            "GeneRatio",
            "BgRatio",
            "pvalue",
            "p.adjust",
            "qvalue",
            "minus_log10_pvalue",
            "minus_log10_p_adjust",
            "geneID",
        ],
    )

    write_csv(
        out_dir / "SupFig5B_GO_down_long.csv",
        down_rows,
        [
            "direction",
            "group",
            "Description",
            "ID",
            "Count",
            "GeneRatio",
            "BgRatio",
            "pvalue",
            "p.adjust",
            "qvalue",
            "minus_log10_pvalue",
            "minus_log10_p_adjust",
            "geneID",
        ],
    )

    all_rows = up_rows + down_rows
    all_rows.sort(
        key=lambda r: (
            ["Up", "Down"].index(r["direction"]),
            sort_group(r["group"]),
            r["pvalue"] if r["pvalue"] is not None else 999,
        )
    )

    write_csv(
        out_dir / "SupFig5B_GO_updown_long.csv",
        all_rows,
        [
            "direction",
            "group",
            "Description",
            "ID",
            "Count",
            "GeneRatio",
            "BgRatio",
            "pvalue",
            "p.adjust",
            "qvalue",
            "minus_log10_pvalue",
            "minus_log10_p_adjust",
            "geneID",
        ],
    )

    top_rows = []
    for direction, rows in [("Up", up_rows), ("Down", down_rows)]:
        for group in ["AZFc", "iNOA_B", "iNOA_S", "KS"]:
            group_rows = [row for row in rows if row["group"] == group]
            group_rows.sort(key=lambda r: r["pvalue"] if r["pvalue"] is not None else 999)
            for rank, row in enumerate(group_rows[:15], start=1):
                top_rows.append(
                    {
                        "direction": direction,
                        "group": group,
                        "rank_in_group": rank,
                        "Description": row["Description"],
                        "ID": row["ID"],
                        "Count": row["Count"],
                        "pvalue": row["pvalue"],
                        "p.adjust": row["p.adjust"],
                        "minus_log10_pvalue": row["minus_log10_pvalue"],
                        "minus_log10_p_adjust": row["minus_log10_p_adjust"],
                        "geneID": row["geneID"],
                    }
                )

    write_csv(
        out_dir / "SupFig5B_GO_updown_top15_by_group.csv",
        top_rows,
        [
            "direction",
            "group",
            "rank_in_group",
            "Description",
            "ID",
            "Count",
            "pvalue",
            "p.adjust",
            "minus_log10_pvalue",
            "minus_log10_p_adjust",
            "geneID",
        ],
    )

    write_csv(
        out_dir / "SupFig5B_GO_up_top15_by_group.csv",
        [row for row in top_rows if row["direction"] == "Up"],
        [
            "direction",
            "group",
            "rank_in_group",
            "Description",
            "ID",
            "Count",
            "pvalue",
            "p.adjust",
            "minus_log10_pvalue",
            "minus_log10_p_adjust",
            "geneID",
        ],
    )

    write_csv(
        out_dir / "SupFig5B_GO_down_top15_by_group.csv",
        [row for row in top_rows if row["direction"] == "Down"],
        [
            "direction",
            "group",
            "rank_in_group",
            "Description",
            "ID",
            "Count",
            "pvalue",
            "p.adjust",
            "minus_log10_pvalue",
            "minus_log10_p_adjust",
            "geneID",
        ],
    )


def main():
    code_dir = Path(__file__).resolve().parent
    out_dir = code_dir.parent / "supfig5_source_data"
    raw_dir = out_dir / "raw_support"

    build_panel_a(raw_dir, out_dir)
    build_panel_b(raw_dir, out_dir)

    print("Supplementary Figure 5 source-data tables generated in:", out_dir)


if __name__ == "__main__":
    main()
