#!/usr/bin/env python3
import csv
import math
import shutil
from pathlib import Path


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


def readable_pathway(pathway):
    name = pathway.replace("HALLMARK_", "")
    parts = [part for part in name.split("_") if part]
    return " ".join(part.capitalize() for part in parts)


def read_gsea_table(path, stage, group):
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            pathway = raw.get("Description", "")
            row = {
                "stage": stage,
                "group": group,
                "comparison": f"{group} vs Ctrl",
                "pathway": pathway,
                "pathway_readable": readable_pathway(pathway),
                "setSize": raw.get("setSize", ""),
                "NES": raw.get("NES", ""),
                "pvalue": raw.get("pvalue", ""),
                "p.adjust": raw.get("p.adjust", ""),
                "qvalues": raw.get("qvalues", ""),
                "rank": raw.get("rank", ""),
                "leading_edge": raw.get("leading_edge", ""),
                "core_enrichment": raw.get("core_enrichment", ""),
                "_NES": to_float(raw.get("NES")),
                "_p_adjust": to_float(raw.get("p.adjust")),
            }
            if row["_NES"] is None or row["_p_adjust"] is None:
                continue
            row["neglog10_fdr"] = min(-math.log10(max(row["_p_adjust"], 1e-300)), 20.0)
            rows.append(row)
    return rows


def select_pathways(rows, max_pathways):
    stats = {}
    for row in rows:
        key = row["pathway_readable"]
        val = stats.setdefault(key, {"min_fdr": float("inf"), "max_abs_nes": 0.0})
        val["min_fdr"] = min(val["min_fdr"], row["_p_adjust"])
        val["max_abs_nes"] = max(val["max_abs_nes"], abs(row["_NES"]))
    ranked = sorted(
        stats.items(),
        key=lambda item: (item[1]["min_fdr"], -item[1]["max_abs_nes"], item[0]),
    )
    return [name for name, _ in ranked[:max_pathways]]


def pathway_order(rows):
    stats = {}
    for row in rows:
        key = row["pathway_readable"]
        val = stats.setdefault(key, {"min_fdr": float("inf"), "sum_abs_nes": 0.0, "n": 0})
        val["min_fdr"] = min(val["min_fdr"], row["_p_adjust"])
        val["sum_abs_nes"] += abs(row["_NES"])
        val["n"] += 1
    ranked = sorted(
        stats.items(),
        key=lambda item: (item[1]["min_fdr"], -(item[1]["sum_abs_nes"] / item[1]["n"]), item[0]),
    )
    return [name for name, _ in ranked]


def make_plot_matrix(rows, comparisons, max_pathways):
    selected = set(select_pathways(rows, max_pathways))
    filtered = [row for row in rows if row["pathway_readable"] in selected]
    ordered_pathways = pathway_order(filtered)

    lookup = {}
    for row in filtered:
        lookup[(row["comparison"], row["pathway_readable"])] = row

    matrix = []
    for pathway in ordered_pathways:
        for comparison in comparisons:
            match = lookup.get((comparison, pathway))
            matrix.append(
                {
                    "comparison": comparison,
                    "pathway_readable": pathway,
                    "pathway": "" if match is None else match["pathway"],
                    "NES": 0.0 if match is None else match["_NES"],
                    "p.adjust": "" if match is None else match["p.adjust"],
                    "pvalue": "" if match is None else match["pvalue"],
                    "neglog10_fdr": 0.0 if match is None else match["neglog10_fdr"],
                    "is_missing": 1 if match is None else 0,
                    "stage": "" if match is None else match["stage"],
                    "group": "" if match is None else match["group"],
                }
            )
    return filtered, matrix


def write_csv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            out = {column: row.get(column, "") for column in columns}
            writer.writerow(out)


def copy_raw_support(raw_map, src_root, dst_root):
    dst_root.mkdir(parents=True, exist_ok=True)
    for label, rel_path in raw_map.items():
        source = src_root / rel_path
        target = dst_root / f"{label}_P20073103_GSEA_enrichment.xls"
        shutil.copy2(source, target)


def main():
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent
    gsea_root = project_root / "原来作图的数据" / "汇总" / "11.生殖细胞组间差异分析" / "GSEA"
    out_dir = project_root / "source data" / "supfig3_source_data"
    raw_support_dir = out_dir / "raw_support"

    late_files = {
        "Late_primary_SPCs_CtrlvsOA": Path("gsea_Late_primary_SPCs_CtrlvsOA/P20073103/gsea_Late_primary_SPCs_CtrlvsOA_P20073103_GSEA_enrichment.xls"),
        "Late_primary_SPCs_CtrlvsAZFc_Del": Path("gsea_Late_primary_SPCs_CtrlvsAZFc_Del/P20073103/P20073103_GSEA_enrichment.xls"),
        "Late_primary_SPCs_CtrlvsKS": Path("gsea_Late_primary_SPCs_CtrlvsKS/P20073103/P20073103_GSEA_enrichment.xls"),
        "Late_primary_SPCs_CtrlvsiNOA_B": Path("gsea_Late_primary_SPCs_CtrlvsiNOA_B/P20073103/P20073103_GSEA_enrichment.xls"),
    }
    round_files = {
        "Round_Spermatids_CtrlvsOA": Path("gsea_Round_Spermatids_CtrlvsOA/P20073103/P20073103_GSEA_enrichment.xls"),
        "Round_Spermatids_CtrlvsAZFc_Del": Path("gsea_Round_Spermatids_CtrlvsAZFc_Del/P20073103/P20073103_GSEA_enrichment.xls"),
        "Round_Spermatids_CtrlvsKS": Path("gsea_Round_Spermatids_CtrlvsKS/P20073103/P20073103_GSEA_enrichment.xls"),
    }

    late_rows = []
    for label, rel_path in late_files.items():
        group = label.replace("Late_primary_SPCs_Ctrlvs", "")
        late_rows.extend(read_gsea_table(gsea_root / rel_path, "Late_primary_SPCs", group))

    round_rows = []
    for label, rel_path in round_files.items():
        group = label.replace("Round_Spermatids_Ctrlvs", "")
        round_rows.extend(read_gsea_table(gsea_root / rel_path, "Round_Spermatids", group))

    late_comparisons = ["OA vs Ctrl", "AZFc_Del vs Ctrl", "KS vs Ctrl", "iNOA_B vs Ctrl"]
    round_comparisons = ["OA vs Ctrl", "AZFc_Del vs Ctrl", "KS vs Ctrl"]

    _, late_matrix = make_plot_matrix(late_rows, late_comparisons, max_pathways=24)
    _, round_matrix = make_plot_matrix(round_rows, round_comparisons, max_pathways=28)

    all_cols = [
        "stage",
        "group",
        "comparison",
        "pathway",
        "pathway_readable",
        "setSize",
        "NES",
        "pvalue",
        "p.adjust",
        "qvalues",
        "rank",
        "leading_edge",
        "core_enrichment",
        "neglog10_fdr",
    ]
    matrix_cols = [
        "stage",
        "group",
        "comparison",
        "pathway_readable",
        "pathway",
        "NES",
        "pvalue",
        "p.adjust",
        "neglog10_fdr",
        "is_missing",
    ]

    late_rows_sorted = sorted(
        late_rows,
        key=lambda row: (
            late_comparisons.index(row["comparison"]),
            row["_p_adjust"],
            row["pathway_readable"],
        ),
    )
    round_rows_sorted = sorted(
        round_rows,
        key=lambda row: (
            round_comparisons.index(row["comparison"]),
            row["_p_adjust"],
            row["pathway_readable"],
        ),
    )

    write_csv(out_dir / "SupFig3A_Late_primary_SPCs_GSEA_all.csv", late_rows_sorted, all_cols)
    write_csv(out_dir / "SupFig3B_Round_Spermatids_GSEA_all.csv", round_rows_sorted, all_cols)
    write_csv(out_dir / "SupFig3A_plot_matrix.csv", late_matrix, matrix_cols)
    write_csv(out_dir / "SupFig3B_plot_matrix.csv", round_matrix, matrix_cols)

    raw_map = {}
    raw_map.update(late_files)
    raw_map.update(round_files)
    copy_raw_support(raw_map, gsea_root, raw_support_dir)

    readme = out_dir / "README.txt"
    readme.write_text(
        "\n".join(
            [
                "Supplementary Figure 3 source data (Panels A-B)",
                "",
                "Files:",
                "- SupFig3A_Late_primary_SPCs_GSEA_all.csv: all Hallmark rows from 4 comparisons (OA/AZFc_Del/KS/iNOA_B vs Ctrl).",
                "- SupFig3B_Round_Spermatids_GSEA_all.csv: all Hallmark rows from 3 comparisons (OA/AZFc_Del/KS vs Ctrl).",
                "- SupFig3A_plot_matrix.csv: panel A bubble-matrix values after pathway selection (max_pathways=24).",
                "- SupFig3B_plot_matrix.csv: panel B bubble-matrix values after pathway selection (max_pathways=28).",
                "- raw_support/*.xls: original GSEA enrichment tables used to build the above CSV files.",
                "",
                "Selection logic follows:",
                "- pathway ranking by minimum adjusted p-value, then by maximum absolute NES.",
                "- point size = -log10(FDR), capped at 20.",
            ]
        ),
        encoding="utf-8",
    )

    print(f"Wrote source data to: {out_dir}")


if __name__ == "__main__":
    main()
