#!/usr/bin/env python3
"""Build SourceData_Fig1.xlsx from fig1_source_data CSVs.
Splits large sheets to respect Excel row limits.
"""
from __future__ import annotations

from pathlib import Path
import math
import pandas as pd

# Excel allows 1,048,576 rows total including header
EXCEL_MAX_ROWS = 1_048_575

root = Path(__file__).resolve().parents[1]
src = root / "fig1_source_data"
out = root / "SourceData_Fig1.xlsx"

sheet_files = [
    ("Fig1B_UMAP_by_group", src / "Fig1B_UMAP_by_group.csv"),
    ("Fig1C_UMAP_by_celltype", src / "Fig1C_UMAP_by_celltype.csv"),
    ("Fig1D_marker_feature_long", src / "Fig1D_marker_feature_long.csv"),
    ("Fig1E_group_count_long", src / "Fig1E_group_composition_count_long.csv"),
    ("Fig1E_group_percent_long", src / "Fig1E_group_composition_percent_long.csv"),
    ("Fig1F_sample_somatic_germ", src / "Fig1F_sample_somatic_germ_counts.csv"),
    ("Fig1F_countxlsx_mapping", src / "Fig1F_countxlsx_mapping.csv"),
    ("Fig1G_dotplot_full", src / "Fig1G_dotplot_full_genes.csv"),
]

with pd.ExcelWriter(out, engine="openpyxl") as writer:
    for sheet, path in sheet_files:
        if not path.exists():
            continue
        df = pd.read_csv(path)
        if len(df) <= EXCEL_MAX_ROWS:
            df.to_excel(writer, sheet_name=sheet, index=False)
        else:
            # split across multiple sheets
            n_parts = math.ceil(len(df) / EXCEL_MAX_ROWS)
            for i in range(n_parts):
                start = i * EXCEL_MAX_ROWS
                end = min((i + 1) * EXCEL_MAX_ROWS, len(df))
                part = df.iloc[start:end]
                part.to_excel(writer, sheet_name=f"{sheet}_p{i+1}", index=False)

    readme = (
        "SourceData_Fig1.xlsx generated from fig1_source_data CSVs.\n"
        "Sheets:\n"
        "- Fig1B_UMAP_by_group: UMAP coords per cell, grouped by phenotype\n"
        "- Fig1C_UMAP_by_celltype: UMAP coords per cell, grouped by cell type (cluster mapping)\n"
        "- Fig1D_marker_feature_long(_p1/_p2...): per-cell marker expression (long)\n"
        "- Fig1E_group_count_long / percent_long: composition by phenotype\n"
        "- Fig1F_sample_somatic_germ: somatic/germ counts per sample (age from count.xlsx if available)\n"
        "- Fig1F_countxlsx_mapping: age/somatic/germ/total values parsed from count.xlsx\n"
        "- Fig1G_dotplot_full: cell type dotplot (avg_expr, pct_expr)\n"
    )
    pd.DataFrame({"README": [readme]}).to_excel(writer, sheet_name="README", index=False)

print(out)
