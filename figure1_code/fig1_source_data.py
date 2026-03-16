#!/usr/bin/env python3
"""Generate Fig1 Source Data tables from available exports.

This script does not fabricate missing annotations. It produces:
- Fig1B UMAP-by-group table (from expression.csv)
- Fig1D marker feature table (from expression.csv)
- Fig1E group composition tables (from view_cellular_composition_of_samples.zip)
- Fig1F somatic/germ counts per sample (from view_cellular_composition_of_samples (1)/count.csv)
- Optional sample-level dotplot tables (from compare_by_expression.tar.gz)
- A cluster->celltype template (if only cluster IDs exist)
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd
import zipfile
import tarfile
import re

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT
OUT = ROOT / "fig1_source_data"
OUT.mkdir(parents=True, exist_ok=True)

# Cluster-to-cell-type mapping (from screenshot iShot_2026-02-12_10.56.46.png)
CLUSTER_TO_CELLTYPE = {
    1: "SSCs",
    2: "STs",
    3: "LCs",
    4: "Spermatids",
    5: "SPCs",
    6: "SPCs",
    7: "Myeloid",
    8: "STs",
    9: "SPGs",
    10: "Spermatids",
    11: "Myeloid",
    12: "SSCs",
    13: "Spermatids",
    14: "Myeloid",
    15: "Spermatids",
    16: "LCs",
    17: "SPCs",
    18: "ECs",
    19: "Spermatids",
    20: "Spermatids",
    21: "Spermatids",
    22: "Spermatids",
    23: "SPCs",
    24: "SPCs",
    25: "STs",
    26: "SPCs",
    27: "Myeloid",
    28: "Spermatids",
    29: "LCs",
    30: "Spermatids",
    31: "LCs",
    32: "SPCs",
    33: "Myeloid",
    34: "Spermatids",
    35: "STs",
    36: "Lym",
    37: "Myeloid",
}


def map_group(sample_id: str) -> str:
    s = str(sample_id)
    if s.startswith("Ctrl"):
        return "Ctrl"
    if s.startswith("OA"):
        return "OA"
    if s.startswith("AZFc_Del"):
        return "AZFc_Del"
    if s.startswith("KS"):
        return "KS"
    if s.startswith("iNOA_B"):
        return "iNOA_B"
    if s.startswith("iNOA_S") or s.startswith("iNOS_S"):
        return "iNOA_S"
    return "Unknown"


def normalize_sample_id(sample_id: str) -> str:
    return str(sample_id).strip().replace(" ", "").replace("-", "_")


def load_count_xlsx() -> pd.DataFrame | None:
    candidates = [
        OUT / "raw_support" / "count.xlsx",
        ROOT / "fig1_source_data" / "raw_support" / "count.xlsx",
        ROOT.parent / "原来作图的数据" / "无精症" / "count.xlsx",
    ]
    count_path = next((p for p in candidates if p.exists()), None)
    if count_path is None:
        return None

    count_df = pd.read_excel(count_path, sheet_name=0)
    rename_map = {
        "Sub-group": "sub_group",
        "Age": "age",
        "germline lineages": "germline_lineages",
        "somatic lineages": "somatic_lineages",
        "Total cell umber": "total_cells",
    }
    count_df = count_df.rename(columns=rename_map)

    required = {"sub_group", "age", "germline_lineages", "somatic_lineages", "total_cells"}
    if not required.issubset(count_df.columns):
        return None

    merged = pd.DataFrame({
        "sample_id": count_df["sub_group"].map(normalize_sample_id),
        "age_from_count": pd.to_numeric(count_df["age"], errors="coerce"),
        "somatic_from_count": pd.to_numeric(count_df["somatic_lineages"], errors="coerce"),
        "germ_from_count": pd.to_numeric(count_df["germline_lineages"], errors="coerce"),
        "total_from_count": pd.to_numeric(count_df["total_cells"], errors="coerce"),
    })
    merged = merged.dropna(subset=["sample_id"]).drop_duplicates(subset=["sample_id"], keep="first")
    return merged


# ---------- Fig1B / Fig1D ----------
expr_path = SRC / "expression.csv"
expr = None
gene_cols = []
if expr_path.exists():
    expr = pd.read_csv(expr_path)
    # normalize column names
    expr = expr.rename(columns={"Sample ID": "sample_id", "cell": "cell_id"})
    if "sample_id" in expr.columns:
        expr["group"] = expr["sample_id"].map(map_group)

    # Fig1B UMAP-by-group
    needed = [c for c in ["cell_id", "sample_id", "group", "UMAP1", "UMAP2"] if c in expr.columns]
    if set(["cell_id", "sample_id", "UMAP1", "UMAP2"]).issubset(expr.columns):
        fig1b = expr[needed].copy()
        fig1b.to_csv(OUT / "Fig1B_UMAP_by_group.csv", index=False)

    # Fig1D marker feature table (long)
    gene_cols = [c for c in expr.columns if c.endswith("normalised expression value")]
    if gene_cols:
        meta_cols = [c for c in ["cell_id", "sample_id", "group", "UMAP1", "UMAP2"] if c in expr.columns]
        long = expr[meta_cols + gene_cols].melt(
            id_vars=meta_cols,
            var_name="gene",
            value_name="expression",
        )
        long["gene"] = long["gene"].str.replace(" normalised expression value", "", regex=False)
        long.to_csv(OUT / "Fig1D_marker_feature_long.csv", index=False)

# ---------- Cell type annotation from raw clusters ----------
def build_cell_annotation(raw_df: pd.DataFrame) -> pd.DataFrame:
    raw_df = raw_df.rename(columns={"cell": "cell_id", "Raw cluster": "raw_cluster"})
    raw_df["cell_type"] = raw_df["raw_cluster"].map(CLUSTER_TO_CELLTYPE)
    if expr is not None:
        raw_df = raw_df.merge(
            expr[["cell_id", "sample_id", "group"]].drop_duplicates(),
            on="cell_id",
            how="left",
        )
    return raw_df


dot_path = SRC / "dot.csv"
raw_cluster_path = SRC / "expression (1).csv"

annotation = None
if dot_path.exists():
    dot_df = pd.read_csv(dot_path)
    annotation = build_cell_annotation(dot_df.copy())
    # keep UMAP/TSNE columns if present
    keep_cols = [c for c in annotation.columns if c in ["cell_id", "raw_cluster", "cell_type", "sample_id", "group", "UMAP1", "UMAP2", "TSNE1", "TSNE2"]]
    annotation[keep_cols].to_csv(OUT / "cell_annotation_from_clusters.csv", index=False)

    # Fig1C UMAP by cell type
    if set(["UMAP1", "UMAP2", "cell_type"]).issubset(annotation.columns):
        fig1c_cols = [c for c in ["cell_id", "cell_type", "sample_id", "group", "UMAP1", "UMAP2"] if c in annotation.columns]
        annotation[fig1c_cols].to_csv(OUT / "Fig1C_UMAP_by_celltype.csv", index=False)

    # Fig1G dotplot (full marker list from dot.csv)
    gene_cols_full = [c for c in dot_df.columns if c.endswith("normalised expression value")]
    # drop GAPDH if present (housekeeping, not a marker)
    gene_cols_full = [c for c in gene_cols_full if not c.startswith("GAPDH ")]
    if gene_cols_full:
        dot_df = dot_df.rename(columns={"cell": "cell_id"})
        dot_df = dot_df.merge(annotation[["cell_id", "cell_type"]], on="cell_id", how="left")
        out_rows = []
        for gcol in gene_cols_full:
            gene = gcol.replace(" normalised expression value", "")
            sub = dot_df[["cell_type", gcol]].dropna()
            if sub.empty:
                continue
            stats = sub.groupby("cell_type")[gcol].agg(
                avg_expr="mean",
                pct_expr=lambda x: (x > 0).mean(),
            ).reset_index()
            stats["gene"] = gene
            out_rows.append(stats)
        if out_rows:
            full_dot = pd.concat(out_rows, ignore_index=True)
            full_dot.to_csv(OUT / "Fig1G_dotplot_full_genes.csv", index=False)

elif raw_cluster_path.exists():
    raw = pd.read_csv(raw_cluster_path)
    annotation = build_cell_annotation(raw)
    annotation.to_csv(OUT / "cell_annotation_from_clusters.csv", index=False)

    # Fig1C UMAP by cell type (if UMAP exists)
    if set(["UMAP1", "UMAP2", "cell_type"]).issubset(annotation.columns):
        fig1c_cols = [c for c in ["cell_id", "cell_type", "sample_id", "group", "UMAP1", "UMAP2"] if c in annotation.columns]
        annotation[fig1c_cols].to_csv(OUT / "Fig1C_UMAP_by_celltype.csv", index=False)

    # Fig1G dotplot (available genes only)
    if expr is not None and gene_cols:
        expr_ct = expr.merge(annotation[["cell_id", "cell_type"]], on="cell_id", how="left")
        out_rows = []
        for gcol in gene_cols:
            gene = gcol.replace(" normalised expression value", "")
            sub = expr_ct[["cell_type", gcol]].dropna()
            if sub.empty:
                continue
            stats = sub.groupby("cell_type")[gcol].agg(
                avg_expr="mean",
                pct_expr=lambda x: (x > 0).mean(),
            ).reset_index()
            stats["gene"] = gene
            out_rows.append(stats)
        if out_rows:
            dot = pd.concat(out_rows, ignore_index=True)
            dot.to_csv(OUT / "Fig1G_dotplot_available_genes.csv", index=False)

# ---------- Fig1E group composition ----------
zip_path = SRC / "view_cellular_composition_of_samples.zip"
if zip_path.exists():
    with zipfile.ZipFile(zip_path) as z:
        if "count.csv" in z.namelist():
            count = pd.read_csv(z.open("count.csv"))
            count.to_csv(OUT / "Fig1E_group_composition_count_wide.csv", index=False)
            count_long = count.melt(id_vars=["celltype"], var_name="group", value_name="n_cells")
            count_long.to_csv(OUT / "Fig1E_group_composition_count_long.csv", index=False)
        if "percent.csv" in z.namelist():
            percent = pd.read_csv(z.open("percent.csv"))
            percent.to_csv(OUT / "Fig1E_group_composition_percent_wide.csv", index=False)
            percent_long = percent.melt(id_vars=["celltype"], var_name="group", value_name="percent")
            percent_long.to_csv(OUT / "Fig1E_group_composition_percent_long.csv", index=False)

# ---------- Fig1F somatic/germ counts per sample ----------
sample_dir = SRC / "view_cellular_composition_of_samples (1)"
count_path = sample_dir / "count.csv"
if count_path.exists():
    sample_count = pd.read_csv(count_path)
    # melt to long
    sample_long = sample_count.melt(id_vars=["celltype"], var_name="sample_id", value_name="n_cells")
    sample_long = sample_long[sample_long["sample_id"] != "All samples"]

    somatic_types = {"ECs", "LCs", "Lym", "Myeloid", "Myoid"}
    germ_types = {"SPCs", "SPGs", "SSCs", "STs", "Spermatids"}

    som = (sample_long[sample_long["celltype"].isin(somatic_types)]
           .groupby("sample_id", as_index=False)["n_cells"].sum()
           .rename(columns={"n_cells": "somatic_cell_number"}))
    germ = (sample_long[sample_long["celltype"].isin(germ_types)]
            .groupby("sample_id", as_index=False)["n_cells"].sum()
            .rename(columns={"n_cells": "germ_cell_number"}))
    total = (sample_long.groupby("sample_id", as_index=False)["n_cells"].sum()
             .rename(columns={"n_cells": "total_cells"}))

    fig1f = som.merge(germ, on="sample_id", how="outer").merge(total, on="sample_id", how="outer")
    fig1f["sample_id"] = fig1f["sample_id"].map(normalize_sample_id)
    fig1f["group"] = fig1f["sample_id"].map(map_group)
    fig1f["age"] = pd.NA

    count_meta = load_count_xlsx()
    if count_meta is not None:
        fig1f = fig1f.merge(count_meta, on="sample_id", how="left")
        fig1f["age"] = fig1f["age_from_count"].combine_first(pd.to_numeric(fig1f["age"], errors="coerce"))
        fig1f["somatic_cell_number"] = fig1f["somatic_from_count"].combine_first(pd.to_numeric(fig1f["somatic_cell_number"], errors="coerce"))
        fig1f["germ_cell_number"] = fig1f["germ_from_count"].combine_first(pd.to_numeric(fig1f["germ_cell_number"], errors="coerce"))
        fig1f["total_cells"] = fig1f["total_from_count"].combine_first(pd.to_numeric(fig1f["total_cells"], errors="coerce"))
        fig1f[
            [
                "sample_id",
                "age_from_count",
                "somatic_from_count",
                "germ_from_count",
                "total_from_count",
            ]
        ].to_csv(OUT / "Fig1F_countxlsx_mapping.csv", index=False)
        fig1f = fig1f.drop(columns=["age_from_count", "somatic_from_count", "germ_from_count", "total_from_count"])

    for col in ["somatic_cell_number", "germ_cell_number", "total_cells", "age"]:
        fig1f[col] = pd.to_numeric(fig1f[col], errors="coerce").astype("Int64")
    fig1f = fig1f.sort_values("sample_id").reset_index(drop=True)
    fig1f.to_csv(OUT / "Fig1F_sample_somatic_germ_counts.csv", index=False)

# ---------- Optional: dotplot by sample (from compare_by_expression.tar.gz) ----------
compare_tar = SRC / "compare_by_expression.tar.gz"
if compare_tar.exists():
    with tarfile.open(compare_tar, "r:gz") as tf:
        members = {m.name: m for m in tf.getmembers()}
        for name in ["dotplot_export.csv", "scaled_dotplot_export.csv", "heatmap_export.csv", "scaled_heatmap_export.csv"]:
            if name in members:
                f = tf.extractfile(members[name])
                if f:
                    df = pd.read_csv(f)
                    df.to_csv(OUT / f"compare_{name}", index=False)

# ---------- Cluster->celltype template (manual fill) ----------
cluster_dir = SRC / "rank_genes_results"
if cluster_dir.exists():
    clusters = []
    for p in cluster_dir.glob("cluster *.csv"):
        m = re.search(r"cluster\s+(\d+)", p.name)
        if m:
            clusters.append(int(m.group(1)))
    if clusters:
        tmpl = pd.DataFrame({"cluster_id": sorted(set(clusters)), "cell_type": pd.NA})
        tmpl.to_csv(OUT / "cluster_to_celltype_template.csv", index=False)

# ---------- README ----------
readme = OUT / "README.txt"
readme.write_text(
    """Fig1 Source Data generated from available exports.

Created files:
- Fig1B_UMAP_by_group.csv (from expression.csv)
- Fig1D_marker_feature_long.csv (from expression.csv)
- cell_annotation_from_clusters.csv (from expression (1).csv + cluster->cell type mapping)
- Fig1C_UMAP_by_celltype.csv (from cell annotations, if UMAP available)
- Fig1G_dotplot_full_genes.csv (cell type dotplot from dot.csv, if available)
- Fig1G_dotplot_available_genes.csv (fallback: cell type dotplot for available genes only)
- Fig1E_group_composition_* (from view_cellular_composition_of_samples.zip)
- Fig1F_sample_somatic_germ_counts.csv (from view_cellular_composition_of_samples (1)/count.csv, with age/count override from count.xlsx if present)
- Fig1F_countxlsx_mapping.csv (sample-level mapping pulled from count.xlsx, if present)
- compare_dotplot_export.csv, compare_scaled_dotplot_export.csv (sample-level; not cell-type-level)
- cluster_to_celltype_template.csv (manual fill if needed)

Missing data (not available in exports):
- Donor ages remain NA only if count.xlsx is unavailable.

Notes:
- Cluster-to-cell-type mapping comes from screenshot iShot_2026-02-12_10.56.46.png.
- Fig1G_dotplot_full_genes.csv uses dot.csv (marker list export). GAPDH is excluded as a housekeeping gene.
- Fig1G_dotplot_available_genes.csv only includes genes present in expression.csv.
"""
)
