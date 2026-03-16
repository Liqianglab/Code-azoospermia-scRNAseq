#!/usr/bin/env python3
"""Generate Supplementary Fig1 Source Data tables from available exports."""

from __future__ import annotations

from pathlib import Path
import math
import pandas as pd

HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]

# Locate data folder (current layout or nested "source data")
DATA_CANDIDATES = [ROOT / "figure1", ROOT, ROOT / "source data"]
SRC = None
for cand in DATA_CANDIDATES:
    if (cand / "dot.csv").exists():
        SRC = cand
        break
if SRC is None:
    raise FileNotFoundError("dot.csv not found under expected locations")

OUT = ROOT / "supfig1_source_data"
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

# -------------------- Sup Fig 1A: Feature plots --------------------
marker_list = [
    "BEX1", "WT1", "PECAM1", "CDH5", "INSL3", "IGF1",
    "CD163", "CD14", "FAM129A", "ACTA2", "UTF1", "PIWIL4",
    "MKI67", "AURKA", "STRA8", "DMC1", "BRCA1", "PSMC3IP",
    "PRM2", "SPEM1",
]

feature_out = OUT / "SupFig1A_featureplot_long.csv"
found_markers = set()
long_dfs = []

def _to_long(df: pd.DataFrame, present: list[str]) -> pd.DataFrame:
    value_cols = [f"{g} normalised expression value" for g in present]
    meta_cols = [c for c in ["cell_id", "raw_cluster", "cell_type", "UMAP1", "UMAP2"] if c in df.columns]
    long = df[meta_cols + value_cols].melt(
        id_vars=meta_cols,
        var_name="gene",
        value_name="expression",
    )
    long["gene"] = long["gene"].str.replace(" normalised expression value", "", regex=False)
    return long

DOT = SRC / "dot.csv"
if DOT.exists():
    dot = pd.read_csv(DOT)
    dot = dot.rename(columns={"cell": "cell_id", "Raw cluster": "raw_cluster"})
    if "raw_cluster" in dot.columns:
        dot["cell_type"] = dot["raw_cluster"].map(CLUSTER_TO_CELLTYPE)

    present = [g for g in marker_list if f"{g} normalised expression value" in dot.columns]
    found_markers.update(present)
    if present:
        long_dfs.append(_to_long(dot, present))

# Extra expression export provided for missing markers (if available)
extra_expr = OUT / "expression.csv"
if extra_expr.exists():
    extra = pd.read_csv(extra_expr)
    extra = extra.rename(columns={"cell": "cell_id", "Raw cluster": "raw_cluster"})
    if "Major cell types" in extra.columns:
        extra["cell_type"] = extra["Major cell types"]
    elif "raw_cluster" in extra.columns:
        extra["cell_type"] = extra["raw_cluster"].map(CLUSTER_TO_CELLTYPE)

    present = [g for g in marker_list if f"{g} normalised expression value" in extra.columns]
    found_markers.update(present)
    if present:
        long_dfs.append(_to_long(extra, present))

# Single-gene export (e.g., insl3.csv) if provided
insl3_candidates = [
    OUT / "insl3.csv",
    ROOT / "insl3.csv",
    ROOT / "figure1" / "insl3.csv",
]
insl3_path = next((p for p in insl3_candidates if p.exists()), None)
if insl3_path:
    insl3 = pd.read_csv(insl3_path)
    insl3 = insl3.rename(columns={"cell": "cell_id", "Raw cluster": "raw_cluster"})
    if "Major cell types" in insl3.columns:
        insl3["cell_type"] = insl3["Major cell types"]
    elif "raw_cluster" in insl3.columns:
        insl3["cell_type"] = insl3["raw_cluster"].map(CLUSTER_TO_CELLTYPE)

    present = [g for g in marker_list if f"{g} normalised expression value" in insl3.columns]
    found_markers.update(present)
    if present:
        long_dfs.append(_to_long(insl3, present))

if long_dfs:
    long = pd.concat(long_dfs, ignore_index=True)
    # Drop duplicates if a cell/gene appears in multiple sources
    dedup_cols = [c for c in ["cell_id", "gene"] if c in long.columns]
    if dedup_cols:
        long = long.drop_duplicates(subset=dedup_cols)
    long.to_csv(feature_out, index=False)

missing_markers = [g for g in marker_list if g not in found_markers]

# -------------------- Sup Fig 1B: Heatmap --------------------
heatmap_dir = SRC / "Average expression Heatmap"
avg_heat = heatmap_dir / "average_expression_export.csv"
scaled_heat = heatmap_dir / "scaled_average_expression_export.csv"


def _cluster_to_celltype(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "cluster" in df.columns:
        df["cell_type"] = df["cluster"].map(CLUSTER_TO_CELLTYPE)
    return df


def _summarize_heatmap(df: pd.DataFrame, out_prefix: str) -> None:
    df = _cluster_to_celltype(df)
    if "cell_type" not in df.columns:
        return
    gene_cols = [c for c in df.columns if c not in ["cluster", "cell_type"]]
    grouped = df.groupby("cell_type")[gene_cols].mean(numeric_only=True)

    # Wide: genes as rows, cell types as columns
    wide = grouped.T.reset_index().rename(columns={"index": "gene"})
    wide.to_csv(OUT / f"{out_prefix}_wide.csv", index=False)

    # Long
    long = wide.melt(id_vars=["gene"], var_name="cell_type", value_name="value")
    long.to_csv(OUT / f"{out_prefix}_long.csv", index=False)


if avg_heat.exists():
    _summarize_heatmap(pd.read_csv(avg_heat), "SupFig1B_heatmap_avgexpr")

if scaled_heat.exists():
    _summarize_heatmap(pd.read_csv(scaled_heat), "SupFig1B_heatmap_scaled")

# -------------------- Sup Fig 1B: GO dotplot --------------------
project_root = ROOT.parents[1]
go_preferred_candidates = [
    OUT / "raw_support" / "LC_GOALL_enrichment_sig.xlsx",
    ROOT / "supfig1_source_data" / "raw_support" / "LC_GOALL_enrichment_sig.xlsx",
    project_root / "原来作图的数据" / "无精症" / "LC_GOALL_enrichment_sig.xlsx",
]
GO_PATH = next((p for p in go_preferred_candidates if p.exists()), None)
if GO_PATH is None:
    go_candidates = sorted(project_root.rglob("GOALL_enrichment_sig.csv"))
    preferred = [p for p in go_candidates if "42." in str(p)]
    GO_PATH = preferred[0] if preferred else (go_candidates[0] if go_candidates else None)
if GO_PATH and GO_PATH.exists():
    if GO_PATH.suffix.lower() == ".xlsx":
        go = pd.read_excel(GO_PATH, sheet_name=0)
    else:
        go = pd.read_csv(GO_PATH)
    go = go.rename(columns={"GROUP": "group"})
    # keep common columns
    keep_cols = [c for c in go.columns if c in [
        "group", "ID", "Description", "GeneRatio", "BgRatio", "pvalue",
        "p.adjust", "qvalue", "geneID", "Count"
    ]]
    go = go[keep_cols]

    # Top terms per group by adjusted p-value
    top_n = 5
    go = go.sort_values(["group", "p.adjust", "pvalue"], ascending=[True, True, True])
    top = go.groupby("group", as_index=False).head(top_n).copy()

    # Avoid -inf for p.adjust==0
    min_nonzero = top.loc[top["p.adjust"] > 0, "p.adjust"].min()
    if pd.notna(min_nonzero):
        top["p.adjust"] = top["p.adjust"].replace(0, min_nonzero)
    top["neg_log10_padj"] = -top["p.adjust"].apply(lambda x: math.log10(x) if x > 0 else 0)

    top.to_csv(OUT / "SupFig1B_GO_top_terms.csv", index=False)

# -------------------- README --------------------
readme = OUT / "README.txt"
readme.write_text(
    """Supplementary Fig1 Source Data generated from available exports.

Created files:
- SupFig1A_featureplot_long.csv: per-cell feature-plot data (UMAP1/UMAP2 + expression) for markers present in dot.csv and optional expression.csv.
- SupFig1B_heatmap_avgexpr_wide/long.csv: average expression (by cell type) from average_expression_export.csv.
- SupFig1B_heatmap_scaled_wide/long.csv: scaled expression (by cell type) from scaled_average_expression_export.csv.
- SupFig1B_GO_top_terms.csv: top GO terms per major cell type (from GOALL_enrichment_sig.csv).

Notes:
- Missing markers not present in dot.csv/expression.csv were not included in SupFig1A_featureplot_long.csv.
- Cluster-to-cell-type mapping uses the screenshot iShot_2026-02-12_10.56.46.png.
"""
)

# Write missing marker list for tracking
(pd.DataFrame({"missing_marker": missing_markers})
   .to_csv(OUT / "SupFig1A_missing_markers.csv", index=False))

print(f"Output folder: {OUT}")
