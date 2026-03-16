#!/usr/bin/env python3
from __future__ import annotations

import math
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parents[1]  # .../source data
ROOT = BASE.parent
OUT_DIR = BASE / "fig7_source_data"
RAW_SUPPORT_DIR = OUT_DIR / "raw_support"

FIG7_AI = ROOT / "Figure 7_20251210.ai"
FIG7_PDF = ROOT / "Figure 7_20251210.pdf"

GROUP_ORDER = ["Ctrl", "OA", "AZFc_Del", "iNOA_B", "iNOA_S", "KS"]
GROUP_ORDER_PANEL_D = ["Ctrl", "OA", "AZFc_Del", "iNOA_B", "KS"]

# Panel A
SCMETA_DIR = ROOT / "原来作图的数据" / "汇总" / "61.生殖亚群的scMetabolism分析"
PANEL_A_FEATURE_FILES = {
    "Glycolysis": SCMETA_DIR / "featureplot" / "P20073103_Glycolysis_Gluconeogenesis_scMeta_fea.pdf",
    "Citrate_cycle_TCA_cycle": SCMETA_DIR / "featureplot" / "P20073103_Citrate_cycle_TCA_cycle_scMeta_fea.pdf",
    "Oxidative_phosphorylation": SCMETA_DIR / "featureplot" / "P20073103_Oxidative_phosphorylation_scMeta_fea.pdf",
}
PANEL_A_GROUP_DIFF = SCMETA_DIR / "Diff" / "P20073103_group.xls"
PANEL_A_GROUP_DIFF_NO_OA = SCMETA_DIR / "Diff" / "P20073103_group去OA.xls"

# Panel B / D
RAW7_DIR = ROOT / "原来作图的数据" / "7"
PANEL_B_ALL_EXPR = RAW7_DIR / "late阶段所有基因表达.csv"
PANEL_B_DIAGRAM = RAW7_DIR / "late阶段_代谢通路图基因表达_提取.csv"
PANEL_B_LOG2FC_TABLE = RAW7_DIR / "plots" / "late_stage_metabolism_log2fc_table.csv"
PANEL_B_ANNO = RAW7_DIR / "Annotation_Row.csv"
PANEL_B_MAIN_SCRIPT = RAW7_DIR / "make_metabolism_plots.R"
PANEL_B_EXTRACT_SCRIPT = RAW7_DIR / "extract_diagram_genes.R"
PANEL_B_HEATMAP_SCRIPT = RAW7_DIR / "make_diagram_log2fc_heatmap.R"
PANEL_D_INPUT = RAW7_DIR / "B3_shengzhi_P20073103.diff_PRO.h5ad - Extraction-late - ODX_expression.csv"
PANEL_D_PLOT = RAW7_DIR / "plots" / "late_stage_oxphos_score_violin.pdf"

# Panel C
PANEL_C_INPUT = ROOT / "原来作图的数据" / "ST细胞提取 - Gene set enrichment1_expression.csv"
PANEL_C_SCORE_COL = "Lactate_Glycolysis_Signature.csv"
PANEL_C_VIOLIN_PDF = ROOT / "原来作图的数据" / "metabolism_lactate_shuttle" / "ST_lactate_export_violin_MCT4_LDHA.pdf"

# Panel E
PANEL_E_DIR = Path("/Users/xq/Desktop/睾丸相关/睾丸补图/figure7")
PANEL_E_SCRIPT = PANEL_E_DIR / "make_dotplots.R"
PANEL_E_INPUTS = [
    ("ST", "ST -ldha-slc16a3.csv", ["SLC16A3", "LDHA"]),
    ("Late_primary_SPCs", "late-spc-slc16a1-ldhb.csv", ["SLC16A1", "LDHB"]),
    ("Round_spermatids", "round-spc-slc16a1-ldhb.csv", ["SLC16A1", "LDHB"]),
]
PANEL_E_PLOTS = [
    PANEL_E_DIR / "plots" / "ST_lactate_export_dotplot.pdf",
    PANEL_E_DIR / "plots" / "Late_primary_SPCs_lactate_uptake_dotplot.pdf",
    PANEL_E_DIR / "plots" / "Round_spermatids_lactate_uptake_dotplot.pdf",
    PANEL_E_DIR / "plots" / "Figure7_lactate_shuttle_dotplots_combined.pdf",
]

# Panel F
PANEL_F_DIRS = {
    "CTRL": Path("/Users/xq/Desktop/未命名文件夹/control-mct1-red-scp3-green-1-20x"),
    "INOA_B": Path("/Users/xq/Desktop/未命名文件夹/inob-mct1-red-scp3-green-2-20x"),
    "KS": Path("/Users/xq/Desktop/未命名文件夹/KS-mct1-red-scp3-green-2-20x"),
    "AZFc_Del": Path("/Users/xq/Desktop/未命名文件夹/AZ-mct1-red-scp3-green-1-20x"),
}


def ensure_dirs() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RAW_SUPPORT_DIR.mkdir(parents=True, exist_ok=True)


def copy_to_raw_support(src: Path, alias: str | None = None) -> str:
    if not src.exists() or not src.is_file():
        return ""
    name = alias if alias else src.name
    dst = RAW_SUPPORT_DIR / name
    shutil.copy2(src, dst)
    return f"raw_support/{name}"


def infer_group_from_cell(cell: str) -> str:
    token = str(cell)
    m = re.match(r"^(Ctrl|OA|AZFc_Del|iNOA_B|iNOA_S|KS)", token)
    return m.group(1) if m else ""


def safe_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def summarize_score(df: pd.DataFrame, group_col: str, score_col: str, group_order: list[str]) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(
            columns=[group_col, "n_cells", "mean", "median", "q1", "q3", "min", "max"]
        )

    g = df.groupby(group_col, observed=False)[score_col]
    summary = (
        g.agg(n_cells="count", mean="mean", median="median", min="min", max="max")
        .reset_index()
        .merge(g.quantile(0.25).reset_index(name="q1"), on=group_col, how="left")
        .merge(g.quantile(0.75).reset_index(name="q3"), on=group_col, how="left")
    )
    summary[group_col] = pd.Categorical(summary[group_col], categories=group_order, ordered=True)
    summary = summary.sort_values(group_col).reset_index(drop=True)
    return summary


def classify_oxphos_complex(gene: str) -> str:
    g = str(gene)
    if g.startswith("COQ"):
        return "CoQ"
    if g.startswith("MT-ND") or g.startswith("NDUF"):
        return "Complex I"
    if g.startswith("SDHAF") or g in {"SDHA", "SDHB", "SDHC", "SDHD"}:
        return "Complex II"
    if g.startswith("UQCR") or g == "MT-CYB":
        return "Complex III"
    if g.startswith("COX") or g.startswith("MT-CO"):
        return "Complex IV"
    if g.startswith("ATP5") or g.startswith("MT-ATP"):
        return "Complex V"
    return "Other"


def build_panel_a() -> dict:
    generated = []
    raw_inputs = []

    # Featureplot registry
    feature_rows = []
    for pathway, fpath in PANEL_A_FEATURE_FILES.items():
        rel = copy_to_raw_support(fpath, f"Fig7A_{fpath.name}")
        if rel:
            raw_inputs.append(rel)
        feature_rows.append(
            {
                "pathway": pathway,
                "featureplot_pdf": str(fpath),
                "exists": fpath.exists(),
                "raw_support_copy": rel,
            }
        )

    pd.DataFrame(feature_rows).to_csv(OUT_DIR / "Fig7A_featureplot_file_registry.csv", index=False)
    generated.append("Fig7A_featureplot_file_registry.csv")

    # Group-level scMeta matrix (proxy numeric support)
    pathway_map = {
        "Glycolysis / Gluconeogenesis": "Glycolysis",
        "Citrate cycle (TCA cycle)": "Citrate_cycle_TCA_cycle",
        "Oxidative phosphorylation": "Oxidative_phosphorylation",
    }

    long_frames = []
    for source_name, table_path, alias in [
        ("with_OA", PANEL_A_GROUP_DIFF, "Fig7A_P20073103_group.xls"),
        ("no_OA", PANEL_A_GROUP_DIFF_NO_OA, "Fig7A_P20073103_group_no_OA.xls"),
    ]:
        if not table_path.exists():
            continue
        rel = copy_to_raw_support(table_path, alias)
        if rel:
            raw_inputs.append(rel)

        df = pd.read_csv(table_path, sep="\t")
        id_col = "id" if "id" in df.columns else df.columns[0]
        df = df[df[id_col].isin(pathway_map.keys())].copy()
        if df.empty:
            continue

        melted = df.melt(id_vars=[id_col], var_name="sample", value_name="score")
        split = melted["sample"].str.extract(
            r"^(?P<cell_type>.+)_(?P<group>Ctrl|OA|AZFc_Del|iNOA_B|iNOA_S|KS)$"
        )
        melted["cell_type"] = split["cell_type"]
        melted["group"] = split["group"]
        melted["score"] = safe_num(melted["score"])
        melted["pathway"] = melted[id_col].map(pathway_map)
        melted["source_table"] = source_name
        melted = melted[["source_table", "pathway", "cell_type", "group", "score"]].dropna(
            subset=["group"]
        )
        long_frames.append(melted)

    if long_frames:
        long_df = pd.concat(long_frames, ignore_index=True)
        group_order_a = ["Ctrl", "OA", "AZFc_Del", "iNOA_B", "iNOA_S", "KS"]
        long_df["group"] = pd.Categorical(long_df["group"], categories=group_order_a, ordered=True)
        long_df = long_df.sort_values(["pathway", "cell_type", "group", "source_table"]).reset_index(
            drop=True
        )
        long_df.to_csv(OUT_DIR / "Fig7A_scmeta_diff_selected_pathways_long.csv", index=False)
        generated.append("Fig7A_scmeta_diff_selected_pathways_long.csv")

        late_df = long_df[long_df["cell_type"] == "Late_primary_SPCs"].copy()
        late_df.to_csv(OUT_DIR / "Fig7A_scmeta_diff_late_primary_pathways.csv", index=False)
        generated.append("Fig7A_scmeta_diff_late_primary_pathways.csv")

    note = pd.DataFrame(
        [
            {
                "panel": "A",
                "note": "Panel A uses exported scMetabolism featureplot PDFs. In this workspace, no direct per-cell score table for these three UMAPs was found; group-level scMeta diff matrices are provided as numeric support.",
            }
        ]
    )
    note.to_csv(OUT_DIR / "Fig7A_source_note.csv", index=False)
    generated.append("Fig7A_source_note.csv")

    return {
        "subpanel": "Metabolism score UMAPs (Glycolysis/TCA/OXPHOS)",
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../fig7_code/fig7_source_data.py"],
        "status": "partial",
        "availability_note": "Featureplot exports available; per-cell score matrix not located in this folder",
    }


def build_panel_b() -> dict:
    generated = []
    raw_inputs = []
    code_refs = ["../fig7_code/fig7_source_data.py"]

    for src, alias in [
        (PANEL_B_ALL_EXPR, "Fig7B_late_stage_all_genes_expression.csv"),
        (PANEL_B_DIAGRAM, "Fig7B_late_stage_diagram_genes_extracted.csv"),
        (PANEL_B_LOG2FC_TABLE, "Fig7B_late_stage_metabolism_log2fc_table.csv"),
        (PANEL_B_ANNO, "Fig7B_Annotation_Row.csv"),
        (PANEL_B_MAIN_SCRIPT, "Fig7B_make_metabolism_plots.R"),
        (PANEL_B_EXTRACT_SCRIPT, "Fig7B_extract_diagram_genes.R"),
        (PANEL_B_HEATMAP_SCRIPT, "Fig7B_make_diagram_log2fc_heatmap.R"),
        (RAW7_DIR / "plots" / "late_stage_diagram_genes_log2fc_heatmap.pdf", "Fig7B_late_stage_diagram_genes_log2fc_heatmap.pdf"),
    ]:
        rel = copy_to_raw_support(src, alias)
        if rel:
            raw_inputs.append(rel)

    for script_ref in [
        "raw_support/Fig7B_make_metabolism_plots.R",
        "raw_support/Fig7B_extract_diagram_genes.R",
        "raw_support/Fig7B_make_diagram_log2fc_heatmap.R",
    ]:
        if (RAW_SUPPORT_DIR / script_ref.split("/", 1)[1]).exists():
            code_refs.append(script_ref)

    eps = 1e-6

    if PANEL_B_ALL_EXPR.exists():
        expr = pd.read_csv(PANEL_B_ALL_EXPR).rename(
            columns={"Gene": "gene", "CTRL": "ctrl", "Disease": "disease"}
        )
        for col in ["ctrl", "disease"]:
            expr[col] = safe_num(expr[col])
        expr["log2fc"] = ((expr["disease"] + eps) / (expr["ctrl"] + eps)).apply(
            lambda x: math.log2(x) if pd.notna(x) and x > 0 else float("nan")
        )
        expr = expr[["gene", "ctrl", "disease", "log2fc"]]
        expr.to_csv(OUT_DIR / "Fig7B_full_expression_log2fc.csv", index=False)
        generated.append("Fig7B_full_expression_log2fc.csv")

    if PANEL_B_DIAGRAM.exists():
        diag = pd.read_csv(PANEL_B_DIAGRAM).copy()
        # Normalize expected columns
        rename_map = {
            "Gene": "gene",
            "CTRL": "ctrl",
            "Disease": "disease",
        }
        diag = diag.rename(columns=rename_map)
        for col in ["ctrl", "disease"]:
            if col in diag.columns:
                diag[col] = safe_num(diag[col])

        if "log2fc" not in diag.columns and {"ctrl", "disease"}.issubset(diag.columns):
            diag["log2fc"] = ((diag["disease"] + eps) / (diag["ctrl"] + eps)).apply(
                lambda x: math.log2(x) if pd.notna(x) and x > 0 else float("nan")
            )

        keep_cols = [
            c
            for c in ["pathway", "subcategory", "gene", "diagram_order", "ctrl", "disease", "log2fc"]
            if c in diag.columns
        ]
        if "diagram_order" in diag.columns:
            diag = diag.sort_values("diagram_order")
        diag = diag[keep_cols].reset_index(drop=True)
        diag.to_csv(OUT_DIR / "Fig7B_diagram_gene_log2fc.csv", index=False)
        generated.append("Fig7B_diagram_gene_log2fc.csv")

        if "pathway" in diag.columns and "gene" in diag.columns:
            ox = diag[diag["pathway"].astype(str).str.contains("Oxidative", case=False, na=False)].copy()
            ox["complex"] = ox["gene"].apply(classify_oxphos_complex)
            if "log2fc" in ox.columns:
                summary = (
                    ox.groupby("complex", as_index=False)
                    .agg(
                        n_genes=("gene", "count"),
                        mean_log2fc=("log2fc", "mean"),
                        median_log2fc=("log2fc", "median"),
                    )
                    .sort_values(
                        "complex",
                        key=lambda s: s.map(
                            {
                                "Complex I": 1,
                                "Complex II": 2,
                                "Complex III": 3,
                                "Complex IV": 4,
                                "Complex V": 5,
                                "CoQ": 6,
                                "Other": 7,
                            }
                        ),
                    )
                    .reset_index(drop=True)
                )
                summary.to_csv(OUT_DIR / "Fig7B_oxphos_complex_summary.csv", index=False)
                generated.append("Fig7B_oxphos_complex_summary.csv")

    note = pd.DataFrame(
        [
            {
                "panel": "B",
                "note": "The final metabolic pathway cartoon is AI-laid out; gene-level log2FC values and scripts used to derive/replot heatmaps are provided.",
            }
        ]
    )
    note.to_csv(OUT_DIR / "Fig7B_source_note.csv", index=False)
    generated.append("Fig7B_source_note.csv")

    return {
        "subpanel": "Metabolic pathway diagram (log2FC in Late_primary_SPCs)",
        "generated": generated,
        "raw": raw_inputs,
        "code": code_refs,
        "status": "available",
        "availability_note": "Diagram gene log2FC table and source scripts are available",
    }


def build_panel_c() -> dict:
    generated = []
    raw_inputs = []

    rel = copy_to_raw_support(PANEL_C_INPUT, "Fig7C_ST_GSEA_expression.csv")
    if rel:
        raw_inputs.append(rel)
    rel = copy_to_raw_support(PANEL_C_VIOLIN_PDF, "Fig7C_ST_lactate_export_violin_MCT4_LDHA.pdf")
    if rel:
        raw_inputs.append(rel)

    if not PANEL_C_INPUT.exists():
        return {
            "subpanel": "STs_Lactate/Glycolysis violin",
            "generated": generated,
            "raw": raw_inputs,
            "code": ["../fig7_code/fig7_source_data.py"],
            "status": "missing",
            "availability_note": "Input CSV not found",
        }

    df = pd.read_csv(PANEL_C_INPUT)
    df["group"] = df["cell"].astype(str).apply(infer_group_from_cell)
    df["score"] = safe_num(df[PANEL_C_SCORE_COL])

    keep_cols = [
        "cell",
        "group",
        "Major cell types",
        "score",
        "UMAP1",
        "UMAP2",
        "TSNE1",
        "TSNE2",
    ]
    cells = df[keep_cols].copy()
    cells = cells[cells["group"] != ""].reset_index(drop=True)
    cells["group"] = pd.Categorical(cells["group"], categories=GROUP_ORDER, ordered=True)
    cells = cells.sort_values(["group", "Major cell types", "cell"]).reset_index(drop=True)
    cells.to_csv(OUT_DIR / "Fig7C_ST_lactate_glycolysis_cells.csv", index=False)
    generated.append("Fig7C_ST_lactate_glycolysis_cells.csv")

    summary = summarize_score(cells, "group", "score", GROUP_ORDER)
    summary.to_csv(OUT_DIR / "Fig7C_ST_lactate_glycolysis_summary.csv", index=False)
    generated.append("Fig7C_ST_lactate_glycolysis_summary.csv")

    subtype_summary = summarize_score(
        cells.rename(columns={"Major cell types": "subtype"}), "subtype", "score", ["ST_a", "ST_b", "ST_c"]
    )
    subtype_summary.to_csv(OUT_DIR / "Fig7C_ST_lactate_glycolysis_by_subtype_summary.csv", index=False)
    generated.append("Fig7C_ST_lactate_glycolysis_by_subtype_summary.csv")

    return {
        "subpanel": "STs_Lactate/Glycolysis violin",
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../fig7_code/fig7_source_data.py"],
        "status": "available",
        "availability_note": "Cell-level score table and group summaries are available",
    }


def build_panel_d() -> dict:
    generated = []
    raw_inputs = []

    rel = copy_to_raw_support(PANEL_D_INPUT, "Fig7D_late_stage_ODX_expression.csv")
    if rel:
        raw_inputs.append(rel)
    rel = copy_to_raw_support(PANEL_D_PLOT, "Fig7D_late_stage_oxphos_score_violin.pdf")
    if rel:
        raw_inputs.append(rel)

    if not PANEL_D_INPUT.exists():
        return {
            "subpanel": "Late_primary_SPCs_Oxidative phosphorylation violin",
            "generated": generated,
            "raw": raw_inputs,
            "code": ["../fig7_code/fig7_source_data.py"],
            "status": "missing",
            "availability_note": "Input CSV not found",
        }

    df = pd.read_csv(PANEL_D_INPUT)
    score_col = "oxidative phosphorylation.csv"
    df = df.rename(columns={"gname": "group", score_col: "score"}).copy()
    df["score"] = safe_num(df["score"])
    df = df[df["group"].isin(GROUP_ORDER_PANEL_D)].copy()

    keep_cols = ["cell", "group", "score", "UMAP1", "UMAP2", "TSNE1", "TSNE2"]
    cells = df[keep_cols].copy()
    cells["group"] = pd.Categorical(cells["group"], categories=GROUP_ORDER_PANEL_D, ordered=True)
    cells = cells.sort_values(["group", "cell"]).reset_index(drop=True)
    cells.to_csv(OUT_DIR / "Fig7D_late_primary_SPCs_oxphos_cells.csv", index=False)
    generated.append("Fig7D_late_primary_SPCs_oxphos_cells.csv")

    summary = summarize_score(cells, "group", "score", GROUP_ORDER_PANEL_D)
    summary.to_csv(OUT_DIR / "Fig7D_late_primary_SPCs_oxphos_summary.csv", index=False)
    generated.append("Fig7D_late_primary_SPCs_oxphos_summary.csv")

    return {
        "subpanel": "Late_primary_SPCs_Oxidative phosphorylation violin",
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../fig7_code/fig7_source_data.py", "raw_support/Fig7B_make_metabolism_plots.R"],
        "status": "available",
        "availability_note": "Original ODX score table and group summaries are available",
    }


def build_dotplot_input(csv_path: Path, cell_type: str, gene_order: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(csv_path)
    expr_cols = [c for c in df.columns if c.endswith(" raw expression value")]
    if not expr_cols:
        raise RuntimeError(f"No expression columns found in {csv_path}")

    long_df = df[["cell", "gname", *expr_cols]].melt(
        id_vars=["cell", "gname"], value_vars=expr_cols, var_name="gene", value_name="raw_expr"
    )
    long_df["gene"] = long_df["gene"].str.replace(" raw expression value", "", regex=False)
    long_df["group"] = long_df["gname"].astype(str)
    long_df["raw_expr"] = safe_num(long_df["raw_expr"]).fillna(0.0)
    long_df["log1p_expr"] = long_df["raw_expr"].apply(lambda v: math.log1p(v))
    long_df["cell_type"] = cell_type
    long_df = long_df[long_df["group"].isin(GROUP_ORDER)].copy()

    summary = (
        long_df.groupby(["cell_type", "group", "gene"], as_index=False)
        .agg(
            pct_expressed=("raw_expr", lambda s: (s > 0).mean() * 100),
            mean_log1p=("log1p_expr", "mean"),
            mean_raw_expr=("raw_expr", "mean"),
            n_cells=("cell", "count"),
        )
    )

    summary["group"] = pd.Categorical(summary["group"], categories=GROUP_ORDER, ordered=True)
    summary["gene"] = pd.Categorical(summary["gene"], categories=gene_order, ordered=True)
    summary = summary.sort_values(["cell_type", "gene", "group"]).reset_index(drop=True)

    long_df["group"] = pd.Categorical(long_df["group"], categories=GROUP_ORDER, ordered=True)
    long_df["gene"] = pd.Categorical(long_df["gene"], categories=gene_order, ordered=True)
    long_df = long_df.sort_values(["cell_type", "gene", "group", "cell"]).reset_index(drop=True)

    return long_df, summary


def build_panel_e() -> dict:
    generated = []
    raw_inputs = []
    code_refs = ["../fig7_code/fig7_source_data.py"]

    rel = copy_to_raw_support(PANEL_E_SCRIPT, "Fig7E_make_dotplots.R")
    if rel:
        raw_inputs.append(rel)
        code_refs.append(rel)

    for plot_path in PANEL_E_PLOTS:
        rel_plot = copy_to_raw_support(plot_path, f"Fig7E_{plot_path.name}")
        if rel_plot:
            raw_inputs.append(rel_plot)

    all_long = []
    all_sum = []

    for cell_type, fname, gene_order in PANEL_E_INPUTS:
        src = PANEL_E_DIR / fname
        rel = copy_to_raw_support(src, f"Fig7E_{fname}")
        if rel:
            raw_inputs.append(rel)
        if not src.exists():
            continue

        long_df, sum_df = build_dotplot_input(src, cell_type, gene_order)
        long_name = f"Fig7E_{cell_type}_dotplot_cells_long.csv"
        sum_name = f"Fig7E_{cell_type}_dotplot_summary.csv"
        long_df.to_csv(OUT_DIR / long_name, index=False)
        sum_df.to_csv(OUT_DIR / sum_name, index=False)
        generated.extend([long_name, sum_name])

        all_long.append(long_df)
        all_sum.append(sum_df)

    if all_long:
        merged_long = pd.concat(all_long, ignore_index=True)
        merged_long.to_csv(OUT_DIR / "Fig7E_dotplot_cells_long_all.csv", index=False)
        generated.append("Fig7E_dotplot_cells_long_all.csv")

    if all_sum:
        merged_sum = pd.concat(all_sum, ignore_index=True)
        merged_sum.to_csv(OUT_DIR / "Fig7E_dotplot_summary_all.csv", index=False)
        generated.append("Fig7E_dotplot_summary_all.csv")

    note = pd.DataFrame(
        [
            {
                "panel": "E",
                "metric_definition": "pct_expressed = mean(raw_expr > 0) * 100; mean_log1p = mean(log1p(raw_expr)); group order follows Ctrl, OA, AZFc_Del, iNOA_B, iNOA_S, KS when present.",
            }
        ]
    )
    note.to_csv(OUT_DIR / "Fig7E_source_note.csv", index=False)
    generated.append("Fig7E_source_note.csv")

    return {
        "subpanel": "Dotplots (STs/Late_primary_SPCs/Round_spermatids)",
        "generated": generated,
        "raw": raw_inputs,
        "code": code_refs,
        "status": "available",
        "availability_note": "Original cell-level inputs and plotting script are available",
    }


def guess_if_role(file_name: str) -> str:
    name = file_name.lower()
    if "merge" in name:
        return "merge"
    if "c1x0" in name:
        return "channel_c1"
    if "c0x0" in name:
        return "channel_c0"
    if "c0-2" in name:
        return "channel_c0_alt"
    return "unknown"


def extract_ai_tif_links() -> list[str]:
    if not FIG7_AI.exists():
        return []
    proc = subprocess.run(
        ["strings", "-n", "8", str(FIG7_AI)],
        capture_output=True,
        text=True,
        check=False,
    )
    lines = proc.stdout.splitlines()

    links: set[str] = set()
    for i, line in enumerate(lines):
        if "/Users/xq/Desktop/" not in line:
            continue
        frag = line[line.find("/Users/xq/Desktop/") :].strip().replace('"', "")
        candidates = [frag]
        if i + 1 < len(lines):
            candidates.append((frag + lines[i + 1].strip()).replace('"', ""))

        for cand in candidates:
            cand = cand.replace("/Users/xq/Desktop//", "/Users/xq/Desktop/")
            m = re.search(r"/Users/xq/Desktop/[^/\n\r\"]+/[^/\n\r\"]+\.tif", cand)
            if m:
                links.add(m.group(0))

    return sorted(links)


def build_panel_f() -> dict:
    generated = []
    raw_inputs = [str(FIG7_AI)]

    ai_links = extract_ai_tif_links()
    ai_rows = []
    for link in ai_links:
        p = Path(link)
        relocated = Path("/Users/xq/Desktop/未命名文件夹") / p.parent.name / p.name
        ai_rows.append(
            {
                "ai_link_path": link,
                "exists_at_ai_path": p.exists(),
                "relocated_candidate": str(relocated),
                "exists_at_relocated_candidate": relocated.exists(),
            }
        )

    ai_df = pd.DataFrame(ai_rows)
    ai_df.to_csv(OUT_DIR / "Fig7F_ai_linked_tif_paths.csv", index=False)
    generated.append("Fig7F_ai_linked_tif_paths.csv")

    inventory_rows = []
    ai_by_name = {}
    for link in ai_links:
        ai_by_name.setdefault(Path(link).name, []).append(link)

    for group, dpath in PANEL_F_DIRS.items():
        for tif in sorted(dpath.glob("*.tif")) if dpath.exists() else []:
            matches = ai_by_name.get(tif.name, [])
            inventory_rows.append(
                {
                    "group": group,
                    "folder": str(dpath),
                    "file_name": tif.name,
                    "file_path": str(tif),
                    "file_role_guess": guess_if_role(tif.name),
                    "in_ai_link_by_name": bool(matches),
                    "ai_link_path_example": matches[0] if matches else "",
                    "exists": tif.exists(),
                }
            )

    inv_df = pd.DataFrame(inventory_rows)
    inv_df.to_csv(OUT_DIR / "Fig7F_image_path_inventory.csv", index=False)
    generated.append("Fig7F_image_path_inventory.csv")

    raw_inputs.extend([str(v) for v in PANEL_F_DIRS.values()])

    return {
        "subpanel": "IF image panel inventory (MCT1/SCP3)",
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../fig7_code/fig7_source_data.py"],
        "status": "partial",
        "availability_note": "AI-linked tif names and relocated local tif paths are inventoried",
    }


def build_panel_g() -> dict:
    note = pd.DataFrame(
        [
            {
                "panel": "G",
                "note": "Panel G is a conceptual schematic (Proposed metabolic uncoupling) manually composed in AI; no standalone numeric source table.",
            }
        ]
    )
    note.to_csv(OUT_DIR / "Fig7G_note.csv", index=False)

    return {
        "subpanel": "Working model schematic",
        "generated": ["Fig7G_note.csv"],
        "raw": [str(FIG7_AI), str(FIG7_PDF)],
        "code": ["none"],
        "status": "available",
        "availability_note": "Conceptual panel (manual drawing in AI)",
    }


def build_mapping_and_availability(panel_meta: dict[str, dict]) -> None:
    mapping_rows = []
    availability_rows = []

    panel_order = ["A", "B", "C", "D", "E", "F", "G"]
    for panel in panel_order:
        meta = panel_meta[panel]
        mapping_rows.append(
            {
                "panel": panel,
                "subpanel": meta["subpanel"],
                "generated_source_data": ";".join(meta["generated"]),
                "raw_input": ";".join(meta["raw"]),
                "code": ";".join(meta["code"]),
                "status": meta["status"],
            }
        )
        availability_rows.append(
            {
                "panel": panel,
                "status": meta["status"],
                "note": meta["availability_note"],
            }
        )

    pd.DataFrame(mapping_rows).to_csv(OUT_DIR / "Fig7_file_mapping.csv", index=False)
    pd.DataFrame(availability_rows).to_csv(OUT_DIR / "Fig7_panels_data_availability.csv", index=False)


def write_readme() -> None:
    csv_files = sorted(p.name for p in OUT_DIR.glob("*.csv"))
    readme = [
        "Figure 7 source data (Panels A-G)",
        "",
        "Generated files:",
    ]
    readme.extend([f"- {name}" for name in csv_files])
    readme.extend(
        [
            "",
            "Notes:",
            "- Panel A uses featureplot exports plus scMeta group-level diff tables (proxy numeric support).",
            "- Panel B pathway cartoon is AI-layout; provided source tables and scripts contain gene-level log2FC values.",
            "- Panel E dotplot metrics follow the original plotting script: pct_expressed and mean_log1p.",
            "- Panel F is an IF image path inventory; raw tifs remain in external Desktop folders.",
            "",
            "Build script:",
            "- ../fig7_code/fig7_source_data.py",
        ]
    )
    (OUT_DIR / "README.txt").write_text("\n".join(readme) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()

    panel_meta = {
        "A": build_panel_a(),
        "B": build_panel_b(),
        "C": build_panel_c(),
        "D": build_panel_d(),
        "E": build_panel_e(),
        "F": build_panel_f(),
        "G": build_panel_g(),
    }

    build_mapping_and_availability(panel_meta)
    write_readme()


if __name__ == "__main__":
    main()
