#!/usr/bin/env python3
from __future__ import annotations

import math
import shutil
import tarfile
from pathlib import Path
from typing import Iterable

import pandas as pd

BASE = Path(__file__).resolve().parents[1]  # .../source data
ROOT = BASE.parent
OUT_DIR = BASE / "fig6_source_data"
RAW_SUPPORT_DIR = OUT_DIR / "raw_support"

RAW_C = BASE / "fig6c.gz"
RAW_D = BASE / "fig6d.gz"
RAW_SUPPORT_C = OUT_DIR / "raw_support_fig6c"
RAW_SUPPORT_D = OUT_DIR / "raw_support_fig6d"

GO_CTRL_VS_DISEASE = BASE / "_fig6_probe_20260220" / "zips" / "LCa_ctrl_vs_disease" / "GO.csv"
GO_DISEASE_VS_CTRL = BASE / "_fig6_probe_20260220" / "zips" / "LCa_disease_vs_ctrl" / "GO.csv"

SIG_MEANS_CTRL = ROOT / "原来作图的数据" / "无精症" / "互作number of interaction" / "significant_means-ctrl.csv"
SIG_MEANS_DISEASE = ROOT / "原来作图的数据" / "无精症" / "互作number of interaction" / "significant_means-disease.csv"
CHAT_R = ROOT / "原来作图的数据" / "无精症" / "互作number of interaction" / "chat.R"

LCA_VOLCANO_GO_PDF = ROOT / "原来作图的数据" / "无精症" / "LCa火山图+Go富集.pdf"
OMIC_BUBBLE_PDF = ROOT / "原来作图的数据" / "无精症" / "OmicStudio_Bubble_Pro.pdf"

GROUP1_DIR = Path("/Users/xq/Desktop/group1_result.zip")
GROUP1_TSV = GROUP1_DIR / "DEGs_result_group1.tsv"
GROUP1_XLS = GROUP1_DIR / "DEGs_result_group1.xls"
GROUP1_VOLCANO_PDF = GROUP1_DIR / "DEGs_result_group1_volcano.pdf"
GROUP1_VOLCANO_SVG = GROUP1_DIR / "DEGs_result_group1_volcano.svg"
GROUP1_VOLCANO_PNG = GROUP1_DIR / "DEGs_result_group1_volcano.png"

LIGAND_CTRL_PDF = ROOT / "原来作图的数据" / "无精症" / "体细胞正常配体.pdf"
LIGAND_DISEASE_PDF = ROOT / "原来作图的数据" / "无精症" / "体细胞疾病配体.pdf"
RECEPTOR_CTRL_PDF = ROOT / "原来作图的数据" / "无精症" / "体细胞正常受体.pdf"
RECEPTOR_DISEASE_PDF = ROOT / "原来作图的数据" / "无精症" / "体细胞疾病受体.pdf"

LC_DIFF_DIRS = [
    ROOT / "原来作图的数据" / "汇总" / "20.LCs_STs_Myeloid_Lym_ECs_diff" / "LCs" / "Diff" / "P20073103",
    ROOT / "原来作图的数据" / "汇总" / "12.非生殖细胞组间差异分析" / "Diff" / "LCs",
]

FIG6F_TERMS = [
    "extracellular matrix structural constituent",
    "collagen-containing extracellular matrix",
    "focal adhesion",
    "cell-substrate junction",
]

FIG6G_TERMS = [
    "RNA catabolic process",
    "negative regulation of cell growth",
    "positive regulation of RNA splicing",
    "histone deacetylase binding",
]

IF_PATTERNS = [
    "/Users/xq/Desktop/未命名文件夹/KS-GAS1-GREEN-CYP11A1-RED-2/*.tif",
    "/Users/xq/Desktop/未命名文件夹/AZ-GAS1-GREEN-CYP11A1-RED-1/*.tif",
    "/Users/xq/Desktop/未命名文件夹/control-CYP11A1-red-GAS1-green-1-20x/*.tif",
    "/Users/xq/Desktop/未命名文件夹/control-CYP11A1-red-GAS1-green-1-40x/*.tif",
    "/Users/xq/Desktop/未命名文件夹/AZ-wt1-red-dhh-green-2-20x/*.tif",
    "/Users/xq/Desktop/未命名文件夹/AZ-wt1-red-dhh-green-2-40x/*.tif",
    "/Users/xq/Desktop/睾丸相关/未命名文件夹/DHH-WT1-1/merge*.tif",
    "/Users/xq/Desktop/睾丸相关/未命名文件夹/DHH-WT1-6*/merge.tif",
]


def ensure_extract(src: Path, dst: Path) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    sentinel = dst / "dotplot_export.csv"
    if sentinel.exists():
        return
    with tarfile.open(src, "r:gz") as tf:
        tf.extractall(dst)


def copy_to_raw_support(src: Path, out_name: str | None = None) -> str:
    RAW_SUPPORT_DIR.mkdir(parents=True, exist_ok=True)
    name = out_name if out_name else src.name
    dst = RAW_SUPPORT_DIR / name
    if src.exists():
        shutil.copy2(src, dst)
    return f"raw_support/{name}"


def safe_float(value) -> float | None:
    try:
        if pd.isna(value):
            return None
        return float(value)
    except Exception:
        return None


def neg_log10(value) -> float | None:
    v = safe_float(value)
    if v is None or v <= 0:
        return None
    return -math.log10(v)


def build_panel_a_edges() -> pd.DataFrame:
    rows = []
    ligand_edges = [
        ("Leydig Cells", "Immune Cells"),
        ("Leydig Cells", "Endothelial Cells"),
        ("Leydig Cells", "Sertoli Cells"),
        ("Leydig Cells", "Myoid Cells"),
    ]
    receptor_edges = [
        ("Immune Cells", "Leydig Cells"),
        ("Endothelial Cells", "Leydig Cells"),
        ("Sertoli Cells", "Leydig Cells"),
        ("Myoid Cells", "Leydig Cells"),
    ]
    for condition in ["Ctrl", "Azoospermia"]:
        for source, target in ligand_edges:
            rows.append(
                {
                    "condition": condition,
                    "network": "LC ligands",
                    "source": source,
                    "target": target,
                }
            )
        for source, target in receptor_edges:
            rows.append(
                {
                    "condition": condition,
                    "network": "LC receptors",
                    "source": source,
                    "target": target,
                }
            )
    return pd.DataFrame(rows)


def build_panel_ab_hedgehog_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    pair_cols = None
    for condition, path in [("Ctrl", SIG_MEANS_CTRL), ("Disease", SIG_MEANS_DISEASE)]:
        df = pd.read_csv(path)
        if pair_cols is None:
            pair_cols = [c for c in df.columns if "|" in c]
        subset = df[
            df["interacting_pair"].astype(str).str.contains(
                "DHH|PTCH|GAS1|CDON|BOC|HHIP|SHH", case=False, regex=True
            )
        ].copy()
        long_rows = []
        for _, row in subset.iterrows():
            for cell_pair in pair_cols:
                value = safe_float(row.get(cell_pair))
                if value is None:
                    continue
                long_rows.append(
                    {
                        "condition": condition,
                        "interacting_pair": row.get("interacting_pair", ""),
                        "cell_pair": cell_pair,
                        "mean_value": value,
                    }
                )
        frames.append(pd.DataFrame(long_rows))

    long_df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if long_df.empty:
        return long_df, long_df

    summary = (
        long_df.groupby(["condition", "interacting_pair"], as_index=False)
        .agg(
            n_edges=("cell_pair", "count"),
            mean_strength=("mean_value", "mean"),
            max_strength=("mean_value", "max"),
        )
        .sort_values(["condition", "n_edges", "mean_strength"], ascending=[True, False, False])
        .reset_index(drop=True)
    )
    return long_df, summary


def canonical_term_match(df: pd.DataFrame, terms: Iterable[str]) -> pd.DataFrame:
    rows = []
    description_lower = df["Description"].astype(str).str.lower()
    for idx, term in enumerate(terms, start=1):
        exact = df[description_lower == term.lower()]
        if exact.empty:
            exact = df[description_lower.str.contains(term.lower(), regex=False)]
        if exact.empty:
            continue
        row = exact.sort_values("pvalue", ascending=True).iloc[0].to_dict()
        row["term_order"] = idx
        rows.append(row)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("term_order").reset_index(drop=True)


def build_panel_fg_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    f_df = pd.read_csv(GO_CTRL_VS_DISEASE).copy()
    g_df = pd.read_csv(GO_DISEASE_VS_CTRL).copy()
    f_df["panel"] = "F"
    f_df["comparison"] = "LCa_ctrl_vs_disease"
    g_df["panel"] = "G"
    g_df["comparison"] = "LCa_disease_vs_ctrl"
    both = pd.concat([f_df, g_df], ignore_index=True)
    both["neg_log10_pvalue"] = both["pvalue"].apply(neg_log10)
    both = both.sort_values(["panel", "pvalue"]).reset_index(drop=True)

    f_sel = canonical_term_match(f_df, FIG6F_TERMS).copy()
    g_sel = canonical_term_match(g_df, FIG6G_TERMS).copy()
    for df, panel, comp in [(f_sel, "F", "LCa_ctrl_vs_disease"), (g_sel, "G", "LCa_disease_vs_ctrl")]:
        if not df.empty:
            df["panel"] = panel
            df["comparison"] = comp
            df["neg_log10_pvalue"] = df["pvalue"].apply(neg_log10)

    return both, f_sel, g_sel


def build_panel_e_group1_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    if not GROUP1_TSV.exists():
        return pd.DataFrame(), pd.DataFrame()

    raw = pd.read_csv(GROUP1_TSV, sep="\t")
    df = raw.rename(
        columns={
            "names": "gene",
            "scores": "scores",
            "logfoldchanges": "logfoldchanges",
            "pvals": "pvals",
            "pvals_adj": "pvals_adj",
            "pct_nz_group": "pct_nz_group",
            "pct_nz_reference": "pct_nz_reference",
            "mean": "mean_expression",
        }
    ).copy()
    if "pvals" in df.columns:
        df["neg_log10_pvals"] = df["pvals"].apply(neg_log10)
    if "pvals_adj" in df.columns:
        df["neg_log10_pvals_adj"] = df["pvals_adj"].apply(neg_log10)

    target_up, target_down = 124, 513
    summary_rows = []

    for lfc_cut in [0.25, 0.5]:
        summary_rows.append(
            {
                "mode": "fixed_cutoff",
                "p_column": "pvals",
                "p_cutoff": 0.05,
                "lfc_abs_cutoff": lfc_cut,
                "up_count": int(((df["logfoldchanges"] > lfc_cut) & (df["pvals"] < 0.05)).sum()),
                "down_count": int(((df["logfoldchanges"] < -lfc_cut) & (df["pvals"] < 0.05)).sum()),
                "target_up": target_up,
                "target_down": target_down,
            }
        )

    # Exact threshold scan against target counts for panel E labels.
    abs_cuts = sorted(df["logfoldchanges"].abs().dropna().round(6).unique().tolist())
    best = None
    exact_hits = []
    for cut in abs_cuts:
        up = int(((df["logfoldchanges"] > cut) & (df["pvals"] < 0.05)).sum())
        down = int(((df["logfoldchanges"] < -cut) & (df["pvals"] < 0.05)).sum())
        dist = abs(up - target_up) + abs(down - target_down)
        if best is None or dist < best[0]:
            best = (dist, float(cut), up, down)
        if up == target_up and down == target_down:
            exact_hits.append((float(cut), up, down))

    if exact_hits:
        exact_cut = exact_hits[0][0]
        summary_rows.append(
            {
                "mode": "exact_target_match",
                "p_column": "pvals",
                "p_cutoff": 0.05,
                "lfc_abs_cutoff": exact_cut,
                "up_count": target_up,
                "down_count": target_down,
                "target_up": target_up,
                "target_down": target_down,
            }
        )

    if best is not None:
        summary_rows.append(
            {
                "mode": "best_scan_match",
                "p_column": "pvals",
                "p_cutoff": 0.05,
                "lfc_abs_cutoff": best[1],
                "up_count": best[2],
                "down_count": best[3],
                "target_up": target_up,
                "target_down": target_down,
                "distance_to_target": best[0],
            }
        )

    # Also report adjusted p-value behavior (for transparency).
    summary_rows.append(
        {
            "mode": "adjp_reference",
            "p_column": "pvals_adj",
            "p_cutoff": 0.05,
            "lfc_abs_cutoff": 0.5,
            "up_count": int(((df["logfoldchanges"] > 0.5) & (df["pvals_adj"] < 0.05)).sum()),
            "down_count": int(((df["logfoldchanges"] < -0.5) & (df["pvals_adj"] < 0.05)).sum()),
            "target_up": target_up,
            "target_down": target_down,
        }
    )

    summary = pd.DataFrame(summary_rows)
    return df, summary


def pick_lfc_p_columns(df: pd.DataFrame) -> tuple[str | None, str | None]:
    lfc_candidates = ["avg_logFC", "avg_log2FC", "Log2FC", "log2FC", "logFC"]
    p_candidates = ["p_val_adj", "Adjusted p value", "p.adjust", "P_val_adj", "p_val", "P_val"]
    lfc = next((c for c in lfc_candidates if c in df.columns), None)
    pcol = next((c for c in p_candidates if c in df.columns), None)
    return lfc, pcol


def build_panel_e_candidate_table() -> pd.DataFrame:
    rows = []
    target_up, target_down = 124, 513
    for diff_dir in LC_DIFF_DIRS:
        if not diff_dir.exists():
            continue
        for path in sorted(diff_dir.glob("*.diffexpressed.xls")):
            df = pd.read_csv(path, sep="\t")
            lfc_col, p_col = pick_lfc_p_columns(df)
            if lfc_col is None or p_col is None:
                continue
            up = int(((df[lfc_col] > 0.5) & (df[p_col] < 0.05)).sum())
            down = int(((df[lfc_col] < -0.5) & (df[p_col] < 0.05)).sum())
            rows.append(
                {
                    "file": str(path.relative_to(ROOT)),
                    "comparison": path.name.replace("P20073103_", "").replace(".diffexpressed.xls", ""),
                    "lfc_column": lfc_col,
                    "p_column": p_col,
                    "up_log2fc_gt_0.5_p_lt_0.05": up,
                    "down_log2fc_lt_-0.5_p_lt_0.05": down,
                    "delta_vs_panelE_up124": up - target_up,
                    "delta_vs_panelE_down513": down - target_down,
                    "abs_total_delta": abs(up - target_up) + abs(down - target_down),
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.sort_values(["abs_total_delta", "file"]).drop_duplicates(
        subset=["comparison", "up_log2fc_gt_0.5_p_lt_0.05", "down_log2fc_lt_-0.5_p_lt_0.05"]
    )
    return out.reset_index(drop=True)


def build_panel_h_inventory() -> pd.DataFrame:
    rows = []
    for pattern in IF_PATTERNS:
        matches = sorted(Path("/").glob(pattern.lstrip("/")))
        rows.append(
            {
                "pattern": pattern,
                "matched_files": len(matches),
                "example_file": str(matches[0]) if matches else "",
            }
        )
    return pd.DataFrame(rows)


def write_csv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def build_mapping_tables(has_group1: bool) -> None:
    rows = [
        {
            "panel": "A",
            "subpanel": "LC ligands / LC receptors schematic network",
            "generated_source_data": "Fig6A_network_edges_manual.csv;Fig6AB_hedgehog_interaction_matrix_long.csv;Fig6AB_hedgehog_interaction_summary.csv",
            "raw_input": "raw_support/Fig6A_LC_ligands_ctrl.pdf;raw_support/Fig6A_LC_ligands_disease.pdf;raw_support/Fig6A_LC_receptors_ctrl.pdf;raw_support/Fig6A_LC_receptors_disease.pdf;raw_support/Fig6B_significant_means_ctrl.csv;raw_support/Fig6B_significant_means_disease.csv",
            "code": "../fig6_code/fig6_source_data.py",
            "status": "partial",
        },
        {
            "panel": "B",
            "subpanel": "Hedgehog signaling interaction bubble",
            "generated_source_data": "Fig6AB_hedgehog_interaction_matrix_long.csv;Fig6AB_hedgehog_interaction_summary.csv",
            "raw_input": "raw_support/Fig6B_significant_means_ctrl.csv;raw_support/Fig6B_significant_means_disease.csv;raw_support/Fig6B_chat_original.R;raw_support/Fig6B_OmicStudio_Bubble_Pro.pdf",
            "code": "../fig6_code/fig6_source_data.py",
            "status": "available",
        },
        {
            "panel": "C",
            "subpanel": "LC_a dotplot",
            "generated_source_data": "Fig6C_dotplot_LCa.csv;Fig6C_heatmap_LCa.csv;Fig6C_violin_cells_LCa.csv;Fig6CD_source_note.csv",
            "raw_input": "../fig6c.gz;raw_support_fig6c/dotplot_export.csv;raw_support_fig6c/heatmap_export.csv;raw_support_fig6c/violin_export.csv",
            "code": "../fig6_code/fig6_cd_prepare_source_data.py;../fig6_code/fig6_cd_dotplot.py;../fig6_code/fig6_source_data.py",
            "status": "available",
        },
        {
            "panel": "D",
            "subpanel": "ST dotplot",
            "generated_source_data": "Fig6D_dotplot_ST_no_OA.csv;Fig6D_dotplot_ST_full_with_OA.csv;Fig6D_heatmap_ST_no_OA.csv;Fig6D_heatmap_ST_full_with_OA.csv;Fig6D_violin_cells_ST.csv;Fig6CD_source_note.csv",
            "raw_input": "../fig6d.gz;raw_support_fig6d/dotplot_export.csv;raw_support_fig6d/heatmap_export.csv;raw_support_fig6d/violin_export.csv",
            "code": "../fig6_code/fig6_cd_prepare_source_data.py;../fig6_code/fig6_cd_dotplot.py;../fig6_code/fig6_source_data.py",
            "status": "available",
        },
        {
            "panel": "E",
            "subpanel": "LC_a volcano",
            "generated_source_data": (
                "Fig6E_volcano_group1_all.csv;Fig6E_volcano_group1_threshold_summary.csv;Fig6E_note.csv;Fig6E_volcano_candidate_counts.csv"
                if has_group1
                else "Fig6E_volcano_candidate_counts.csv;Fig6E_note.csv"
            ),
            "raw_input": (
                "raw_support/Fig6E_group1_DEGs_result.tsv;raw_support/Fig6E_group1_DEGs_result.xls;raw_support/Fig6E_group1_volcano.pdf;raw_support/Fig6E_group1_volcano.svg;raw_support/Fig6E_LCa_volcano_go_panel.pdf"
                if has_group1
                else "raw_support/Fig6E_LCa_volcano_go_panel.pdf;raw_support/Fig6E_candidates_P20073103_AZFc_DelvsCtrl.diffexpressed.xls;raw_support/Fig6E_candidates_P20073103_iNOA_SvsCtrl.diffexpressed.xls;raw_support/Fig6E_candidates_P20073103_iNOA_BvsCtrl.diffexpressed.xls;raw_support/Fig6E_candidates_P20073103_KSvsCtrl.diffexpressed.xls;raw_support/Fig6E_candidates_P20073103_OAvsCtrl.diffexpressed.xls"
            ),
            "code": "../fig6_code/fig6_source_data.py",
            "status": "available" if has_group1 else "partial",
        },
        {
            "panel": "F",
            "subpanel": "GO terms enriched in azoospermia",
            "generated_source_data": "Fig6FG_GO_all_long.csv;Fig6F_GO_terms_selected.csv",
            "raw_input": "raw_support/Fig6F_GO_LCa_ctrl_vs_disease.csv",
            "code": "../fig6_code/fig6_source_data.py",
            "status": "available",
        },
        {
            "panel": "G",
            "subpanel": "GO terms enriched in ctrl",
            "generated_source_data": "Fig6FG_GO_all_long.csv;Fig6G_GO_terms_selected.csv",
            "raw_input": "raw_support/Fig6G_GO_LCa_disease_vs_ctrl.csv",
            "code": "../fig6_code/fig6_source_data.py",
            "status": "available",
        },
        {
            "panel": "H",
            "subpanel": "IF image panel inventory",
            "generated_source_data": "Fig6H_image_path_inventory.csv",
            "raw_input": "external Desktop tif paths from Figure 6 AI links (see Fig6H_image_path_inventory.csv)",
            "code": "../fig6_code/fig6_source_data.py",
            "status": "partial",
        },
        {
            "panel": "I",
            "subpanel": "Working model schematic",
            "generated_source_data": "Fig6I_note.csv",
            "raw_input": "none (manual schematic in AI)",
            "code": "none",
            "status": "available",
        },
    ]
    mapping = pd.DataFrame(rows)
    write_csv(OUT_DIR / "Fig6_file_mapping.csv", mapping)

    availability = pd.DataFrame(
        [
            {"panel": r["panel"], "status": r["status"], "note": r["subpanel"]}
            for r in rows
        ]
    )
    write_csv(OUT_DIR / "Fig6_panels_data_availability.csv", availability)


def write_readme(has_group1: bool) -> None:
    text = """Figure 6 source data (Panels A-I)

Generated files:
- Fig6A_network_edges_manual.csv
- Fig6AB_hedgehog_interaction_matrix_long.csv
- Fig6AB_hedgehog_interaction_summary.csv
- Fig6C_dotplot_LCa.csv
- Fig6C_heatmap_LCa.csv
- Fig6C_violin_cells_LCa.csv
- Fig6D_dotplot_ST_full_with_OA.csv
- Fig6D_dotplot_ST_no_OA.csv
- Fig6D_heatmap_ST_full_with_OA.csv
- Fig6D_heatmap_ST_no_OA.csv
- Fig6D_violin_cells_ST.csv
- Fig6CD_source_note.csv
- Fig6E_volcano_group1_all.csv
- Fig6E_volcano_group1_threshold_summary.csv
- Fig6E_volcano_candidate_counts.csv
- Fig6E_note.csv
- Fig6FG_GO_all_long.csv
- Fig6F_GO_terms_selected.csv
- Fig6G_GO_terms_selected.csv
- Fig6H_image_path_inventory.csv
- Fig6I_note.csv
- Fig6_file_mapping.csv
- Fig6_panels_data_availability.csv

Notes:
- Panels C/D now use original archives: ../fig6c.gz and ../fig6d.gz.
- Panel E uses group1_result DEG table when available; threshold summary is in Fig6E_volcano_group1_threshold_summary.csv.
- If group1_result is unavailable, fallback candidate LCs diff tables are summarized in Fig6E_volcano_candidate_counts.csv.
- Panels F/G selected term bars are extracted from LCa_ctrl_vs_disease and LCa_disease_vs_ctrl GO.csv files.

Build scripts:
- ../fig6_code/fig6_cd_prepare_source_data.py
- ../fig6_code/fig6_cd_dotplot.py
- ../fig6_code/fig6_source_data.py
"""
    if not has_group1:
        text = text.replace(
            "- Fig6E_volcano_group1_all.csv\n- Fig6E_volcano_group1_threshold_summary.csv\n",
            "",
        )
    (OUT_DIR / "README.txt").write_text(text, encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RAW_SUPPORT_DIR.mkdir(parents=True, exist_ok=True)

    ensure_extract(RAW_C, RAW_SUPPORT_C)
    ensure_extract(RAW_D, RAW_SUPPORT_D)

    # Reuse existing C/D builder to keep output compatibility.
    from fig6_cd_prepare_source_data import main as build_cd_main

    build_cd_main()

    # Keep local copies of key raw inputs.
    copy_to_raw_support(GO_CTRL_VS_DISEASE, "Fig6F_GO_LCa_ctrl_vs_disease.csv")
    copy_to_raw_support(GO_DISEASE_VS_CTRL, "Fig6G_GO_LCa_disease_vs_ctrl.csv")
    copy_to_raw_support(SIG_MEANS_CTRL, "Fig6B_significant_means_ctrl.csv")
    copy_to_raw_support(SIG_MEANS_DISEASE, "Fig6B_significant_means_disease.csv")
    copy_to_raw_support(CHAT_R, "Fig6B_chat_original.R")
    copy_to_raw_support(LCA_VOLCANO_GO_PDF, "Fig6E_LCa_volcano_go_panel.pdf")
    copy_to_raw_support(OMIC_BUBBLE_PDF, "Fig6B_OmicStudio_Bubble_Pro.pdf")
    copy_to_raw_support(LIGAND_CTRL_PDF, "Fig6A_LC_ligands_ctrl.pdf")
    copy_to_raw_support(LIGAND_DISEASE_PDF, "Fig6A_LC_ligands_disease.pdf")
    copy_to_raw_support(RECEPTOR_CTRL_PDF, "Fig6A_LC_receptors_ctrl.pdf")
    copy_to_raw_support(RECEPTOR_DISEASE_PDF, "Fig6A_LC_receptors_disease.pdf")
    has_group1 = GROUP1_TSV.exists()
    if has_group1:
        copy_to_raw_support(GROUP1_TSV, "Fig6E_group1_DEGs_result.tsv")
        copy_to_raw_support(GROUP1_XLS, "Fig6E_group1_DEGs_result.xls")
        copy_to_raw_support(GROUP1_VOLCANO_PDF, "Fig6E_group1_volcano.pdf")
        copy_to_raw_support(GROUP1_VOLCANO_SVG, "Fig6E_group1_volcano.svg")
        copy_to_raw_support(GROUP1_VOLCANO_PNG, "Fig6E_group1_volcano.png")

    # Panel A
    a_edges = build_panel_a_edges()
    write_csv(OUT_DIR / "Fig6A_network_edges_manual.csv", a_edges)

    # Panel A/B
    ab_long, ab_summary = build_panel_ab_hedgehog_tables()
    write_csv(OUT_DIR / "Fig6AB_hedgehog_interaction_matrix_long.csv", ab_long)
    write_csv(OUT_DIR / "Fig6AB_hedgehog_interaction_summary.csv", ab_summary)

    # Panel E
    e_group1, e_group1_summary = build_panel_e_group1_tables()
    if not e_group1.empty:
        write_csv(OUT_DIR / "Fig6E_volcano_group1_all.csv", e_group1)
    if not e_group1_summary.empty:
        write_csv(OUT_DIR / "Fig6E_volcano_group1_threshold_summary.csv", e_group1_summary)

    e_candidates = build_panel_e_candidate_table()
    write_csv(OUT_DIR / "Fig6E_volcano_candidate_counts.csv", e_candidates)

    if not e_group1_summary.empty and (e_group1_summary["mode"] == "exact_target_match").any():
        exact = e_group1_summary[e_group1_summary["mode"] == "exact_target_match"].iloc[0]
        e_note = pd.DataFrame(
            [
                {
                    "target_panelE": "Sig_Up(124), Sig_Down(513)",
                    "used_file": str(GROUP1_TSV),
                    "p_column": "pvals",
                    "p_cutoff": 0.05,
                    "lfc_column": "logfoldchanges",
                    "lfc_abs_cutoff": float(exact["lfc_abs_cutoff"]),
                    "up_count": int(exact["up_count"]),
                    "down_count": int(exact["down_count"]),
                    "comment": "Exact target matched from group1_result DEG table.",
                }
            ]
        )
    elif not e_candidates.empty:
        best = e_candidates.iloc[0]
        e_note = pd.DataFrame(
            [
                {
                    "target_panelE": "Sig_Up(124), Sig_Down(513)",
                    "used_file": best["file"],
                    "p_column": best["p_column"],
                    "p_cutoff": 0.05,
                    "lfc_column": best["lfc_column"],
                    "lfc_abs_cutoff": 0.5,
                    "up_count": int(best["up_log2fc_gt_0.5_p_lt_0.05"]),
                    "down_count": int(best["down_log2fc_lt_-0.5_p_lt_0.05"]),
                    "comment": "Fallback candidate from current LCs diffexpressed tables (no exact match).",
                }
            ]
        )
    else:
        e_note = pd.DataFrame(
            [
                {
                    "target_panelE": "Sig_Up(124), Sig_Down(513)",
                    "used_file": "",
                    "p_column": "",
                    "p_cutoff": "",
                    "lfc_column": "",
                    "lfc_abs_cutoff": "",
                    "up_count": "",
                    "down_count": "",
                    "comment": "No readable candidate diffexpressed table found.",
                }
            ]
        )
    write_csv(OUT_DIR / "Fig6E_note.csv", e_note)

    # Copy top candidate raw files for panel E for traceability.
    for path in [
        ROOT / "原来作图的数据" / "汇总" / "20.LCs_STs_Myeloid_Lym_ECs_diff" / "LCs" / "Diff" / "P20073103" / "P20073103_AZFc_DelvsCtrl.diffexpressed.xls",
        ROOT / "原来作图的数据" / "汇总" / "20.LCs_STs_Myeloid_Lym_ECs_diff" / "LCs" / "Diff" / "P20073103" / "P20073103_iNOA_SvsCtrl.diffexpressed.xls",
        ROOT / "原来作图的数据" / "汇总" / "20.LCs_STs_Myeloid_Lym_ECs_diff" / "LCs" / "Diff" / "P20073103" / "P20073103_iNOA_BvsCtrl.diffexpressed.xls",
        ROOT / "原来作图的数据" / "汇总" / "20.LCs_STs_Myeloid_Lym_ECs_diff" / "LCs" / "Diff" / "P20073103" / "P20073103_KSvsCtrl.diffexpressed.xls",
        ROOT / "原来作图的数据" / "汇总" / "20.LCs_STs_Myeloid_Lym_ECs_diff" / "LCs" / "Diff" / "P20073103" / "P20073103_OAvsCtrl.diffexpressed.xls",
    ]:
        copy_to_raw_support(path, f"Fig6E_candidates_{path.name}")

    # Panel F/G
    fg_all, f_selected, g_selected = build_panel_fg_tables()
    write_csv(OUT_DIR / "Fig6FG_GO_all_long.csv", fg_all)
    write_csv(OUT_DIR / "Fig6F_GO_terms_selected.csv", f_selected)
    write_csv(OUT_DIR / "Fig6G_GO_terms_selected.csv", g_selected)

    # Panel H
    h_inventory = build_panel_h_inventory()
    write_csv(OUT_DIR / "Fig6H_image_path_inventory.csv", h_inventory)

    # Panel I
    i_note = pd.DataFrame(
        [
            {
                "panel": "I",
                "type": "schematic",
                "raw_numeric_source": "none",
                "note": "Conceptual model manually composed in Figure 6 AI.",
            }
        ]
    )
    write_csv(OUT_DIR / "Fig6I_note.csv", i_note)

    build_mapping_tables(has_group1=has_group1)
    write_readme(has_group1=has_group1)


if __name__ == "__main__":
    main()
