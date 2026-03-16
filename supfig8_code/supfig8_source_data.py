#!/usr/bin/env python3
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

BASE = Path(__file__).resolve().parents[1]  # .../source data
ROOT = BASE.parent
OUT_DIR = BASE / "supfig8_source_data"
RAW_SUPPORT_DIR = OUT_DIR / "raw_support"

SUPFIG8_AI = ROOT / "Supplementary Figure 8.ai"
SUPFIG8_PDF = ROOT / "Supplementary Figure 8.pdf"
SUPFIG8_PNG = ROOT / "Supplementary Figure 8.png"

PACKAGE_DIR = ROOT / "原来作图的数据" / "Integrated_NicheAxis_NCstyle_Package_v2"
NC_FIGS_DIR = PACKAGE_DIR / "NCstyle_figs"
SCORES_XLSX = PACKAGE_DIR / "DonorLevel_IntegratedScores_withDHH_v2.xlsx"

PANEL_A_FILES = [
    PACKAGE_DIR / "Fig9A_SomaticScores_heatmap_RdBu_groupbar_orderPC1v2_framed_v2.pdf",
    PACKAGE_DIR / "Fig9A_SomaticScores_heatmap_groupedWithinPC1v2_v2.pdf",
    NC_FIGS_DIR / "Fig9A_SomaticScores_heatmap_groupbar_orderPC1v2.pdf",
    NC_FIGS_DIR / "Fig9A_SomaticScores_heatmap_groupbar_orderPC1v2.png",
    NC_FIGS_DIR / "Fig9A_SomaticScores_heatmap_groupbar_orderPC1v2.svg",
]

PANEL_B_FILES = [
    NC_FIGS_DIR / "FigS_PCA_PC1v2_PC2v2_colored.pdf",
    NC_FIGS_DIR / "FigS_PCA_PC1v2_PC2v2_colored.png",
    NC_FIGS_DIR / "FigS_PCA_PC1v2_PC2v2_colored_labeled.pdf",
    NC_FIGS_DIR / "FigS_PCA_PC1v2_PC2v2_colored_labeled.png",
]

PANEL_C_FILES = [
    NC_FIGS_DIR / "Fig9B_PC1v2_vs_GermMaturity_colored.pdf",
    NC_FIGS_DIR / "Fig9B_PC1v2_vs_GermMaturity_colored.png",
]

PANEL_D_FILES = [
    NC_FIGS_DIR / "Fig9C_PC1v2_vs_GermLatestStage_colored.pdf",
    NC_FIGS_DIR / "Fig9C_PC1v2_vs_GermLatestStage_colored.png",
]

PANEL_E_FILES = [
    PACKAGE_DIR / "Fig9H_ScoreCorrelations_Spearman_RdBu_allCells_stars_noHh.pdf",
    PACKAGE_DIR / "Fig9H_ScoreCorrelations_Spearman_RdBu_FDRsig_noHh.pdf",
    NC_FIGS_DIR / "FigS_ScoreCorrelation_Spearman_NCstyle.pdf",
    NC_FIGS_DIR / "FigS_ScoreCorrelation_Spearman_NCstyle.png",
]

GROUP_ORDER = ["Ctrl", "OA", "AZFc_Del", "iNOA_B", "iNOA_S", "KS"]

A_SCORE_MAP = [
    ("ST Maturity", "ST_Maturity"),
    ("BTB Integrity", "BTB_Integrity"),
    ("ST Stress IEG", "ST_Stress_IEG"),
    ("ST SASP", "ST_SASP"),
    ("ST Lactate Export", "ST_Lactate_Export"),
    ("ST DHH Ligand", "ST_DHH_Ligand"),
    ("LC Mature", "LC_Mature"),
    ("LC Immature", "LC_Immature"),
    ("Hh Response", "Hh_Response"),
]

E_SCORE_MAP = [
    ("ST Maturity", "ST_Maturity"),
    ("BTB Integrity", "BTB_Integrity"),
    ("ST Stress IEG", "ST_Stress_IEG"),
    ("ST SASP", "ST_SASP"),
    ("ST Lactate Export", "ST_Lactate_Export"),
    ("ST DHH Ligand", "ST_DHH_Ligand"),
    ("LC Mature", "LC_Mature"),
    ("LC Immature", "LC_Immature"),
    ("Germ Maturity Index", "Germ_Maturity_Index"),
    ("Germ Latest Stage", "Germ_Latest_Stage"),
]


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


def zscore_series(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    mu = x.mean(skipna=True)
    sd = x.std(skipna=True, ddof=0)
    if pd.isna(sd) or float(sd) == 0.0:
        return pd.Series(np.zeros(len(x)), index=s.index)
    return (x - mu) / sd


def star_from_q(q: float | None) -> str:
    if q is None or pd.isna(q):
        return ""
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


def add_main_and_shared_raw() -> list[str]:
    raw: list[str] = []
    for src, alias in [
        (SUPFIG8_AI, "SupFig8_main.ai"),
        (SUPFIG8_PDF, "SupFig8_main.pdf"),
        (SUPFIG8_PNG, "SupFig8_main.png"),
        (SCORES_XLSX, "SupFig8_DonorLevel_IntegratedScores_withDHH_v2.xlsx"),
    ]:
        rel = copy_to_raw_support(src, alias)
        if rel:
            raw.append(rel)
    return raw


def prepare_base_df() -> pd.DataFrame:
    df = pd.read_excel(SCORES_XLSX)
    df["Group"] = pd.Categorical(df["Group"], categories=GROUP_ORDER, ordered=True)
    return df


def build_panel_a(df: pd.DataFrame, shared_raw: list[str]) -> dict:
    generated: list[str] = []
    raw = list(shared_raw)
    for src in PANEL_A_FILES:
        rel = copy_to_raw_support(src, f"SupFig8A_{src.name}")
        if rel:
            raw.append(rel)

    sample_order = (
        df[["Sample_ID", "Group", "PC1_v2", "PC2_v2"]]
        .sort_values("PC1_v2", ascending=True)
        .reset_index(drop=True)
    )
    sample_order["sample_rank_pc1v2"] = np.arange(1, len(sample_order) + 1)
    sample_order.to_csv(OUT_DIR / "SupFig8A_sample_order_by_PC1v2.csv", index=False)
    generated.append("SupFig8A_sample_order_by_PC1v2.csv")

    row_order = [x[0] for x in A_SCORE_MAP]
    row_df = pd.DataFrame(
        {
            "row_order": np.arange(1, len(row_order) + 1),
            "score_display": row_order,
            "score_column": [x[1] for x in A_SCORE_MAP],
        }
    )
    row_df.to_csv(OUT_DIR / "SupFig8A_score_row_order.csv", index=False)
    generated.append("SupFig8A_score_row_order.csv")

    df_ord = sample_order[["Sample_ID"]].merge(df, on="Sample_ID", how="left")

    long_rows: list[dict] = []
    z_wide = pd.DataFrame(index=row_order, columns=df_ord["Sample_ID"].tolist(), dtype=float)
    raw_wide = pd.DataFrame(index=row_order, columns=df_ord["Sample_ID"].tolist(), dtype=float)

    for disp, col in A_SCORE_MAP:
        vals = pd.to_numeric(df_ord[col], errors="coerce")
        zvals = zscore_series(vals)
        raw_wide.loc[disp, :] = vals.values
        z_wide.loc[disp, :] = zvals.values

        for sid, grp, raw_v, z_v in zip(df_ord["Sample_ID"], df_ord["Group"], vals, zvals):
            long_rows.append(
                {
                    "Sample_ID": sid,
                    "Group": grp,
                    "score_display": disp,
                    "score_column": col,
                    "value_raw": raw_v,
                    "value_zscore_rowwise": z_v,
                }
            )

    long_df = pd.DataFrame(long_rows)
    long_df.to_csv(OUT_DIR / "SupFig8A_somatic_scores_long.csv", index=False)
    generated.append("SupFig8A_somatic_scores_long.csv")

    z_wide.insert(0, "score_display", z_wide.index)
    z_wide.to_csv(OUT_DIR / "SupFig8A_somatic_scores_heatmap_z_matrix.csv", index=False)
    generated.append("SupFig8A_somatic_scores_heatmap_z_matrix.csv")

    raw_wide.insert(0, "score_display", raw_wide.index)
    raw_wide.to_csv(OUT_DIR / "SupFig8A_somatic_scores_heatmap_raw_matrix.csv", index=False)
    generated.append("SupFig8A_somatic_scores_heatmap_raw_matrix.csv")

    group_bar = sample_order[["Sample_ID", "Group", "sample_rank_pc1v2"]].copy()
    group_bar.to_csv(OUT_DIR / "SupFig8A_group_bar_annotation.csv", index=False)
    generated.append("SupFig8A_group_bar_annotation.csv")

    note = pd.DataFrame(
        [
            {
                "panel": "A",
                "note": "Heatmap uses donor-level scores row-wise z-scored across all donors; samples ordered by PC1_v2 ascending.",
            }
        ]
    )
    note.to_csv(OUT_DIR / "SupFig8A_source_note.csv", index=False)
    generated.append("SupFig8A_source_note.csv")

    return {
        "subpanel": "Somatic niche score heatmap with group bar",
        "generated": generated,
        "raw": raw,
        "code": ["../supfig8_code/supfig8_source_data.py", "../supfig8_code/supfig8_plot_panels.py"],
        "status": "available",
        "availability_note": "Direct donor-level score table and matching heatmap exports are available.",
    }


def build_panel_b(df: pd.DataFrame, shared_raw: list[str]) -> dict:
    generated: list[str] = []
    raw = list(shared_raw)
    for src in PANEL_B_FILES:
        rel = copy_to_raw_support(src, f"SupFig8B_{src.name}")
        if rel:
            raw.append(rel)

    out = df[["Sample_ID", "Group", "PC1", "PC2", "PC1_v2", "PC2_v2"]].copy()
    out = out.sort_values("Sample_ID").reset_index(drop=True)
    out.to_csv(OUT_DIR / "SupFig8B_pca_points.csv", index=False)
    generated.append("SupFig8B_pca_points.csv")

    note = pd.DataFrame(
        [
            {
                "panel": "B",
                "note": "PCA scatter in Supplementary Figure 8 uses PC1_v2 and PC2_v2 coordinates from donor-level integrated scores.",
            }
        ]
    )
    note.to_csv(OUT_DIR / "SupFig8B_source_note.csv", index=False)
    generated.append("SupFig8B_source_note.csv")

    return {
        "subpanel": "PCA scatter (PC1 vs PC2)",
        "generated": generated,
        "raw": raw,
        "code": ["../supfig8_code/supfig8_source_data.py", "../supfig8_code/supfig8_plot_panels.py"],
        "status": "available",
        "availability_note": "Exact donor-level PCA coordinates are available.",
    }


def build_panel_c(df: pd.DataFrame, shared_raw: list[str]) -> dict:
    generated: list[str] = []
    raw = list(shared_raw)
    for src in PANEL_C_FILES:
        rel = copy_to_raw_support(src, f"SupFig8C_{src.name}")
        if rel:
            raw.append(rel)

    out = df[["Sample_ID", "Group", "PC1_v2", "Germ_Maturity_Index"]].copy()
    out = out.rename(columns={"PC1_v2": "x_niche_collapse_axis", "Germ_Maturity_Index": "y_germ_maturity_index"})
    out.to_csv(OUT_DIR / "SupFig8C_niche_vs_germ_maturity.csv", index=False)
    generated.append("SupFig8C_niche_vs_germ_maturity.csv")

    rho, pval = spearmanr(out["x_niche_collapse_axis"], out["y_germ_maturity_index"], nan_policy="omit")
    fit = np.polyfit(out["x_niche_collapse_axis"], out["y_germ_maturity_index"], 1)
    summary = pd.DataFrame(
        [
            {
                "panel": "C",
                "x": "Niche collapse axis (PC1v2)",
                "y": "Germ maturity index",
                "n_donors": int(out.shape[0]),
                "spearman_rho": float(rho),
                "p_value": float(pval),
                "line_slope": float(fit[0]),
                "line_intercept": float(fit[1]),
            }
        ]
    )
    summary.to_csv(OUT_DIR / "SupFig8C_spearman_summary.csv", index=False)
    generated.append("SupFig8C_spearman_summary.csv")

    return {
        "subpanel": "PC1v2 vs germ maturity index",
        "generated": generated,
        "raw": raw,
        "code": ["../supfig8_code/supfig8_source_data.py", "../supfig8_code/supfig8_plot_panels.py"],
        "status": "available",
        "availability_note": "Exact rho and p-value can be recomputed from donor-level table.",
    }


def build_panel_d(df: pd.DataFrame, shared_raw: list[str]) -> dict:
    generated: list[str] = []
    raw = list(shared_raw)
    for src in PANEL_D_FILES:
        rel = copy_to_raw_support(src, f"SupFig8D_{src.name}")
        if rel:
            raw.append(rel)

    out = df[["Sample_ID", "Group", "PC1_v2", "Germ_Latest_Stage"]].copy()
    out = out.rename(columns={"PC1_v2": "x_niche_collapse_axis", "Germ_Latest_Stage": "y_germline_latest_stage_0to7"})
    out.to_csv(OUT_DIR / "SupFig8D_niche_vs_germ_latest_stage.csv", index=False)
    generated.append("SupFig8D_niche_vs_germ_latest_stage.csv")

    rho, pval = spearmanr(out["x_niche_collapse_axis"], out["y_germline_latest_stage_0to7"], nan_policy="omit")
    fit = np.polyfit(out["x_niche_collapse_axis"], out["y_germline_latest_stage_0to7"], 1)
    summary = pd.DataFrame(
        [
            {
                "panel": "D",
                "x": "Niche collapse axis (PC1v2)",
                "y": "Germline latest stage (0-7)",
                "n_donors": int(out.shape[0]),
                "spearman_rho": float(rho),
                "p_value": float(pval),
                "line_slope": float(fit[0]),
                "line_intercept": float(fit[1]),
            }
        ]
    )
    summary.to_csv(OUT_DIR / "SupFig8D_spearman_summary.csv", index=False)
    generated.append("SupFig8D_spearman_summary.csv")

    return {
        "subpanel": "PC1v2 vs germline latest stage",
        "generated": generated,
        "raw": raw,
        "code": ["../supfig8_code/supfig8_source_data.py", "../supfig8_code/supfig8_plot_panels.py"],
        "status": "available",
        "availability_note": "Exact rho and p-value can be recomputed from donor-level table.",
    }


def build_panel_e(df: pd.DataFrame, shared_raw: list[str]) -> dict:
    generated: list[str] = []
    raw = list(shared_raw)
    for src in PANEL_E_FILES:
        rel = copy_to_raw_support(src, f"SupFig8E_{src.name}")
        if rel:
            raw.append(rel)

    score_names = [x[0] for x in E_SCORE_MAP]
    cols = [x[1] for x in E_SCORE_MAP]
    work = df[cols].copy()
    work.columns = score_names

    n = len(score_names)
    corr = np.full((n, n), np.nan)
    pval = np.full((n, n), np.nan)

    for i in range(n):
        for j in range(n):
            rho, pv = spearmanr(work.iloc[:, i], work.iloc[:, j], nan_policy="omit")
            corr[i, j] = rho
            pval[i, j] = pv

    corr_df = pd.DataFrame(corr, index=score_names, columns=score_names)
    pval_df = pd.DataFrame(pval, index=score_names, columns=score_names)

    tri = np.triu_indices(n, 1)
    _, qvals, _, _ = multipletests(pval[tri], method="fdr_bh")

    qmat = np.full((n, n), np.nan)
    qmat[tri] = qvals
    qmat[(tri[1], tri[0])] = qvals
    np.fill_diagonal(qmat, 0.0)
    q_df = pd.DataFrame(qmat, index=score_names, columns=score_names)

    stars = np.full((n, n), "", dtype=object)
    for i in range(n):
        for j in range(n):
            if i == j:
                stars[i, j] = ""
            else:
                stars[i, j] = star_from_q(qmat[i, j])
    stars_df = pd.DataFrame(stars, index=score_names, columns=score_names)

    text = np.empty((n, n), dtype=object)
    for i in range(n):
        for j in range(n):
            value = corr[i, j]
            if np.isnan(value):
                text[i, j] = ""
            elif i == j:
                text[i, j] = f"{value:.2f}"
            else:
                text[i, j] = f"{value:.2f}{stars[i, j]}"
    text_df = pd.DataFrame(text, index=score_names, columns=score_names)

    corr_df.to_csv(OUT_DIR / "SupFig8E_score_correlation_spearman_matrix.csv")
    generated.append("SupFig8E_score_correlation_spearman_matrix.csv")

    pval_df.to_csv(OUT_DIR / "SupFig8E_score_correlation_pvalue_matrix.csv")
    generated.append("SupFig8E_score_correlation_pvalue_matrix.csv")

    q_df.to_csv(OUT_DIR / "SupFig8E_score_correlation_fdr_bh_matrix.csv")
    generated.append("SupFig8E_score_correlation_fdr_bh_matrix.csv")

    stars_df.to_csv(OUT_DIR / "SupFig8E_score_correlation_significance_stars_matrix.csv")
    generated.append("SupFig8E_score_correlation_significance_stars_matrix.csv")

    text_df.to_csv(OUT_DIR / "SupFig8E_score_correlation_text_matrix.csv")
    generated.append("SupFig8E_score_correlation_text_matrix.csv")

    long_rows: list[dict] = []
    for i, row_name in enumerate(score_names):
        for j, col_name in enumerate(score_names):
            long_rows.append(
                {
                    "row_score": row_name,
                    "col_score": col_name,
                    "spearman_rho": corr[i, j],
                    "p_value": pval[i, j],
                    "fdr_bh_q": qmat[i, j],
                    "stars": stars[i, j],
                    "text_label": text[i, j],
                }
            )
    long_df = pd.DataFrame(long_rows)
    long_df.to_csv(OUT_DIR / "SupFig8E_score_correlation_long.csv", index=False)
    generated.append("SupFig8E_score_correlation_long.csv")

    legend = pd.DataFrame(
        [
            {"symbol": "*", "rule": "q < 0.05 (Benjamini-Hochberg FDR)"},
            {"symbol": "**", "rule": "q < 0.01 (Benjamini-Hochberg FDR)"},
            {"symbol": "***", "rule": "q < 0.001 (Benjamini-Hochberg FDR)"},
        ]
    )
    legend.to_csv(OUT_DIR / "SupFig8E_significance_legend.csv", index=False)
    generated.append("SupFig8E_significance_legend.csv")

    return {
        "subpanel": "Donor-level score correlation heatmap",
        "generated": generated,
        "raw": raw,
        "code": ["../supfig8_code/supfig8_source_data.py", "../supfig8_code/supfig8_plot_panels.py"],
        "status": "available",
        "availability_note": "Spearman matrix and BH-FDR significance stars are directly reproducible from donor-level score table.",
    }


def write_file_mapping(panel_map: dict[str, dict]) -> None:
    rows = []
    for panel in ["A", "B", "C", "D", "E"]:
        info = panel_map[panel]
        rows.append(
            {
                "panel": panel,
                "subpanel": info["subpanel"],
                "generated_source_data": ";".join(info["generated"]),
                "raw_input": ";".join(info["raw"]),
                "code": ";".join(info["code"]),
                "status": info["status"],
            }
        )
    pd.DataFrame(rows).to_csv(OUT_DIR / "SupFig8_file_mapping.csv", index=False)


def write_availability(panel_map: dict[str, dict]) -> None:
    rows = []
    for panel in ["A", "B", "C", "D", "E"]:
        rows.append({"panel": panel, "status": panel_map[panel]["status"], "note": panel_map[panel]["availability_note"]})
    pd.DataFrame(rows).to_csv(OUT_DIR / "SupFig8_panels_data_availability.csv", index=False)


def write_readme() -> None:
    files = sorted(p.name for p in OUT_DIR.iterdir() if p.is_file())
    lines = [
        "Supplementary Figure 8 source data (Panels A-E)",
        "",
        "Generated files:",
    ]
    for f in files:
        lines.append(f"- {f}")
    lines += [
        "",
        "Notes:",
        "- Panel A: donor-level somatic scores with row-wise z-score matrix (ordered by PC1_v2).",
        "- Panels C/D: Spearman statistics are recomputed directly from donor-level table.",
        "- Panel E: correlation matrix stars use Benjamini-Hochberg FDR (q-value).",
        "",
        "Build script:",
        "- ../supfig8_code/supfig8_source_data.py",
        "- ../supfig8_code/supfig8_plot_panels.py",
    ]
    (OUT_DIR / "README.txt").write_text("\n".join(lines), encoding="utf-8")


def write_panel_confidence(panel_map: dict[str, dict]) -> None:
    rows = []
    note_map = {
        "A": "Panel A 的 donor-level 体细胞niche score原始表与对应热图导出文件均已定位。",
        "B": "Panel B 的 PCA 坐标（PC1_v2/PC2_v2）可直接回溯到 donor-level 原始表。",
        "C": "Panel C 的 Spearman rho/p 可由 donor-level 原始表精确重算（与图注数值一致）。",
        "D": "Panel D 的 Spearman rho/p 可由 donor-level 原始表精确重算（与图注数值一致）。",
        "E": "Panel E 的相关矩阵与显著性星号可由 donor-level 原始表按 BH-FDR 精确重算。",
    }
    evidence_map = {
        "A": "Supplementary Figure 8.pdf; source data/supfig8_source_data/SupFig8A_somatic_scores_heatmap_z_matrix.csv",
        "B": "Supplementary Figure 8.pdf; source data/supfig8_source_data/SupFig8B_pca_points.csv",
        "C": "Supplementary Figure 8.pdf; source data/supfig8_source_data/SupFig8C_spearman_summary.csv",
        "D": "Supplementary Figure 8.pdf; source data/supfig8_source_data/SupFig8D_spearman_summary.csv",
        "E": "Supplementary Figure 8.pdf; source data/supfig8_source_data/SupFig8E_score_correlation_text_matrix.csv",
    }

    for panel in ["A", "B", "C", "D", "E"]:
        info = panel_map[panel]
        source_paths = [f"source data/supfig8_source_data/{x}" for x in info["generated"]]
        raw_paths = [f"source data/supfig8_source_data/{x}" for x in info["raw"]]
        source_data = "; ".join(source_paths + raw_paths)
        code = "; ".join(
            [f"source data/supfig8_code/{Path(x).name}" for x in info["code"] if x.startswith("../supfig8_code/")]
        )
        rows.append(
            {
                "panel": panel,
                "source_data": source_data,
                "code": code,
                "confidence": "High" if info["status"] == "available" else "Medium",
                "evidence": evidence_map[panel],
                "notes": note_map[panel],
            }
        )

    csv_path = BASE / "SupFig8_panel_source_code_confidence.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)

    md_lines = [
        "# Supplementary Figure 8 panel-source_code-confidence checklist",
        "",
        "## Summary",
        "- Main layout file: `Supplementary Figure 8.ai`",
        "- Latest exported figure: `Supplementary Figure 8.png`",
        "- Unified source-data workspace: `source data/supfig8_source_data`",
        "",
        "## Panel mapping",
        "",
        "| Panel | Source data | Code | Confidence | Evidence | Notes |",
        "|---|---|---|---|---|---|",
    ]
    for row in rows:
        md_lines.append(
            f"| {row['panel']} | {row['source_data']} | {row['code']} | {row['confidence']} | {row['evidence']} | {row['notes']} |"
        )
    md_lines += [
        "",
        "## Confidence rule used",
        "- High: has direct donor-level raw numeric table and can be recomputed exactly.",
        "- Medium: has strong source evidence, but missing key numeric raw input.",
    ]
    (BASE / "SupFig8_panel_source_code_confidence.md").write_text("\n".join(md_lines), encoding="utf-8")


def main() -> None:
    ensure_dirs()
    shared_raw = add_main_and_shared_raw()
    df = prepare_base_df()

    panel_map = {
        "A": build_panel_a(df, shared_raw),
        "B": build_panel_b(df, shared_raw),
        "C": build_panel_c(df, shared_raw),
        "D": build_panel_d(df, shared_raw),
        "E": build_panel_e(df, shared_raw),
    }

    write_file_mapping(panel_map)
    write_availability(panel_map)
    write_readme()
    write_panel_confidence(panel_map)


if __name__ == "__main__":
    main()
