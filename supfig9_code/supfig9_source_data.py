#!/usr/bin/env python3
from __future__ import annotations

import gzip
import math
import re
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

BASE = Path(__file__).resolve().parents[1]  # .../source data
ROOT = BASE.parent
OUT_DIR = BASE / "supfig9_source_data"
RAW_SUPPORT_DIR = OUT_DIR / "raw_support"

SUPFIG9_AI = ROOT / "Supplementary Figure 9.ai"
SUPFIG9_PDF = ROOT / "Supplementary Figure 9.pdf"
SUPFIG9_PNG = ROOT / "Supplementary Figure 9.png"

EXT_DIR = ROOT / "原来作图的数据" / "External validition"
DIV_DIR = EXT_DIR / "ExternalValidation_divergencePlots_CONTINUOUS_bundle"
GENE_DIR = EXT_DIR / "ExternalValidation_GeneScores8_cellPvalues_bundle"
HK_DIR = EXT_DIR / "ExternalValidation_GSE149512_metabolic_uncoupling_pvaluesCD_bundle_v4"
RAW149_DIR = EXT_DIR / "GSE149512_RAW"
SIG_DIR = EXT_DIR / "各种基因集"

PANEL_A_CELL = DIV_DIR / "GSE235321_externalValidation_cellLevel_pseudotime_entropy_moduleScores.csv"
PANEL_A_SMOOTH = DIV_DIR / "GSE235321_developmental_divergence_entropy_smoothedCurve_FINAL_CONTINUOUS.csv"
PANEL_A_SCRIPT = DIV_DIR / "external_validation_divergence_plot_unified_CONTINUOUS.py"
PANEL_A_README = DIV_DIR / "README_divergence_CONTINUOUS.txt"
PANEL_A_PDF = DIV_DIR / "ExternalValidation_GSE235321_developmental_divergence_entropy_FINAL_CONTINUOUS.pdf"
PANEL_A_PNG = DIV_DIR / "ExternalValidation_GSE235321_developmental_divergence_entropy_FINAL_CONTINUOUS.png"

GSE149_CELL = GENE_DIR / "GSE149512_cellLevel_geneScores8_withCellPvalues_data.csv"
GSE149_PVAL = GENE_DIR / "GSE149512_geneScores8_cellLevel_pvalues.csv"
GSE149_SCRIPT = GENE_DIR / "external_validation_geneScores8_cellPvalues.py"
GSE149_README = GENE_DIR / "README_geneScores8_cellPvalues.txt"
GSE149_SIG_MANIFEST = GENE_DIR / "signature_manifest_externalValidation.csv"
GSE149_PDF = GENE_DIR / "ExternalValidation_GSE149512_GeneScores_8panels_cellPvalues.pdf"
GSE149_PNG = GENE_DIR / "ExternalValidation_GSE149512_GeneScores_8panels_cellPvalues.png"

HK_DONOR_CSV = HK_DIR / "GSE149512_donorLevel_export_uptake_uncoupling_spermatidFraction.csv"
HK_DONOR_XLSX = HK_DIR / "GSE149512_donorLevel_export_uptake_uncoupling_spermatidFraction.xlsx"
HK_PVALUE_CSV = HK_DIR / "GSE149512_metabolic_uncoupling_panelsCD_pvalues.csv"
HK_SCRIPT = HK_DIR / "external_validation_metabolic_uncoupling_multiplot_v4_pvaluesCD.py"
HK_MULTIPANEL_PDF = HK_DIR / "ExternalValidation_GSE149512_metabolic_uncoupling_MULTIPANEL_v4_pvaluesCD.pdf"
HK_MULTIPANEL_PNG = HK_DIR / "ExternalValidation_GSE149512_metabolic_uncoupling_MULTIPANEL_v4_pvaluesCD.png"
HK_RESTAT_PDF = HK_DIR / "ExternalValidation_GSE149512_export_vs_uptake_restat.pdf"

SIG_ST_EXPORT = SIG_DIR / "ST_Lactate_Production_Export_Score.csv"
SIG_GERM_UPTAKE = SIG_DIR / "Germ_Lactate_Uptake_Oxidation_Score.csv"

GROUP_ORDER_149 = ["Control", "AZFa", "iNOA", "KS"]
GROUP_ORDER_HJK = ["Normal", "AZFa", "iNOA", "KS"]

PANEL_BG = [
    ("B", "ST maturity", "ST_maturity", lambda d: d["celltype"] == "ST"),
    ("C", "LC immaturity", "LCimmature_z", lambda d: d["celltype"] == "LC"),
    ("D", "Somatic cytokine signaling", "Cytokine_z", lambda d: d["celltype"] != "Germ"),
    ("E", "Somatic SASP", "SASP", lambda d: d["celltype"] != "Germ"),
    ("F", "ST BTB integrity", "BTB", lambda d: d["celltype"] == "ST"),
    ("G", "Germ apoptosis", "GermApoptosis_z", lambda d: d["celltype"] == "Germ"),
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


def star(p: float | None) -> str:
    if p is None or pd.isna(p):
        return ""
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def zscore(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    mu = x.mean(skipna=True)
    sd = x.std(skipna=True, ddof=0)
    if pd.isna(sd) or float(sd) == 0.0:
        return pd.Series(np.zeros(len(x)), index=s.index)
    return (x - mu) / sd


def compute_two_signatures_from_matrix(path: Path, cells_keep: set[str], st_genes: set[str], uptake_genes: set[str]) -> pd.DataFrame:
    union = st_genes | uptake_genes
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n")
        cols = [c.strip('"') for c in header.split(",")]
        cell_cols = cols[1:]

        idx_map = {c: i for i, c in enumerate(cell_cols)}
        keep_cells = [c for c in cell_cols if c in cells_keep]
        keep_idx = np.array([idx_map[c] for c in keep_cells], dtype=int)

        st_sum = np.zeros(len(keep_cells), dtype=np.float64)
        up_sum = np.zeros(len(keep_cells), dtype=np.float64)
        st_n = 0
        up_n = 0

        for line in f:
            if not line:
                continue
            try:
                gene_raw, rest = line.rstrip("\n").split(",", 1)
            except ValueError:
                continue
            gene = gene_raw.strip('"')
            if gene not in union:
                continue

            vals = np.fromstring(rest, sep=",", dtype=np.float64)
            vals = vals[keep_idx]
            logvals = np.log1p(vals)

            if gene in st_genes:
                st_sum += logvals
                st_n += 1
            if gene in uptake_genes:
                up_sum += logvals
                up_n += 1

    out = pd.DataFrame(
        {
            "cell_id": keep_cells,
            "st_export_raw": st_sum / st_n if st_n > 0 else np.nan,
            "germ_uptake_raw": up_sum / up_n if up_n > 0 else np.nan,
        }
    )
    return out


def add_shared_raw() -> list[str]:
    raw = []
    for src, alias in [
        (SUPFIG9_AI, "SupFig9_main.ai"),
        (SUPFIG9_PDF, "SupFig9_main.pdf"),
        (SUPFIG9_PNG, "SupFig9_main.png"),
    ]:
        rel = copy_to_raw_support(src, alias)
        if rel:
            raw.append(rel)
    return raw


def build_panel_a(shared_raw: list[str]) -> dict:
    generated: list[str] = []
    raw = list(shared_raw)
    for src in [PANEL_A_CELL, PANEL_A_SMOOTH, PANEL_A_SCRIPT, PANEL_A_README, PANEL_A_PDF, PANEL_A_PNG]:
        rel = copy_to_raw_support(src, f"SupFig9A_{src.name}")
        if rel:
            raw.append(rel)

    if not PANEL_A_CELL.exists():
        pd.DataFrame([{"panel": "A", "status": "missing", "note": "GSE235321 pseudotime entropy file missing"}]).to_csv(
            OUT_DIR / "SupFig9A_source_note.csv", index=False
        )
        generated.append("SupFig9A_source_note.csv")
        return {
            "subpanel": "Developmental divergence (GSE235321)",
            "generated": generated,
            "raw": raw,
            "code": ["../supfig9_code/supfig9_source_data.py", "../supfig9_code/supfig9_plot_panels.py"],
            "status": "missing",
            "availability_note": "GSE235321 pseudotime entropy source file not found",
        }

    df = pd.read_csv(PANEL_A_CELL)
    df = df.rename(columns={"Cell": "cell_id", "Group": "group", "Pseudotime2": "pseudotime", "MixEntropy": "mix_entropy", "PlotStage": "plot_stage"})
    keep = [c for c in ["cell_id", "group", "Disease", "Stage", "plot_stage", "pseudotime", "mix_entropy", "Batch", "OXPHOS", "Apoptosis"] if c in df.columns]
    df = df[keep].copy()
    df.to_csv(OUT_DIR / "SupFig9A_celllevel_pseudotime_entropy.csv", index=False)
    generated.append("SupFig9A_celllevel_pseudotime_entropy.csv")

    if PANEL_A_SMOOTH.exists():
        sm = pd.read_csv(PANEL_A_SMOOTH)
        # normalize column names for downstream plotting
        colmap = {c: c for c in sm.columns}
        if "pseudotime" not in sm.columns:
            for c in sm.columns:
                if "pseudo" in c.lower():
                    colmap[c] = "pseudotime"
                    break
        if "entropy_smooth" not in sm.columns:
            for c in sm.columns:
                if "entropy" in c.lower():
                    colmap[c] = "entropy_smooth"
                    break
        sm = sm.rename(columns=colmap)
        sm.to_csv(OUT_DIR / "SupFig9A_smoothed_curve.csv", index=False)
        generated.append("SupFig9A_smoothed_curve.csv")
    else:
        sm = pd.DataFrame(columns=["pseudotime", "entropy_smooth"])

    # Stage ranges and arrest point inside Late SPC
    ranges = []
    for stage in ["SPG (Stem)", "Early SPC", "Late SPC", "Round Sperm", "Elongated Sperm"]:
        sub = df[df["plot_stage"] == stage]
        if sub.empty:
            continue
        ranges.append({"plot_stage": stage, "x_min": float(sub["pseudotime"].min()), "x_max": float(sub["pseudotime"].max())})
    stage_df = pd.DataFrame(ranges)
    stage_df.to_csv(OUT_DIR / "SupFig9A_stage_ranges.csv", index=False)
    generated.append("SupFig9A_stage_ranges.csv")

    arrest = np.nan
    if not sm.empty and {"pseudotime", "entropy_smooth"}.issubset(sm.columns) and not stage_df.empty:
        late = stage_df[stage_df["plot_stage"] == "Late SPC"]
        if not late.empty:
            x0, x1 = float(late["x_min"].iloc[0]), float(late["x_max"].iloc[0])
            m = (sm["pseudotime"] >= x0) & (sm["pseudotime"] <= x1)
            if m.any():
                smm = sm.loc[m].sort_values("entropy_smooth")
                if not smm.empty:
                    arrest = float(smm.iloc[0]["pseudotime"])
    pd.DataFrame([{"panel": "A", "arrest_point_pseudotime": arrest}]).to_csv(
        OUT_DIR / "SupFig9A_arrest_point.csv", index=False
    )
    generated.append("SupFig9A_arrest_point.csv")

    grp = df.groupby("group", observed=False).size().reset_index(name="n_cells")
    grp.to_csv(OUT_DIR / "SupFig9A_group_cell_counts.csv", index=False)
    generated.append("SupFig9A_group_cell_counts.csv")

    return {
        "subpanel": "Developmental divergence (GSE235321)",
        "generated": generated,
        "raw": raw,
        "code": ["../supfig9_code/supfig9_source_data.py", "../supfig9_code/supfig9_plot_panels.py"],
        "status": "available",
        "availability_note": "Direct cell-level pseudotime entropy table and smoothed curve are available.",
    }


def build_panels_b_to_g(shared_raw: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if not GSE149_CELL.exists() or not GSE149_PVAL.exists():
        for letter, title, _, _ in PANEL_BG:
            pd.DataFrame([{"panel": letter, "status": "missing", "note": "GSE149512 cell-level table or pvalue file missing"}]).to_csv(
                OUT_DIR / f"SupFig9{letter}_source_note.csv", index=False
            )
            out[letter] = {
                "subpanel": title,
                "generated": [f"SupFig9{letter}_source_note.csv"],
                "raw": list(shared_raw),
                "code": ["../supfig9_code/supfig9_source_data.py", "../supfig9_code/supfig9_plot_panels.py"],
                "status": "missing",
                "availability_note": "GSE149512 source table missing",
            }
        return out

    shared_bg_raw = list(shared_raw)
    for src in [GSE149_CELL, GSE149_PVAL, GSE149_SCRIPT, GSE149_README, GSE149_SIG_MANIFEST, GSE149_PDF, GSE149_PNG]:
        rel = copy_to_raw_support(src, f"SupFig9BG_{src.name}")
        if rel:
            shared_bg_raw.append(rel)

    df = pd.read_csv(GSE149_CELL)
    pval = pd.read_csv(GSE149_PVAL)

    # Global long table for B-G
    long_rows = []
    for letter, title, metric, mask_fn in PANEL_BG:
        sub = df[mask_fn(df)].copy()
        keep = ["cell_id", "sample_code", "group", "celltype", "germ_stage", metric]
        sub = sub[keep].rename(columns={metric: "value"})
        sub["panel"] = letter
        sub["panel_title"] = title
        sub["metric"] = metric
        long_rows.append(sub)
    long_df = pd.concat(long_rows, ignore_index=True)
    long_df.to_csv(OUT_DIR / "SupFig9BG_gene_score_cells_long.csv", index=False)

    pval_all = pval.copy()
    pval_all["star"] = pval_all["p_value"].apply(star)
    pval_all.to_csv(OUT_DIR / "SupFig9BG_cell_level_pvalues_full.csv", index=False)

    for letter, title, metric, mask_fn in PANEL_BG:
        generated = [
            "SupFig9BG_gene_score_cells_long.csv",
            "SupFig9BG_cell_level_pvalues_full.csv",
        ]

        sub = df[mask_fn(df)].copy()
        keep = ["cell_id", "sample_code", "group", "celltype", "germ_stage", metric]
        sub = sub[keep].rename(columns={metric: "value"})
        sub["group"] = pd.Categorical(sub["group"], categories=GROUP_ORDER_149, ordered=True)
        sub = sub.sort_values(["group", "sample_code", "cell_id"]).reset_index(drop=True)
        panel_cells = f"SupFig9{letter}_{metric}_cells.csv"
        sub.to_csv(OUT_DIR / panel_cells, index=False)
        generated.append(panel_cells)

        donor = (
            sub.groupby(["sample_code", "group"], observed=False)["value"]
            .agg(n_cells="count", mean="mean", median="median", q1=lambda x: x.quantile(0.25), q3=lambda x: x.quantile(0.75), min="min", max="max")
            .reset_index()
            .sort_values(["group", "sample_code"])
            .reset_index(drop=True)
        )
        panel_donor = f"SupFig9{letter}_{metric}_donor_summary.csv"
        donor.to_csv(OUT_DIR / panel_donor, index=False)
        generated.append(panel_donor)

        pv = pval_all[pval_all["metric"] == metric].copy()
        pv["panel"] = letter
        pv = pv[["panel", "dataset", "metric", "comparison", "n_control", "n_group", "p_value", "star"]]
        panel_pv = f"SupFig9{letter}_{metric}_pvalues.csv"
        pv.to_csv(OUT_DIR / panel_pv, index=False)
        generated.append(panel_pv)

        out[letter] = {
            "subpanel": title,
            "generated": generated,
            "raw": shared_bg_raw,
            "code": ["../supfig9_code/supfig9_source_data.py", "../supfig9_code/supfig9_plot_panels.py"],
            "status": "available",
            "availability_note": "Direct GSE149512 cell-level scores and cell-level MWU p-values are available.",
        }

    return out


def build_panels_h_to_k(shared_raw: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    raw = list(shared_raw)
    for src in [GSE149_CELL, SIG_ST_EXPORT, SIG_GERM_UPTAKE]:
        rel = copy_to_raw_support(src, f"SupFig9HK_{src.name}")
        if rel:
            raw.append(rel)

    # Prefer direct donor-level H-K bundle if available.
    for src in [HK_DONOR_CSV, HK_DONOR_XLSX, HK_PVALUE_CSV, HK_SCRIPT, HK_MULTIPANEL_PDF, HK_MULTIPANEL_PNG, HK_RESTAT_PDF]:
        rel = copy_to_raw_support(src, f"SupFig9HK_{src.name}")
        if rel:
            raw.append(rel)

    if HK_DONOR_CSV.exists():
        donor_raw = pd.read_csv(HK_DONOR_CSV)
        required = {"sample_code", "Class", "export_ST_z", "uptake_germ_z", "uncoupling_z", "spermatid_fraction"}
        missing = sorted(required - set(donor_raw.columns))
        if not missing:
            generated_hk: list[str] = []

            # Keep a direct archival copy from the original bundle.
            donor_raw.to_csv(OUT_DIR / "SupFig9HK_direct_donor_table.csv", index=False)
            generated_hk.append("SupFig9HK_direct_donor_table.csv")

            if HK_PVALUE_CSV.exists():
                pd.read_csv(HK_PVALUE_CSV).to_csv(OUT_DIR / "SupFig9HK_group_mwu_pvalues.csv", index=False)
                generated_hk.append("SupFig9HK_group_mwu_pvalues.csv")

            donor = donor_raw.copy()
            donor["group_plot"] = donor["Class"].astype(str)
            donor["group"] = donor["group"] if "group" in donor.columns else donor["Class"]

            donor["st_export_z"] = pd.to_numeric(donor["export_ST_z"], errors="coerce")
            donor["germ_uptake_z"] = pd.to_numeric(donor["uptake_germ_z"], errors="coerce")
            donor["uncoupling_index"] = pd.to_numeric(donor["uncoupling_z"], errors="coerce")
            donor["postmeiotic_fraction"] = pd.to_numeric(donor["spermatid_fraction"], errors="coerce")

            if "n_germ" in donor.columns:
                donor["n_germ"] = pd.to_numeric(donor["n_germ"], errors="coerce")
            elif "n_Germ" in donor.columns:
                donor["n_germ"] = pd.to_numeric(donor["n_Germ"], errors="coerce")
            else:
                donor["n_germ"] = np.nan

            if "n_spermatid" in donor.columns:
                donor["n_postmeiotic"] = pd.to_numeric(donor["n_spermatid"], errors="coerce")
            elif "n_postmeiotic" in donor.columns:
                donor["n_postmeiotic"] = pd.to_numeric(donor["n_postmeiotic"], errors="coerce")
            else:
                donor["n_postmeiotic"] = np.nan

            donor["group_plot"] = pd.Categorical(donor["group_plot"], categories=GROUP_ORDER_HJK, ordered=True)
            donor = donor.sort_values(["group_plot", "sample_code"]).reset_index(drop=True)
            donor.to_csv(OUT_DIR / "SupFig9HK_donor_metrics.csv", index=False)
            generated_hk.append("SupFig9HK_donor_metrics.csv")

            h = donor[["sample_code", "group_plot", "st_export_z", "germ_uptake_z"]].copy()
            h.to_csv(OUT_DIR / "SupFig9H_coupling_scatter_donor_points.csv", index=False)
            generated_hk.append("SupFig9H_coupling_scatter_donor_points.csv")

            i_df = donor[["sample_code", "group_plot", "st_export_z"]].copy()
            i_df.to_csv(OUT_DIR / "SupFig9I_st_lactate_export_by_group.csv", index=False)
            generated_hk.append("SupFig9I_st_lactate_export_by_group.csv")

            j = donor[["sample_code", "group_plot", "uncoupling_index", "postmeiotic_fraction", "n_germ", "n_postmeiotic"]].copy()
            j.to_csv(OUT_DIR / "SupFig9J_uncoupling_vs_postmeiotic.csv", index=False)
            generated_hk.append("SupFig9J_uncoupling_vs_postmeiotic.csv")

            k_df = donor[["sample_code", "group_plot", "germ_uptake_z"]].copy()
            k_df.to_csv(OUT_DIR / "SupFig9K_germ_uptake_by_group.csv", index=False)
            generated_hk.append("SupFig9K_germ_uptake_by_group.csv")

            rho_h, p_h = spearmanr(h["st_export_z"], h["germ_uptake_z"], nan_policy="omit")
            rho_j, p_j = spearmanr(j["uncoupling_index"], j["postmeiotic_fraction"], nan_policy="omit")

            normal_unc = j.loc[j["group_plot"] == "Normal", "uncoupling_index"].dropna().to_numpy()
            azo_unc = j.loc[j["group_plot"] != "Normal", "uncoupling_index"].dropna().to_numpy()
            p_unc_mwu = mannwhitneyu(normal_unc, azo_unc, alternative="two-sided").pvalue if (len(normal_unc) > 0 and len(azo_unc) > 0) else np.nan

            rep_idx = np.repeat(np.arange(len(j)), pd.to_numeric(j["n_germ"], errors="coerce").fillna(0).astype(int).values)
            if len(rep_idx) > 1:
                rho_j_w, p_j_w = spearmanr(j.iloc[rep_idx]["uncoupling_index"], j.iloc[rep_idx]["postmeiotic_fraction"], nan_policy="omit")
            else:
                rho_j_w, p_j_w = np.nan, np.nan

            pd.DataFrame(
                [
                    {"panel": "H", "mode": "donor_level", "n": int(h.shape[0]), "spearman_rho": rho_h, "p_value": p_h},
                    {"panel": "H", "mode": "figure_annotation", "n": int(h.shape[0]), "spearman_rho": rho_h, "p_value": p_unc_mwu},
                    {"panel": "J", "mode": "donor_level", "n": int(j.shape[0]), "spearman_rho": rho_j, "p_value": p_j},
                    {"panel": "J", "mode": "figure_annotation", "n": int(j.shape[0]), "spearman_rho": rho_j, "p_value": p_unc_mwu},
                    {
                        "panel": "J",
                        "mode": "weighted_by_n_germ",
                        "n": int(len(rep_idx)),
                        "spearman_rho": rho_j_w,
                        "p_value": p_j_w,
                    },
                    {"panel": "HK", "mode": "uncoupling_normal_vs_azo_mwu", "n": int(len(normal_unc) + len(azo_unc)), "spearman_rho": np.nan, "p_value": p_unc_mwu},
                ]
            ).to_csv(OUT_DIR / "SupFig9HJ_correlation_summary.csv", index=False)
            generated_hk.append("SupFig9HJ_correlation_summary.csv")

            pd.DataFrame(
                [
                    {
                        "panel": "H-K",
                        "note": "Direct donor-level external-validation H-K table is available and used from ExternalValidation_GSE149512_metabolic_uncoupling_pvaluesCD_bundle_v4.",
                    }
                ]
            ).to_csv(OUT_DIR / "SupFig9HK_source_note.csv", index=False)
            generated_hk.append("SupFig9HK_source_note.csv")

            out["H"] = {
                "subpanel": "Coupling between lactate export and uptake",
                "generated": generated_hk,
                "raw": raw,
                "code": ["../supfig9_code/supfig9_source_data.py", "../supfig9_code/supfig9_plot_panels.py"],
                "status": "available",
                "availability_note": "H donor-level table is directly available in external-validation metabolic uncoupling bundle.",
            }
            out["I"] = {
                "subpanel": "ST lactate export by group",
                "generated": generated_hk,
                "raw": raw,
                "code": ["../supfig9_code/supfig9_source_data.py", "../supfig9_code/supfig9_plot_panels.py"],
                "status": "available",
                "availability_note": "I donor-level table is directly available in external-validation metabolic uncoupling bundle.",
            }
            out["J"] = {
                "subpanel": "Uncoupling predicts post-meiotic output",
                "generated": generated_hk,
                "raw": raw,
                "code": ["../supfig9_code/supfig9_source_data.py", "../supfig9_code/supfig9_plot_panels.py"],
                "status": "available",
                "availability_note": "J donor-level table is directly available; figure annotation p-value uses Normal vs Azoospermia MWU from the same table.",
            }
            out["K"] = {
                "subpanel": "Germ lactate uptake/oxidation by group",
                "generated": generated_hk,
                "raw": raw,
                "code": ["../supfig9_code/supfig9_source_data.py", "../supfig9_code/supfig9_plot_panels.py"],
                "status": "available",
                "availability_note": "K donor-level table is directly available in external-validation metabolic uncoupling bundle.",
            }
            return out

    if not GSE149_CELL.exists() or not RAW149_DIR.exists():
        for letter, title in [("H", "Coupling between lactate export and uptake"), ("I", "ST lactate export by group"), ("J", "Uncoupling predicts post-meiotic output"), ("K", "Germ lactate uptake/oxidation by group")]:
            pd.DataFrame([{"panel": letter, "status": "missing", "note": "GSE149512 table or raw matrices missing"}]).to_csv(
                OUT_DIR / f"SupFig9{letter}_source_note.csv", index=False
            )
            out[letter] = {
                "subpanel": title,
                "generated": [f"SupFig9{letter}_source_note.csv"],
                "raw": raw,
                "code": ["../supfig9_code/supfig9_source_data.py", "../supfig9_code/supfig9_plot_panels.py"],
                "status": "missing",
                "availability_note": "Raw matrices for lactate signatures are missing",
            }
        return out

    df = pd.read_csv(GSE149_CELL, usecols=["cell_id", "sample_code", "group", "celltype", "germ_stage"])
    df["group_plot"] = df["group"].replace({"Control": "Normal"})

    st_genes = set(pd.read_csv(SIG_ST_EXPORT).iloc[:, 0].astype(str).str.strip())
    uptake_genes = set(pd.read_csv(SIG_GERM_UPTAKE).iloc[:, 0].astype(str).str.strip())

    sig_rows = []
    for g in sorted(st_genes):
        sig_rows.append({"signature": "ST_Lactate_Production_Export_Score", "gene": g})
    for g in sorted(uptake_genes):
        sig_rows.append({"signature": "Germ_Lactate_Uptake_Oxidation_Score", "gene": g})
    pd.DataFrame(sig_rows).to_csv(OUT_DIR / "SupFig9HK_signature_gene_sets_long.csv", index=False)

    sample_to_file = {}
    for p in RAW149_DIR.glob("*_LZ*matrix.csv.gz"):
        m = re.search(r"_(LZ\d+)matrix\.csv\.gz$", p.name)
        if m:
            sample_to_file[m.group(1)] = p

    registry_rows = []
    for sample in sorted(df["sample_code"].unique()):
        p = sample_to_file.get(sample)
        registry_rows.append(
            {
                "sample_code": sample,
                "matrix_file": str(p) if p else "",
                "exists": bool(p and p.exists()),
            }
        )
    reg_df = pd.DataFrame(registry_rows)
    reg_df.to_csv(OUT_DIR / "SupFig9HK_raw_matrix_registry.csv", index=False)

    parts = []
    time_rows = []
    for sample in sorted(df["sample_code"].unique()):
        matrix = sample_to_file.get(sample)
        if matrix is None:
            continue
        cells_keep = set(df.loc[df["sample_code"] == sample, "cell_id"])
        t0 = time.time()
        sc = compute_two_signatures_from_matrix(matrix, cells_keep, st_genes, uptake_genes)
        elapsed = time.time() - t0
        sc["sample_code"] = sample
        parts.append(sc)
        time_rows.append({"sample_code": sample, "n_cells_scored": int(sc.shape[0]), "seconds": round(elapsed, 3)})

    if not parts:
        # Create placeholder notes
        pd.DataFrame([{"panel": "H-K", "status": "missing", "note": "No sample matrix could be parsed."}]).to_csv(
            OUT_DIR / "SupFig9HK_source_note.csv", index=False
        )
        for letter, title in [("H", "Coupling between lactate export and uptake"), ("I", "ST lactate export by group"), ("J", "Uncoupling predicts post-meiotic output"), ("K", "Germ lactate uptake/oxidation by group")]:
            out[letter] = {
                "subpanel": title,
                "generated": ["SupFig9HK_source_note.csv", "SupFig9HK_raw_matrix_registry.csv"],
                "raw": raw,
                "code": ["../supfig9_code/supfig9_source_data.py", "../supfig9_code/supfig9_plot_panels.py"],
                "status": "missing",
                "availability_note": "Failed to compute per-cell lactate signatures from raw matrices",
            }
        return out

    score_df = pd.concat(parts, ignore_index=True)
    score_df.to_csv(OUT_DIR / "SupFig9HK_cell_scores_raw.csv", index=False)
    pd.DataFrame(time_rows).to_csv(OUT_DIR / "SupFig9HK_matrix_scoring_runtime.csv", index=False)

    merged = df.merge(score_df, on=["cell_id", "sample_code"], how="left")
    merged.to_csv(OUT_DIR / "SupFig9HK_cell_scores_with_metadata.csv", index=False)

    st_donor = (
        merged[merged["celltype"] == "ST"]
        .groupby(["sample_code", "group", "group_plot"], observed=False)["st_export_raw"]
        .median()
        .reset_index(name="st_export_raw")
    )
    germ_donor = (
        merged[merged["celltype"] == "Germ"]
        .groupby(["sample_code", "group", "group_plot"], observed=False)["germ_uptake_raw"]
        .median()
        .reset_index(name="germ_uptake_raw")
    )
    germ_n = (
        merged[merged["celltype"] == "Germ"]
        .groupby(["sample_code", "group", "group_plot"], observed=False)
        .size()
        .reset_index(name="n_germ")
    )
    post_n = (
        merged[(merged["celltype"] == "Germ") & (merged["germ_stage"].isin(["Round_Sperm", "Elongated_Sperm"]))]
        .groupby(["sample_code", "group", "group_plot"], observed=False)
        .size()
        .reset_index(name="n_postmeiotic")
    )

    donor = (
        st_donor.merge(germ_donor, on=["sample_code", "group", "group_plot"], how="outer")
        .merge(germ_n, on=["sample_code", "group", "group_plot"], how="outer")
        .merge(post_n, on=["sample_code", "group", "group_plot"], how="left")
        .fillna({"n_postmeiotic": 0})
    )

    donor["postmeiotic_fraction"] = donor["n_postmeiotic"] / donor["n_germ"]
    donor["st_export_z"] = zscore(donor["st_export_raw"])
    donor["germ_uptake_z"] = zscore(donor["germ_uptake_raw"])
    donor["uncoupling_index"] = donor["st_export_z"] - donor["germ_uptake_z"]

    donor["group_plot"] = pd.Categorical(donor["group_plot"], categories=GROUP_ORDER_HJK, ordered=True)
    donor = donor.sort_values(["group_plot", "sample_code"]).reset_index(drop=True)
    donor.to_csv(OUT_DIR / "SupFig9HK_donor_metrics.csv", index=False)

    h = donor[["sample_code", "group_plot", "st_export_z", "germ_uptake_z"]].copy()
    h.to_csv(OUT_DIR / "SupFig9H_coupling_scatter_donor_points.csv", index=False)

    i_df = donor[["sample_code", "group_plot", "st_export_z"]].copy()
    i_df.to_csv(OUT_DIR / "SupFig9I_st_lactate_export_by_group.csv", index=False)

    j = donor[["sample_code", "group_plot", "uncoupling_index", "postmeiotic_fraction", "n_germ", "n_postmeiotic"]].copy()
    j.to_csv(OUT_DIR / "SupFig9J_uncoupling_vs_postmeiotic.csv", index=False)

    k_df = donor[["sample_code", "group_plot", "germ_uptake_z"]].copy()
    k_df.to_csv(OUT_DIR / "SupFig9K_germ_uptake_by_group.csv", index=False)

    # Correlation summaries
    rho_h, p_h = spearmanr(h["st_export_z"], h["germ_uptake_z"], nan_policy="omit")
    rho_j, p_j = spearmanr(j["uncoupling_index"], j["postmeiotic_fraction"], nan_policy="omit")

    rep_idx = np.repeat(np.arange(len(j)), j["n_germ"].fillna(0).astype(int).values)
    if len(rep_idx) > 1:
        rho_j_w, p_j_w = spearmanr(j.iloc[rep_idx]["uncoupling_index"], j.iloc[rep_idx]["postmeiotic_fraction"], nan_policy="omit")
    else:
        rho_j_w, p_j_w = np.nan, np.nan

    pd.DataFrame(
        [
            {"panel": "H", "mode": "donor_level", "n": int(h.shape[0]), "spearman_rho": rho_h, "p_value": p_h},
            {"panel": "J", "mode": "donor_level", "n": int(j.shape[0]), "spearman_rho": rho_j, "p_value": p_j},
            {
                "panel": "J",
                "mode": "weighted_by_n_germ",
                "n": int(len(rep_idx)),
                "spearman_rho": rho_j_w,
                "p_value": p_j_w,
            },
        ]
    ).to_csv(OUT_DIR / "SupFig9HJ_correlation_summary.csv", index=False)

    pd.DataFrame(
        [
            {
                "panel": "H-K",
                "note": "Direct precomputed donor-level H-K table was not found in workspace. Donor metrics were reconstructed from GSE149512 raw matrices using provided signature gene lists.",
            }
        ]
    ).to_csv(OUT_DIR / "SupFig9HK_source_note.csv", index=False)

    generated_hk = [
        "SupFig9HK_signature_gene_sets_long.csv",
        "SupFig9HK_raw_matrix_registry.csv",
        "SupFig9HK_cell_scores_raw.csv",
        "SupFig9HK_matrix_scoring_runtime.csv",
        "SupFig9HK_cell_scores_with_metadata.csv",
        "SupFig9HK_donor_metrics.csv",
        "SupFig9H_coupling_scatter_donor_points.csv",
        "SupFig9I_st_lactate_export_by_group.csv",
        "SupFig9J_uncoupling_vs_postmeiotic.csv",
        "SupFig9K_germ_uptake_by_group.csv",
        "SupFig9HJ_correlation_summary.csv",
        "SupFig9HK_source_note.csv",
    ]

    out["H"] = {
        "subpanel": "Coupling between lactate export and uptake",
        "generated": generated_hk,
        "raw": raw,
        "code": ["../supfig9_code/supfig9_source_data.py", "../supfig9_code/supfig9_plot_panels.py"],
        "status": "partial",
        "availability_note": "H donor points are reconstructed from raw matrices + signature genes; no direct original donor summary file found.",
    }
    out["I"] = {
        "subpanel": "ST lactate export by group",
        "generated": generated_hk,
        "raw": raw,
        "code": ["../supfig9_code/supfig9_source_data.py", "../supfig9_code/supfig9_plot_panels.py"],
        "status": "partial",
        "availability_note": "I box data are reconstructed from raw matrices + signature genes; no direct original donor summary file found.",
    }
    out["J"] = {
        "subpanel": "Uncoupling predicts post-meiotic output",
        "generated": generated_hk,
        "raw": raw,
        "code": ["../supfig9_code/supfig9_source_data.py", "../supfig9_code/supfig9_plot_panels.py"],
        "status": "partial",
        "availability_note": "J correlation uses reconstructed donor metrics and post-meiotic fractions from germ-stage counts.",
    }
    out["K"] = {
        "subpanel": "Germ lactate uptake/oxidation by group",
        "generated": generated_hk,
        "raw": raw,
        "code": ["../supfig9_code/supfig9_source_data.py", "../supfig9_code/supfig9_plot_panels.py"],
        "status": "partial",
        "availability_note": "K box data are reconstructed from raw matrices + signature genes; no direct original donor summary file found.",
    }
    return out


def write_file_mapping(panel_map: dict[str, dict]) -> None:
    rows = []
    for panel in list("ABCDEFGHIJK"):
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
    pd.DataFrame(rows).to_csv(OUT_DIR / "SupFig9_file_mapping.csv", index=False)


def write_availability(panel_map: dict[str, dict]) -> None:
    rows = []
    for panel in list("ABCDEFGHIJK"):
        rows.append({"panel": panel, "status": panel_map[panel]["status"], "note": panel_map[panel]["availability_note"]})
    pd.DataFrame(rows).to_csv(OUT_DIR / "SupFig9_panels_data_availability.csv", index=False)


def write_readme(panel_map: dict[str, dict]) -> None:
    listed = set()
    for panel in list("ABCDEFGHIJK"):
        listed.update(panel_map[panel]["generated"])
    listed.update({"SupFig9_file_mapping.csv", "SupFig9_panels_data_availability.csv", "README.txt"})
    files = sorted(listed)
    hk_note = "- Panels H-K: reconstructed donor-level metrics from GSE149512 raw matrices + signature gene lists."
    if panel_map.get("H", {}).get("status") == "available":
        hk_note = "- Panels H-K: directly from ExternalValidation_GSE149512_metabolic_uncoupling_pvaluesCD_bundle_v4 donor-level table."
    lines = [
        "Supplementary Figure 9 source data (Panels A-K)",
        "",
        "Generated files:",
    ]
    for f in files:
        lines.append(f"- {f}")
    lines += [
        "",
        "Notes:",
        "- Panel A: directly from ExternalValidation_divergencePlots_CONTINUOUS bundle (GSE235321).",
        "- Panels B-G: directly from GSE149512 cell-level 8-gene-score table and MWU p-value table.",
        hk_note,
        "",
        "Build script:",
        "- ../supfig9_code/supfig9_source_data.py",
        "- ../supfig9_code/supfig9_plot_panels.py",
    ]
    (OUT_DIR / "README.txt").write_text("\n".join(lines), encoding="utf-8")


def write_confidence(panel_map: dict[str, dict]) -> None:
    evidence = {
        "A": "Supplementary Figure 9.pdf; source data/supfig9_source_data/SupFig9A_celllevel_pseudotime_entropy.csv",
        "B": "Supplementary Figure 9.pdf; source data/supfig9_source_data/SupFig9B_ST_maturity_cells.csv",
        "C": "Supplementary Figure 9.pdf; source data/supfig9_source_data/SupFig9C_LCimmature_z_cells.csv",
        "D": "Supplementary Figure 9.pdf; source data/supfig9_source_data/SupFig9D_Cytokine_z_cells.csv",
        "E": "Supplementary Figure 9.pdf; source data/supfig9_source_data/SupFig9E_SASP_cells.csv",
        "F": "Supplementary Figure 9.pdf; source data/supfig9_source_data/SupFig9F_BTB_cells.csv",
        "G": "Supplementary Figure 9.pdf; source data/supfig9_source_data/SupFig9G_GermApoptosis_z_cells.csv",
        "H": "Supplementary Figure 9.pdf; source data/supfig9_source_data/SupFig9H_coupling_scatter_donor_points.csv",
        "I": "Supplementary Figure 9.pdf; source data/supfig9_source_data/SupFig9I_st_lactate_export_by_group.csv",
        "J": "Supplementary Figure 9.pdf; source data/supfig9_source_data/SupFig9J_uncoupling_vs_postmeiotic.csv",
        "K": "Supplementary Figure 9.pdf; source data/supfig9_source_data/SupFig9K_germ_uptake_by_group.csv",
    }

    notes = {
        "A": "Panel A 原始输入和绘图脚本均在 divergence bundle 中，可直接回溯。",
        "B": "Panel B 直接来自 GSE149512 cell-level score 表与对应 p-value 表。",
        "C": "Panel C 直接来自 GSE149512 cell-level score 表与对应 p-value 表。",
        "D": "Panel D 直接来自 GSE149512 cell-level score 表与对应 p-value 表。",
        "E": "Panel E 直接来自 GSE149512 cell-level score 表与对应 p-value 表。",
        "F": "Panel F 直接来自 GSE149512 cell-level score 表与对应 p-value 表。",
        "G": "Panel G 直接来自 GSE149512 cell-level score 表与对应 p-value 表。",
        "H": "Panel H donor 点位由 raw matrix + signature 基因集重建；未发现原始 donor 汇总表。",
        "I": "Panel I donor 点位由 raw matrix + signature 基因集重建；未发现原始 donor 汇总表。",
        "J": "Panel J donor 指标与 post-meiotic fraction 由重建表计算；图中相关性P值计算口径可能不同。",
        "K": "Panel K donor 点位由 raw matrix + signature 基因集重建；未发现原始 donor 汇总表。",
    }
    if all(panel_map[p]["status"] == "available" for p in ["H", "I", "J", "K"]):
        notes["H"] = "Panel H 直接来自 donor-level metabolic uncoupling bundle（GSE149512），无需重建。"
        notes["I"] = "Panel I 直接来自 donor-level metabolic uncoupling bundle（GSE149512），无需重建。"
        notes["J"] = "Panel J 直接来自 donor-level metabolic uncoupling bundle（GSE149512）；图中标注 p 值使用同表的 Normal vs Azoospermia MWU。"
        notes["K"] = "Panel K 直接来自 donor-level metabolic uncoupling bundle（GSE149512），无需重建。"

    rows = []
    for panel in list("ABCDEFGHIJK"):
        info = panel_map[panel]
        generated = [f"source data/supfig9_source_data/{x}" for x in info["generated"]]
        raws = [f"source data/supfig9_source_data/{x}" for x in info["raw"]]
        rows.append(
            {
                "panel": panel,
                "source_data": "; ".join(generated + raws),
                "code": "source data/supfig9_code/supfig9_source_data.py; source data/supfig9_code/supfig9_plot_panels.py",
                "confidence": "High" if info["status"] == "available" else "Medium",
                "evidence": evidence[panel],
                "notes": notes[panel],
            }
        )

    pd.DataFrame(rows).to_csv(BASE / "SupFig9_panel_source_code_confidence.csv", index=False)

    md = [
        "# Supplementary Figure 9 panel-source_code-confidence checklist",
        "",
        "## Summary",
        "- Main layout file: `Supplementary Figure 9.ai`",
        "- Latest exported figure: `Supplementary Figure 9.png`",
        "- Unified source-data workspace: `source data/supfig9_source_data`",
        "",
        "## Panel mapping",
        "",
        "| Panel | Source data | Code | Confidence | Evidence | Notes |",
        "|---|---|---|---|---|---|",
    ]
    for r in rows:
        md.append(
            f"| {r['panel']} | {r['source_data']} | {r['code']} | {r['confidence']} | {r['evidence']} | {r['notes']} |"
        )
    md += [
        "",
        "## Confidence rule used",
        "- High: direct raw numeric source located and panel can be reproduced directly.",
        "- Medium: panel data reconstructed from upstream raw inputs due missing direct donor-level final table.",
    ]
    (BASE / "SupFig9_panel_source_code_confidence.md").write_text("\n".join(md), encoding="utf-8")


def main() -> None:
    ensure_dirs()
    shared_raw = add_shared_raw()

    panel_map: dict[str, dict] = {}
    panel_map["A"] = build_panel_a(shared_raw)
    panel_map.update(build_panels_b_to_g(shared_raw))
    panel_map.update(build_panels_h_to_k(shared_raw))

    write_file_mapping(panel_map)
    write_availability(panel_map)
    write_readme(panel_map)
    write_confidence(panel_map)


if __name__ == "__main__":
    main()
