#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
External validation: 8-panel gene score violin plots with *cell-level* p-values.

Outputs (in current working directory unless you change paths):
- ExternalValidation_GSE149512_GeneScores_8panels_cellPvalues.pdf/png
- ExternalValidation_GSE235321_GeneScores_8panels_cellPvalues.pdf/png
- GSE149512_cellLevel_geneScores8_withCellPvalues_data.csv
- GSE235321_cellLevel_geneScores8_withCellPvalues_data.csv
- GSE149512_geneScores8_cellLevel_pvalues.csv
- GSE235321_geneScores8_cellLevel_pvalues.csv

Notes
-----
1) P-values are computed using Mann–Whitney U test on *cells* (pseudo-replication risk).
   Donor median points are overlaid for visualization, but p-values are per-cell as requested.

2) Font: script requests Arial; if Arial is not installed, matplotlib will fall back
   to Liberation Sans / DejaVu Sans.
"""

from __future__ import annotations

import os
import gzip
import time
from typing import Dict, List, Optional, Set

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# -----------------------------
# Global plotting style
# -----------------------------
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42  # embed TTF
plt.rcParams["ps.fonttype"] = 42


# -----------------------------
# Helpers
# -----------------------------
def mannwhitney_p(control: np.ndarray, other: np.ndarray) -> float:
    c = np.asarray(control, dtype=float)
    o = np.asarray(other, dtype=float)
    c = c[~np.isnan(c)]
    o = o[~np.isnan(o)]
    if len(c) == 0 or len(o) == 0:
        return np.nan
    try:
        return mannwhitneyu(c, o, alternative="two-sided", method="auto").pvalue
    except TypeError:
        return mannwhitneyu(c, o, alternative="two-sided").pvalue


def format_p(p: float) -> str:
    if np.isnan(p):
        return "NA"
    if p == 0:
        return "P<1e-300"
    if p < 1e-3:
        return f"P={p:.1e}"
    return f"P={p:.3f}"


def star(p: float) -> str:
    if np.isnan(p):
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


def zscore_by_reference(df: pd.DataFrame, value_col: str, ref_mask: pd.Series):
    ref = df.loc[ref_mask, value_col].astype(float)
    mu = ref.mean()
    sd = ref.std(ddof=0)
    if sd == 0 or np.isnan(sd):
        return df[value_col] * np.nan, mu, sd
    return (df[value_col] - mu) / sd, mu, sd


def compute_raw_scores_from_matrix_fast(
    file_path: str,
    cells_keep: Optional[Set[str]],
    gene_to_sets: Dict[str, List[str]],
) -> pd.DataFrame:
    """
    Fast module scores from gzipped CSV (gene x cell counts).
    Computes mean(log1p(count)) across genes for each signature.

    The matrices in this project have quoted headers and gene symbols, e.g. "IL6".
    """
    sig_names = sorted({s for sets in gene_to_sets.values() for s in sets})
    union_genes = set(gene_to_sets.keys())

    with gzip.open(file_path, "rt") as f:
        header = f.readline().rstrip("\n")
        cols_raw = header.split(",")
        cols = [c.strip('"') for c in cols_raw]
        cell_cols_all = cols[1:]

        if cells_keep is None:
            keep_idx = None
            cell_cols = cell_cols_all
        else:
            keep_set = set(cells_keep)
            col_to_idx = {c: i for i, c in enumerate(cell_cols_all)}
            keep = [c for c in cell_cols_all if c in keep_set]
            keep_idx = np.array([col_to_idx[c] for c in keep], dtype=int)
            cell_cols = keep

        sums = {s: np.zeros(len(cell_cols), dtype=np.float64) for s in sig_names}
        counts = {s: 0 for s in sig_names}

        for line in f:
            if not line:
                continue
            try:
                gene_raw, rest = line.rstrip("\n").split(",", 1)
            except ValueError:
                continue
            gene = gene_raw.strip('"')
            if gene not in union_genes:
                continue

            vals = np.fromstring(rest, sep=",", dtype=np.float64)
            if keep_idx is not None:
                vals = vals[keep_idx]
            logvals = np.log1p(vals)

            for s in gene_to_sets[gene]:
                sums[s] += logvals
                counts[s] += 1

        out = pd.DataFrame(index=cell_cols)
        for s in sig_names:
            out[s] = sums[s] / counts[s] if counts[s] > 0 else np.nan
        out.index.name = "cell_id"
        return out


def compute_raw_scores_from_tsvgz_fast(
    file_path: str,
    cells_keep: Optional[Set[str]],
    gene_to_sets: Dict[str, List[str]],
    sep: str = "\t",
) -> pd.DataFrame:
    """
    Fast module scores from gzipped TSV (gene x cell expression, already log-scale).
    Computes mean(expression) across genes for each signature.
    """
    sig_names = sorted({s for sets in gene_to_sets.values() for s in sets})
    union_genes = set(gene_to_sets.keys())

    with gzip.open(file_path, "rt") as f:
        header = f.readline().rstrip("\n")
        cols = header.split(sep)
        cell_cols_all = cols[1:]

        if cells_keep is None:
            keep_idx = None
            cell_cols = cell_cols_all
        else:
            keep_set = set(cells_keep)
            col_to_idx = {c: i for i, c in enumerate(cell_cols_all)}
            keep = [c for c in cell_cols_all if c in keep_set]
            keep_idx = np.array([col_to_idx[c] for c in keep], dtype=int)
            cell_cols = keep

        sums = {s: np.zeros(len(cell_cols), dtype=np.float64) for s in sig_names}
        counts = {s: 0 for s in sig_names}

        for line in f:
            if not line:
                continue
            parts = line.rstrip("\n").split(sep, 1)
            if len(parts) < 2:
                continue
            gene, rest = parts[0], parts[1]
            if gene not in union_genes:
                continue

            vals = np.fromstring(rest, sep=sep, dtype=np.float64)
            if keep_idx is not None:
                vals = vals[keep_idx]

            for s in gene_to_sets[gene]:
                sums[s] += vals
                counts[s] += 1

        out = pd.DataFrame(index=cell_cols)
        for s in sig_names:
            out[s] = sums[s] / counts[s] if counts[s] > 0 else np.nan
        out.index.name = "cell_id"
        return out


def make_violin_panel(
    ax,
    data: pd.DataFrame,
    value_col: str,
    groups: List[str],
    palette: Dict[str, str],
    title: str,
    pvals_df: pd.DataFrame,
    dataset: str,
    donor_col: str = "sample_code",
    control_name: str = "Control",
    comparison_fmt: str = "{g} vs {control}",
):
    arrays_by_group = {}
    donor_meds = {}

    for g in groups:
        vals = data.loc[data["group"] == g, value_col].astype(float).dropna().values
        arrays_by_group[g] = vals
        if donor_col in data.columns:
            donor_meds[g] = (
                data.loc[data["group"] == g].groupby(donor_col)[value_col].median().dropna()
            )
        else:
            donor_meds[g] = pd.Series(dtype=float)

    positions_all = np.arange(1, len(groups) + 1)
    group_to_pos = {g: positions_all[i] for i, g in enumerate(groups)}

    groups_nonempty = [g for g in groups if len(arrays_by_group[g]) > 0]
    arrays = [arrays_by_group[g] for g in groups_nonempty]
    positions = [group_to_pos[g] for g in groups_nonempty]

    if len(groups_nonempty) > 0:
        parts = ax.violinplot(
            arrays,
            positions=positions,
            showmeans=False,
            showmedians=False,
            showextrema=False,
            widths=0.8,
        )
        for i, body in enumerate(parts["bodies"]):
            g = groups_nonempty[i]
            body.set_facecolor(palette[g])
            body.set_edgecolor("black")
            body.set_linewidth(1.5)
            body.set_alpha(1.0)

        ax.boxplot(
            arrays,
            positions=positions,
            widths=0.22,
            patch_artist=True,
            showfliers=False,
            boxprops=dict(facecolor="white", edgecolor="black", linewidth=1.2),
            medianprops=dict(color="black", linewidth=1.2),
            whiskerprops=dict(color="black", linewidth=1.2),
            capprops=dict(color="black", linewidth=1.2),
        )

        rng = np.random.default_rng(0)
        for g in groups_nonempty:
            meds = donor_meds[g]
            if len(meds) == 0:
                continue
            x = group_to_pos[g] + rng.uniform(-0.08, 0.08, size=len(meds))
            ax.scatter(
                x,
                meds.values,
                s=45,
                color=palette[g],
                edgecolor="black",
                linewidth=0.8,
                zorder=5,
            )

        all_vals = np.concatenate([arrays_by_group[g] for g in groups_nonempty])
        y_min = np.nanmin(all_vals)
        y_max = np.nanmax(all_vals)
        if y_min == y_max:
            y_min -= 1
            y_max += 1
        margin = 0.12 * (y_max - y_min)
        ax.set_ylim(y_min - 0.05 * (y_max - y_min), y_max + margin)
        y_text = y_max + 0.02 * (y_max - y_min) + 0.2 * margin
    else:
        y_text = 0.0

    for g in groups:
        if g == control_name:
            continue
        comp = comparison_fmt.format(g=g, control=control_name)
        row = pvals_df[
            (pvals_df["dataset"] == dataset)
            & (pvals_df["metric"] == value_col)
            & (pvals_df["comparison"] == comp)
        ]
        if len(row) == 0:
            continue
        p = row["p_value"].iloc[0]
        txt = f"{star(p)}\n{format_p(p)}"
        ax.text(group_to_pos[g], y_text, txt, ha="center", va="bottom", fontsize=8)

    ax.set_title(title, fontsize=12)
    ax.set_xticks(positions_all)
    ax.set_xticklabels(groups, rotation=45, ha="right")
    ax.set_ylabel("Module score", fontsize=11)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=10)
    ax.grid(False)


# -----------------------------
# Main
# -----------------------------
def main():
    # Gene sets (small lists)
    immune_genes = pd.read_csv("immune.csv").iloc[:, 0].astype(str).str.strip().tolist()
    apopt_genes = pd.read_csv("Germ_Apoptosis_Score.csv").iloc[:, 0].astype(str).str.strip().tolist()
    lcimm_genes = pd.read_csv("LC_Immature_Progenitor_Score.csv").iloc[:, 0].astype(str).str.strip().tolist()

    gene_to_sets: Dict[str, List[str]] = {}
    for g in immune_genes:
        gene_to_sets.setdefault(g, []).append("Cytokine_raw")
    for g in apopt_genes:
        gene_to_sets.setdefault(g, []).append("GermApoptosis_raw")
    for g in lcimm_genes:
        gene_to_sets.setdefault(g, []).append("LCimmature_raw")

    # -------------------------
    # GSE149512
    # -------------------------
    df149 = pd.read_csv("GSE149512_cellLevel_moduleScores_6metrics.csv")
    mapping = pd.read_csv("GSE149512_sample_mapping_corrected.csv")
    mapping149 = mapping[mapping["sample_code"].isin(df149["sample_code"].unique())].copy()

    # compute additional raw scores from each sample matrix
    t0 = time.time()
    score_parts = []
    for sample in sorted(df149["sample_code"].unique()):
        gsm = mapping149.loc[mapping149["sample_code"] == sample, "GSM"].iloc[0]
        mat_path = f"{gsm}_{sample}matrix.csv.gz"
        cells_keep = set(df149.loc[df149.sample_code == sample, "cell_id"])
        sc = compute_raw_scores_from_matrix_fast(mat_path, cells_keep, gene_to_sets)
        sc["sample_code"] = sample
        score_parts.append(sc)
    scores149 = pd.concat(score_parts, axis=0)
    print(f"[GSE149512] additional scores computed in {time.time()-t0:.1f}s")

    df149_full = (
        df149.set_index("cell_id")
        .join(scores149.drop(columns=["sample_code"]), how="left")
        .reset_index()
    )

    # z-scores (reference: Control cells within compartment)
    df149_full["Cytokine_z"], _, _ = zscore_by_reference(
        df149_full, "Cytokine_raw", (df149_full["celltype"] != "Germ") & (df149_full["group"] == "Control")
    )
    df149_full["GermApoptosis_z"], _, _ = zscore_by_reference(
        df149_full, "GermApoptosis_raw", (df149_full["celltype"] == "Germ") & (df149_full["group"] == "Control")
    )
    df149_full["LCimmature_z"], _, _ = zscore_by_reference(
        df149_full, "LCimmature_raw", (df149_full["celltype"] == "LC") & (df149_full["group"] == "Control")
    )

    # panels definition
    groups149 = ["Control", "AZFa", "iNOA", "KS"]
    palette149 = {"Control": "#1B9E77", "AZFa": "#E41A1C", "iNOA": "#6A3D9A", "KS": "#A6761D"}
    panels149 = [
        ("ST maturity", "ST_maturity", lambda d: d["celltype"] == "ST"),
        ("LC immaturity", "LCimmature_z", lambda d: d["celltype"] == "LC"),
        ("ST BTB integrity", "BTB", lambda d: d["celltype"] == "ST"),
        ("Somatic cytokine signaling", "Cytokine_z", lambda d: d["celltype"] != "Germ"),
        ("Somatic SASP", "SASP", lambda d: d["celltype"] != "Germ"),
        ("Germ apoptosis", "GermApoptosis_z", lambda d: d["celltype"] == "Germ"),
        ("ST lactate/glycolysis", "ST_lactate", lambda d: d["celltype"] == "ST"),
        ("SPC OXPHOS", "OXPHOS", lambda d: d["celltype"] == "Germ"),
    ]

    # p-values table (cell-level)
    pvals_records = []
    for title, col, mask_fn in panels149:
        sub = df149_full[mask_fn(df149_full)]
        ctrl = sub[sub["group"] == "Control"][col]
        for g in groups149[1:]:
            oth = sub[sub["group"] == g][col]
            p = mannwhitney_p(ctrl, oth)
            pvals_records.append(
                dict(
                    dataset="GSE149512",
                    panel=title,
                    metric=col,
                    comparison=f"{g} vs Control",
                    n_control=int(ctrl.notna().sum()),
                    n_group=int(oth.notna().sum()),
                    p_value=p,
                )
            )
    pvals149 = pd.DataFrame(pvals_records)

    # save intermediate data
    df149_full.to_csv("GSE149512_cellLevel_geneScores8_withCellPvalues_data.csv", index=False)
    pvals149.to_csv("GSE149512_geneScores8_cellLevel_pvalues.csv", index=False)

    # plot
    fig, axes = plt.subplots(3, 3, figsize=(14, 10))
    fig.suptitle("Gene score (GSE149512)", fontsize=18, fontweight="bold", y=0.98)
    axes_flat = axes.flatten()

    for idx, (title, col, mask_fn) in enumerate(panels149):
        ax = axes_flat[idx]
        sub = df149_full[mask_fn(df149_full)].copy()
        make_violin_panel(
            ax,
            sub,
            col,
            groups149,
            palette149,
            title,
            pvals149,
            dataset="GSE149512",
            donor_col="sample_code",
            control_name="Control",
            comparison_fmt="{g} vs {control}",
        )

    for j in range(len(panels149), 9):
        axes_flat[j].axis("off")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig("ExternalValidation_GSE149512_GeneScores_8panels_cellPvalues.pdf")
    fig.savefig("ExternalValidation_GSE149512_GeneScores_8panels_cellPvalues.png", dpi=300)
    plt.close(fig)

    # -------------------------
    # GSE235321
    # -------------------------
    df235 = pd.read_csv("GSE235321_cellLevel_moduleScores_6metrics.csv")

    # infer donor/sample from cell_id: use first two underscore-separated fields
    df235["sample_code"] = df235["cell_id"].astype(str).str.split("_").str[0:2].str.join("_")

    score_parts = []
    for grp in ["Normal", "NOA1", "NOA2"]:
        expr_path = f"GSE235321_exp.log2TPM.{grp}.tsv.gz"
        cells_keep = set(df235.loc[df235.group == grp, "cell_id"])
        sc = compute_raw_scores_from_tsvgz_fast(expr_path, cells_keep, gene_to_sets, sep="\t")
        sc["group"] = grp
        score_parts.append(sc)
    scores235 = pd.concat(score_parts, axis=0)

    df235_full = df235.set_index("cell_id").join(scores235.drop(columns=["group"]), how="left").reset_index()

    # z-scores (reference: Normal cells within compartment)
    df235_full["Cytokine_z"], _, _ = zscore_by_reference(
        df235_full, "Cytokine_raw", (df235_full["celltype"] != "Germ") & (df235_full["group"] == "Normal")
    )
    df235_full["GermApoptosis_z"], _, _ = zscore_by_reference(
        df235_full, "GermApoptosis_raw", (df235_full["celltype"] == "Germ") & (df235_full["group"] == "Normal")
    )
    df235_full["LCimmature_z"], _, _ = zscore_by_reference(
        df235_full, "LCimmature_raw", (df235_full["celltype"] == "LC") & (df235_full["group"] == "Normal")
    )

    groups235 = ["Normal", "NOA1", "NOA2"]
    palette235 = {"Normal": "#1B9E77", "NOA1": "#E41A1C", "NOA2": "#6A3D9A"}
    panels235 = [
        ("ST maturity", "ST_maturity", lambda d: d["celltype"] == "ST"),
        ("LC immaturity", "LCimmature_z", lambda d: d["celltype"] == "LC"),
        ("ST BTB integrity", "BTB", lambda d: d["celltype"] == "ST"),
        ("Somatic cytokine signaling", "Cytokine_z", lambda d: d["celltype"] != "Germ"),
        ("Somatic SASP", "SASP", lambda d: d["celltype"] != "Germ"),
        ("Germ apoptosis", "GermApoptosis_z", lambda d: d["celltype"] == "Germ"),
        ("ST lactate/glycolysis", "ST_lactate", lambda d: d["celltype"] == "ST"),
        ("SPC OXPHOS", "OXPHOS", lambda d: d["celltype"] == "Germ"),
    ]

    pvals_records = []
    for title, col, mask_fn in panels235:
        sub = df235_full[mask_fn(df235_full)]
        ctrl = sub[sub["group"] == "Normal"][col]
        for g in groups235[1:]:
            oth = sub[sub["group"] == g][col]
            p = mannwhitney_p(ctrl, oth)
            pvals_records.append(
                dict(
                    dataset="GSE235321",
                    panel=title,
                    metric=col,
                    comparison=f"{g} vs Normal",
                    n_control=int(ctrl.notna().sum()),
                    n_group=int(oth.notna().sum()),
                    p_value=p,
                )
            )
    pvals235 = pd.DataFrame(pvals_records)

    df235_full.to_csv("GSE235321_cellLevel_geneScores8_withCellPvalues_data.csv", index=False)
    pvals235.to_csv("GSE235321_geneScores8_cellLevel_pvalues.csv", index=False)

    fig, axes = plt.subplots(3, 3, figsize=(14, 10))
    fig.suptitle("Gene score (GSE235321)", fontsize=18, fontweight="bold", y=0.98)
    axes_flat = axes.flatten()

    for idx, (title, col, mask_fn) in enumerate(panels235):
        ax = axes_flat[idx]
        sub = df235_full[mask_fn(df235_full)].copy()
        make_violin_panel(
            ax,
            sub,
            col,
            groups235,
            palette235,
            title,
            pvals235,
            dataset="GSE235321",
            donor_col="sample_code",
            control_name="Normal",
            comparison_fmt="{g} vs {control}",
        )

    for j in range(len(panels235), 9):
        axes_flat[j].axis("off")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig("ExternalValidation_GSE235321_GeneScores_8panels_cellPvalues.pdf")
    fig.savefig("ExternalValidation_GSE235321_GeneScores_8panels_cellPvalues.png", dpi=300)
    plt.close(fig)

    print("Done.")


if __name__ == "__main__":
    main()
