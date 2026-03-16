#!/usr/bin/env python3
from __future__ import annotations

import math
import re
import shutil
import subprocess
import tarfile
import zipfile
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parents[1]  # .../source data
ROOT = BASE.parent
OUT_DIR = BASE / "fig8_source_data"
RAW_SUPPORT_DIR = OUT_DIR / "raw_support"

FIG8_AI = ROOT / "Figure 8_20251210.ai"
FIG8_PDF = ROOT / "Figure 8_20251210.pdf"
FIG8_PNG = ROOT / "Figure 8_20251210.png"

# Panel A
CELL_ANNOT = BASE / "figure1" / "fig1_source_data" / "cell_annotation_from_clusters.csv"
INTEGRATION_MULTI_DIR = (
    ROOT / "原来作图的数据" / "汇总" / "P20073103_B3_rmRP_integration_report" / "results" / "multi"
)
PANEL_A_FILES = [
    INTEGRATION_MULTI_DIR / "integration_group_CellsPerCluster.xls",
    INTEGRATION_MULTI_DIR / "integration_group_PercentPerCluster.xls",
    INTEGRATION_MULTI_DIR / "integration_sample_CellsPerCluster.xls",
    INTEGRATION_MULTI_DIR / "integration_sample_PercentPerCluster.xls",
    INTEGRATION_MULTI_DIR / "integration_group_celltype_Summary.txt",
    INTEGRATION_MULTI_DIR / "integration_sample_celltype_Summary.txt",
    INTEGRATION_MULTI_DIR / "integration_umap_samples.pdf",
    INTEGRATION_MULTI_DIR / "integration_umap_ugroups.pdf",
]

# Panel B
PANEL_B_INPUT = ROOT / "原来作图的数据" / "5" / "生精细胞相关表达" / "1.csv"
PANEL_B_SCRIPT = ROOT / "原来作图的数据" / "5" / "生精细胞相关表达" / "比例图.R"
SUB_GERM_DIR = (
    ROOT
    / "原来作图的数据"
    / "汇总"
    / "P20073103_B3_rmRP_integration_report"
    / "results"
    / "sub"
    / "GermCells_nobatch"
)
PANEL_B_FILES = [
    SUB_GERM_DIR / "sub_group_CellsPerCluster.xls",
    SUB_GERM_DIR / "sub_group_PercentPerCluster.xls",
    SUB_GERM_DIR / "sub_group_celltype_Summary.txt",
    SUB_GERM_DIR / "sub_umap_samples.pdf",
]

# Panel C/D
PANEL_CD_DIR = ROOT / "原来作图的数据" / "汇总" / "37.shengzhi_diff_gene_up_down"
PANEL_C_PDF = PANEL_CD_DIR / "P20073103_CtrlvsOA_down_DEG_heatmap.pdf"  # down in Ctrl => up in OA
PANEL_D_PDF = PANEL_CD_DIR / "P20073103_CtrlvsOA_up_DEG_heatmap.pdf"    # up in Ctrl => down in OA
PANEL_C_PNG = PANEL_CD_DIR / "P20073103_CtrlvsOA_down_DEG_heatmap.png"
PANEL_D_PNG = PANEL_CD_DIR / "P20073103_CtrlvsOA_up_DEG_heatmap.png"
DIFF_DIR = ROOT / "原来作图的数据" / "汇总" / "11.生殖细胞组间差异分析" / "Diff"
SUBTYPE_ORDER = [
    "Early_primary_SPCs",
    "Elongated_Spermatids",
    "Late_primary_SPCs",
    "Round_Spermatids",
    "Sperm",
    "SPGs",
    "SSCs",
]

# Panel E/F
ENRICH_BASE = ROOT / "原来作图的数据" / "汇总" / "11.生殖细胞组间差异分析" / "Enrich"
ROUND_ENRICH_DIR = ENRICH_BASE / "enrich_Round_Spermatids_CtrlvsOA" / "P20073103"
ELONG_ENRICH_DIR = ENRICH_BASE / "enrich_Elongated_Spermatids_CtrlvsOA" / "P20073103"

ROUND_UP_FULL = ROUND_ENRICH_DIR / "P20073103_up.GOALL_enrichment.xls"
ROUND_DOWN_FULL = ROUND_ENRICH_DIR / "P20073103_down.GOALL_enrichment.xls"
ROUND_UP_SIG = ROUND_ENRICH_DIR / "P20073103_up.GOALL_enrichment_sig.xls"
ROUND_DOWN_SIG = ROUND_ENRICH_DIR / "P20073103_down.GOALL_enrichment_sig.xls"

ELONG_UP_FULL = ELONG_ENRICH_DIR / "P20073103_up.GOALL_enrichment.xls"
ELONG_DOWN_FULL = ELONG_ENRICH_DIR / "P20073103_down.GOALL_enrichment.xls"
ELONG_UP_SIG = ELONG_ENRICH_DIR / "P20073103_up.GOALL_enrichment_sig.xls"
ELONG_DOWN_SIG = ELONG_ENRICH_DIR / "P20073103_down.GOALL_enrichment_sig.xls"

PANEL_E_DOWN_IN_OA_TERMS = [
    "Regulation of RNA splicing",
    "RNA splicing",
    "Regulation of mRNA processing",
    "Regulation of cell-substrate adhesion",
    "ATP metabolic process",
]
PANEL_E_UP_IN_OA_TERMS = [
    "Cilium movement",
    "Fertilization",
    "Spermatid development",
    "Sperm-egg recognition",
    "Germ cell development",
]
PANEL_F_DOWN_IN_OA_TERMS = [
    "Chromosome segregation",
    "RNA splicing",
    "Aerobic respiration",
    "Oxidative phosphorylation",
    "ATP metabolic process",
]
PANEL_F_UP_IN_OA_TERMS = [
    "Motile cilium",
    "Sperm flagellum",
    "Microtubule-based movement",
    "Cilium movement",
    "Fertilization",
]

# Panel G/H candidate support
INFERTILITY_GENE_LIST = ROOT / "原来作图的数据" / "5" / "不育基因表达" / "不育相关基因.csv"
AVG_EXPR_GERM = ROOT / "原来作图的数据" / "average_expression_matrix生精.csv"
PANEL_G_HEATMAP_PDF_CANDIDATES = [
    ROOT / "原来作图的数据" / "无精症" / "heatmap_export(1).pdf",
    ROOT / "原来作图的数据" / "无精症" / "heatmap_export.pdf",
]
PANEL_G_DRGS_XLSX = ROOT / "原来作图的数据" / "2" / "reactome" / "DRGs all-ST.xlsx"
PANEL_G_DRGS_CSV = ROOT / "原来作图的数据" / "2" / "reactome" / "DRGs all-ST.csv"
PANEL_G_ST_MATRIX_XLSX = (
    ROOT / "原来作图的数据" / "2" / "reactome" / "未命名文件夹" / "average_expression_matrix_ST.xlsx"
)
PANEL_G_ST_MATRIX_CSV = (
    ROOT / "原来作图的数据" / "2" / "reactome" / "未命名文件夹" / "average_expression_matrix_ST.csv"
)
PANEL_GH_ARCHIVE = OUT_DIR / "raw_support" / "OA-CTRL热图.gz"

PANEL_H_GENE_LIST = ["CCDC65", "DNAH8", "PIH1D3", "CATSPER3"]
PANEL_H_VIOLIN_DIR_1 = ROOT / "原来作图的数据" / "无精症" / "显著 2"
PANEL_H_VIOLIN_DIR_2 = ROOT / "原来作图的数据" / "无精症" / "violin_plots2"
PANEL_H_VIOLIN_DIR_3 = ROOT / "原来作图的数据" / "5" / "不育基因表达"
PANEL_H_DIRECT_CSV = OUT_DIR / "raw_support" / "fig8h.csv"

# Panel I
PANEL_I_INPUT = ROOT / "原来作图的数据" / "5" / "体细胞表达" / "2.csv"
PANEL_I_SCRIPT = ROOT / "原来作图的数据" / "5" / "体细胞表达" / "比例图.R"

# Panel J
PANEL_J_ZIP = ROOT / "原来作图的数据" / "无精症" / "OA-CTRL.zip"
PANEL_J_MEMBERS = {
    "mature_score": "OA-CTRL/genescore-mature.csv",
    "immature_score": "OA-CTRL/genescore-immature.csv",
    "mature_expression": "OA-CTRL/expression.csv",
    "immature_expression": "OA-CTRL/expression-immature.csv",
}

GROUP_ORDER_CTRL_OA = ["Ctrl", "OA"]
GROUP_ORDER_ALL = ["Ctrl", "OA", "AZFc_Del", "iNOA_B", "iNOA_S", "KS"]
MYOID_RAW_CLUSTERS = {11, 14, 27, 33}


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


def safe_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def neg_log10(value) -> float | None:
    try:
        v = float(value)
    except Exception:
        return None
    if not math.isfinite(v) or v <= 0:
        return None
    return -math.log10(v)


def summarize_score(df: pd.DataFrame, group_col: str, score_col: str, order: list[str]) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=[group_col, "n_cells", "mean", "median", "q1", "q3", "min", "max"])
    g = df.groupby(group_col, observed=False)[score_col]
    out = (
        g.agg(n_cells="count", mean="mean", median="median", min="min", max="max")
        .reset_index()
        .merge(g.quantile(0.25).reset_index(name="q1"), on=group_col, how="left")
        .merge(g.quantile(0.75).reset_index(name="q3"), on=group_col, how="left")
    )
    out[group_col] = pd.Categorical(out[group_col], categories=order, ordered=True)
    return out.sort_values(group_col).reset_index(drop=True)


def pdftotext_lines(pdf_path: Path) -> list[str]:
    if not pdf_path.exists():
        return []
    try:
        proc = subprocess.run(["pdftotext", str(pdf_path), "-"], capture_output=True, text=True, check=True)
    except Exception:
        return []
    out: list[str] = []
    for line in proc.stdout.splitlines():
        text = re.sub(r"\s+", " ", line.strip())
        if text:
            out.append(text)
    return out


def extract_heatmap_numbers(pdf_path: Path) -> list[int]:
    nums: list[int] = []
    for line in pdftotext_lines(pdf_path):
        if re.fullmatch(r"\d+", line):
            nums.append(int(line))
    return nums


def extract_panel_g_genes_and_groups(pdf_path: Path) -> tuple[list[str], list[str]]:
    lines = pdftotext_lines(pdf_path)
    genes: list[str] = []
    groups: list[str] = []

    group_alias_map = {
        "ECM_Remodeling": "ECM_Remodeling",
        "Oxidative_Stress": "Oxidative_Stress",
        "Obstruction_Genes": "Obstruction_Genes",
        "Mito_Respiration": "Mito_Respiration",
    }

    for raw in lines:
        line = raw.strip()
        if not line:
            continue

        if line in group_alias_map and group_alias_map[line] not in groups:
            groups.append(group_alias_map[line])
            continue
        norm = line.replace("−", "-")
        if "Loss of epididymal" in norm and "proteins" in norm:
            canon = "Loss of epididymal-maturation proteins"
            if canon not in groups:
                groups.append(canon)
            continue

        # Gene-like tokens in this exported heatmap text are all-uppercase alnum strings.
        if re.fullmatch(r"[A-Z0-9-]+", line) and len(line) >= 4:
            if line not in {"HEATMAP"} and line not in genes:
                genes.append(line)

    return genes, groups


def pick_terms(df: pd.DataFrame, terms: list[str]) -> pd.DataFrame:
    if df.empty or "Description" not in df.columns:
        return pd.DataFrame(columns=[
            "term_order",
            "term_requested",
            "Description",
            "match_type",
            "found",
            "pvalue",
            "neg_log10_pvalue",
            "Count",
            "GeneRatio",
            "BgRatio",
        ])

    rows: list[dict] = []
    desc = df["Description"].astype(str)
    desc_low = desc.str.lower()

    for idx, term in enumerate(terms, start=1):
        exact = df[desc_low == term.lower()]
        match_type = "exact"
        hit = exact
        if hit.empty:
            hit = df[desc_low.str.contains(term.lower(), regex=False)]
            match_type = "contains"

        if hit.empty:
            rows.append(
                {
                    "term_order": idx,
                    "term_requested": term,
                    "Description": "",
                    "match_type": "not_found",
                    "found": False,
                    "pvalue": None,
                    "neg_log10_pvalue": None,
                    "Count": None,
                    "GeneRatio": "",
                    "BgRatio": "",
                }
            )
            continue

        row = hit.copy()
        row["pvalue"] = safe_num(row["pvalue"])
        best = row.sort_values("pvalue", ascending=True).iloc[0]
        rows.append(
            {
                "term_order": idx,
                "term_requested": term,
                "Description": best.get("Description", ""),
                "match_type": match_type,
                "found": True,
                "pvalue": best.get("pvalue", None),
                "neg_log10_pvalue": neg_log10(best.get("pvalue", None)),
                "Count": best.get("Count", None),
                "GeneRatio": best.get("GeneRatio", ""),
                "BgRatio": best.get("BgRatio", ""),
            }
        )

    out = pd.DataFrame(rows)
    return out.sort_values("term_order").reset_index(drop=True)


def read_zip_csv(zip_path: Path, member: str) -> pd.DataFrame:
    if not zip_path.exists():
        return pd.DataFrame()
    try:
        with zipfile.ZipFile(zip_path) as zf:
            with zf.open(member) as f:
                return pd.read_csv(f)
    except Exception:
        return pd.DataFrame()


def read_tar_gz_csvs(tar_gz_path: Path) -> dict[str, pd.DataFrame]:
    tables: dict[str, pd.DataFrame] = {}
    if not tar_gz_path.exists():
        return tables
    try:
        with tarfile.open(tar_gz_path, "r:gz") as tf:
            for member in tf.getmembers():
                if not member.isfile():
                    continue
                if not member.name.lower().endswith(".csv"):
                    continue
                base = Path(member.name).name
                f = tf.extractfile(member)
                if f is None:
                    continue
                try:
                    tables[base] = pd.read_csv(f)
                except Exception:
                    continue
    except Exception:
        return {}
    return tables


def map_panel_a_celltype(cell_type: str, raw_cluster: int) -> str:
    if cell_type in {"SPCs", "SPGs", "SSCs", "Spermatids", "LCs", "ECs", "STs"}:
        return cell_type
    if cell_type == "Lym":
        return "Immune cells"
    if cell_type == "Myeloid":
        return "Myoid" if int(raw_cluster) in MYOID_RAW_CLUSTERS else "Immune cells"
    return str(cell_type)


def build_panel_a() -> dict:
    generated: list[str] = []
    raw_inputs: list[str] = []

    for src in [FIG8_AI, FIG8_PDF, FIG8_PNG, CELL_ANNOT, *PANEL_A_FILES]:
        rel = copy_to_raw_support(src, f"Fig8A_{src.name}")
        if rel:
            raw_inputs.append(rel)

    if not CELL_ANNOT.exists():
        note = pd.DataFrame(
            [{"panel": "A", "status": "missing", "note": "cell_annotation_from_clusters.csv not found"}]
        )
        note.to_csv(OUT_DIR / "Fig8A_source_note.csv", index=False)
        generated.append("Fig8A_source_note.csv")
        return {
            "subpanel": "UMAP of samples and cell types (Ctrl vs OA)",
            "generated": generated,
            "raw": raw_inputs,
            "code": ["../fig8_code/fig8_source_data.py"],
            "status": "missing",
            "availability_note": "cell-level UMAP table missing",
        }

    cols = ["cell_id", "raw_cluster", "UMAP1", "UMAP2", "cell_type", "sample_id", "group"]
    df = pd.read_csv(CELL_ANNOT, usecols=cols)
    df = df[df["group"].isin(GROUP_ORDER_CTRL_OA)].copy()
    df["raw_cluster"] = safe_num(df["raw_cluster"]).fillna(-1).astype(int)
    df["panel_celltype"] = [map_panel_a_celltype(ct, rc) for ct, rc in zip(df["cell_type"], df["raw_cluster"])]

    germ_types = {"SPCs", "SPGs", "SSCs", "Spermatids"}
    df["panel_major"] = df["panel_celltype"].apply(lambda x: "Germ cells" if x in germ_types else "Somatic cells")

    panel_order = [
        "SPCs",
        "SPGs",
        "SSCs",
        "Spermatids",
        "Immune cells",
        "Myoid",
        "LCs",
        "ECs",
        "STs",
    ]
    df["panel_celltype"] = pd.Categorical(df["panel_celltype"], categories=panel_order, ordered=True)

    keep = ["cell_id", "group", "sample_id", "raw_cluster", "UMAP1", "UMAP2", "cell_type", "panel_major", "panel_celltype"]
    df = df[keep].sort_values(["group", "panel_celltype", "cell_id"]).reset_index(drop=True)
    df.to_csv(OUT_DIR / "Fig8A_ctrl_oa_umap_cells.csv", index=False)
    generated.append("Fig8A_ctrl_oa_umap_cells.csv")

    legend_group = df.groupby("group", observed=False).size().reset_index(name="n_cells")
    legend_group["group"] = pd.Categorical(legend_group["group"], categories=GROUP_ORDER_CTRL_OA, ordered=True)
    legend_group = legend_group.sort_values("group").reset_index(drop=True)
    legend_group.to_csv(OUT_DIR / "Fig8A_legend_group_counts.csv", index=False)
    generated.append("Fig8A_legend_group_counts.csv")

    legend_celltype = df.groupby("panel_celltype", observed=False).size().reset_index(name="n_cells")
    legend_celltype["panel_celltype"] = pd.Categorical(
        legend_celltype["panel_celltype"], categories=panel_order, ordered=True
    )
    legend_celltype = legend_celltype.sort_values("panel_celltype").reset_index(drop=True)
    legend_celltype.to_csv(OUT_DIR / "Fig8A_legend_celltype_counts.csv", index=False)
    generated.append("Fig8A_legend_celltype_counts.csv")

    group_cell = (
        df.groupby(["group", "panel_celltype"], observed=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["group", "panel_celltype"])
        .reset_index(drop=True)
    )
    group_cell.to_csv(OUT_DIR / "Fig8A_group_celltype_counts_long.csv", index=False)
    generated.append("Fig8A_group_celltype_counts_long.csv")

    note = pd.DataFrame(
        [
            {
                "panel": "A",
                "note": "UMAP coordinates are from cell_annotation_from_clusters.csv; Myeloid raw clusters were split into Myoid (11/14/27/33) and Immune cells for Figure 8 legend compatibility.",
            }
        ]
    )
    note.to_csv(OUT_DIR / "Fig8A_source_note.csv", index=False)
    generated.append("Fig8A_source_note.csv")

    return {
        "subpanel": "UMAP of samples and cell types (Ctrl vs OA)",
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../fig8_code/fig8_source_data.py"],
        "status": "available",
        "availability_note": "cell-level UMAP coordinates and legend counts are available",
    }


def build_panel_b() -> dict:
    generated: list[str] = []
    raw_inputs: list[str] = []
    code_refs = ["../fig8_code/fig8_source_data.py"]

    for src in [PANEL_B_INPUT, PANEL_B_SCRIPT, *PANEL_B_FILES]:
        rel = copy_to_raw_support(src, f"Fig8B_{src.name}")
        if rel:
            raw_inputs.append(rel)
            if src == PANEL_B_SCRIPT:
                code_refs.append(rel)

    if not PANEL_B_INPUT.exists():
        note = pd.DataFrame(
            [{"panel": "B", "status": "missing", "note": "1.csv not found in 生精细胞相关表达"}]
        )
        note.to_csv(OUT_DIR / "Fig8B_source_note.csv", index=False)
        generated.append("Fig8B_source_note.csv")
        return {
            "subpanel": "Germ-cell composition stacked bar (Ctrl vs OA)",
            "generated": generated,
            "raw": raw_inputs,
            "code": code_refs,
            "status": "missing",
            "availability_note": "germ composition input missing",
        }

    inp = pd.read_csv(PANEL_B_INPUT).rename(columns={"Cluster": "cluster", "Sample": "group", "Number": "count"})
    inp = inp[inp["group"].isin(GROUP_ORDER_CTRL_OA)].copy()
    inp["count"] = safe_num(inp["count"]).fillna(0).astype(int)
    cluster_order = ["Sperm", "Elongated_Spermatids", "Round_Spermatids", "Late_primary_SPCs", "Early_primary_SPCs", "SPGs", "SSCs"]
    inp["cluster"] = pd.Categorical(inp["cluster"], categories=cluster_order, ordered=True)
    inp["group"] = pd.Categorical(inp["group"], categories=GROUP_ORDER_CTRL_OA, ordered=True)
    inp = inp.sort_values(["group", "cluster"]).reset_index(drop=True)
    inp.to_csv(OUT_DIR / "Fig8B_germ_composition_input.csv", index=False)
    generated.append("Fig8B_germ_composition_input.csv")

    total = inp.groupby("group", observed=False)["count"].transform("sum")
    pct = inp.copy()
    pct["percent"] = pct["count"] / total
    pct.to_csv(OUT_DIR / "Fig8B_germ_composition_percent_long.csv", index=False)
    generated.append("Fig8B_germ_composition_percent_long.csv")

    wide = pct.pivot(index="group", columns="cluster", values="percent").reset_index()
    wide.to_csv(OUT_DIR / "Fig8B_germ_composition_percent_wide.csv", index=False)
    generated.append("Fig8B_germ_composition_percent_wide.csv")

    if (SUB_GERM_DIR / "sub_group_CellsPerCluster.xls").exists():
        sub_cnt = pd.read_csv(SUB_GERM_DIR / "sub_group_CellsPerCluster.xls", sep="\t")
        long = sub_cnt.melt(id_vars=["cluster"], var_name="group", value_name="count")
        long = long[long["group"].isin(GROUP_ORDER_CTRL_OA)].copy()
        long["count"] = safe_num(long["count"]).fillna(0)
        long.to_csv(OUT_DIR / "Fig8B_subgroup_counts_from_integration_summary.csv", index=False)
        generated.append("Fig8B_subgroup_counts_from_integration_summary.csv")

    return {
        "subpanel": "Germ-cell composition stacked bar (Ctrl vs OA)",
        "generated": generated,
        "raw": raw_inputs,
        "code": code_refs,
        "status": "available",
        "availability_note": "direct input count table and plotting script are available",
    }


def build_overlap_matrix(gene_sets: dict[str, set[str]], order: list[str]) -> pd.DataFrame:
    rows: list[dict] = []
    for a in order:
        row = {"subtype": a}
        for b in order:
            row[b] = len(gene_sets.get(a, set()) & gene_sets.get(b, set()))
        rows.append(row)
    return pd.DataFrame(rows)


def build_panel_cd() -> tuple[dict, dict]:
    generated: list[str] = []
    raw_inputs: list[str] = []

    for src in [PANEL_C_PDF, PANEL_C_PNG, PANEL_D_PDF, PANEL_D_PNG]:
        rel = copy_to_raw_support(src, f"Fig8CD_{src.name}")
        if rel:
            raw_inputs.append(rel)

    diff_tables: dict[str, Path] = {}
    file_rows = []
    for subtype in SUBTYPE_ORDER:
        f = DIFF_DIR / f"diff_{subtype}" / "P20073103_CtrlvsOA.diffexpressed.xls"
        diff_tables[subtype] = f
        rel = copy_to_raw_support(f, f"Fig8CD_{subtype}_CtrlvsOA.diffexpressed.xls")
        if rel:
            raw_inputs.append(rel)
        file_rows.append({"subtype": subtype, "file": str(f), "exists": f.exists(), "raw_support_copy": rel})
    pd.DataFrame(file_rows).to_csv(OUT_DIR / "Fig8CD_ctrlvsOA_diff_file_registry.csv", index=False)
    generated.append("Fig8CD_ctrlvsOA_diff_file_registry.csv")

    gene_sets_up_oa: dict[str, set[str]] = {}
    gene_sets_down_oa: dict[str, set[str]] = {}
    up_rows: list[pd.DataFrame] = []
    down_rows: list[pd.DataFrame] = []

    for subtype in SUBTYPE_ORDER:
        path = diff_tables[subtype]
        if not path.exists():
            gene_sets_up_oa[subtype] = set()
            gene_sets_down_oa[subtype] = set()
            continue

        df = pd.read_csv(path, sep="\t")
        for col in ["avg_logFC", "p_val_adj", "p_val", "Ctrl", "OA"]:
            if col in df.columns:
                df[col] = safe_num(df[col])

        up = df[(df["p_val_adj"] < 0.05) & (df["avg_logFC"] < 0)].copy()   # up in OA
        down = df[(df["p_val_adj"] < 0.05) & (df["avg_logFC"] > 0)].copy()  # down in OA

        gene_sets_up_oa[subtype] = set(up["gene_id"].astype(str))
        gene_sets_down_oa[subtype] = set(down["gene_id"].astype(str))

        keep_cols = [c for c in ["gene_id", "avg_logFC", "p_val", "p_val_adj", "Ctrl", "OA", "pct.1", "pct.2"] if c in up.columns]
        if keep_cols:
            up2 = up[keep_cols].copy()
            up2.insert(0, "subtype", subtype)
            up_rows.append(up2)

            down2 = down[keep_cols].copy()
            down2.insert(0, "subtype", subtype)
            down_rows.append(down2)

    if up_rows:
        up_df = pd.concat(up_rows, ignore_index=True)
        up_df = up_df.sort_values(["subtype", "gene_id"]).reset_index(drop=True)
    else:
        up_df = pd.DataFrame(columns=["subtype", "gene_id"])
    up_df.to_csv(OUT_DIR / "Fig8C_up_in_OA_genes_from_diff_long.csv", index=False)
    generated.append("Fig8C_up_in_OA_genes_from_diff_long.csv")

    if down_rows:
        down_df = pd.concat(down_rows, ignore_index=True)
        down_df = down_df.sort_values(["subtype", "gene_id"]).reset_index(drop=True)
    else:
        down_df = pd.DataFrame(columns=["subtype", "gene_id"])
    down_df.to_csv(OUT_DIR / "Fig8D_down_in_OA_genes_from_diff_long.csv", index=False)
    generated.append("Fig8D_down_in_OA_genes_from_diff_long.csv")

    overlap_up = build_overlap_matrix(gene_sets_up_oa, SUBTYPE_ORDER)
    overlap_up.to_csv(OUT_DIR / "Fig8CD_overlap_matrix_up_in_OA.csv", index=False)
    generated.append("Fig8CD_overlap_matrix_up_in_OA.csv")

    overlap_down = build_overlap_matrix(gene_sets_down_oa, SUBTYPE_ORDER)
    overlap_down.to_csv(OUT_DIR / "Fig8CD_overlap_matrix_down_in_OA.csv", index=False)
    generated.append("Fig8CD_overlap_matrix_down_in_OA.csv")

    summary_rows = []
    up_intersection = set.intersection(*[gene_sets_up_oa[s] for s in SUBTYPE_ORDER]) if SUBTYPE_ORDER else set()
    down_intersection = set.intersection(*[gene_sets_down_oa[s] for s in SUBTYPE_ORDER]) if SUBTYPE_ORDER else set()
    for subtype in SUBTYPE_ORDER:
        other_up = set().union(*[gene_sets_up_oa[s] for s in SUBTYPE_ORDER if s != subtype])
        other_down = set().union(*[gene_sets_down_oa[s] for s in SUBTYPE_ORDER if s != subtype])
        summary_rows.append(
            {
                "subtype": subtype,
                "n_up_in_OA": len(gene_sets_up_oa[subtype]),
                "n_up_in_OA_unique": len(gene_sets_up_oa[subtype] - other_up),
                "n_down_in_OA": len(gene_sets_down_oa[subtype]),
                "n_down_in_OA_unique": len(gene_sets_down_oa[subtype] - other_down),
            }
        )
    summary_rows.append(
        {
            "subtype": "Shared_all_7_subtypes",
            "n_up_in_OA": len(up_intersection),
            "n_up_in_OA_unique": None,
            "n_down_in_OA": len(down_intersection),
            "n_down_in_OA_unique": None,
        }
    )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_DIR / "Fig8CD_set_size_summary_from_diff.csv", index=False)
    generated.append("Fig8CD_set_size_summary_from_diff.csv")

    # Direct panel numbers from original DEG heatmap PDFs.
    c_nums = extract_heatmap_numbers(PANEL_C_PDF)
    d_nums = extract_heatmap_numbers(PANEL_D_PDF)

    pdf_rows = []
    for panel, nums, meaning in [
        ("C", c_nums, "Upregulated in OA (from CtrlvsOA_down_DEG_heatmap.pdf)"),
        ("D", d_nums, "Downregulated in OA (from CtrlvsOA_up_DEG_heatmap.pdf)"),
    ]:
        shared = nums[0] if len(nums) >= 1 else None
        spec = nums[1:8] if len(nums) >= 8 else []
        pdf_rows.append(
            {
                "panel": panel,
                "meaning": meaning,
                "shared_count": shared,
                "numbers_raw": ",".join(str(x) for x in nums),
            }
        )
        for subtype, n in zip(SUBTYPE_ORDER, spec):
            pdf_rows.append(
                {
                    "panel": panel,
                    "meaning": f"{meaning} - subtype-specific",
                    "shared_count": None,
                    "subtype": subtype,
                    "subtype_specific_count": n,
                    "numbers_raw": ",".join(str(x) for x in nums),
                }
            )

    pdf_df = pd.DataFrame(pdf_rows)
    pdf_df.to_csv(OUT_DIR / "Fig8CD_panel_numbers_extracted_from_heatmap_pdf.csv", index=False)
    generated.append("Fig8CD_panel_numbers_extracted_from_heatmap_pdf.csv")

    note = pd.DataFrame(
        [
            {
                "panel": "C/D",
                "pdf_panelC_shared": c_nums[0] if len(c_nums) >= 1 else None,
                "pdf_panelD_shared": d_nums[0] if len(d_nums) >= 1 else None,
                "derived_shared_up_in_OA_from_diff": len(up_intersection),
                "derived_shared_down_in_OA_from_diff": len(down_intersection),
                "note": "Original panel C/D counts are extracted directly from summary DEG heatmap PDFs. Diff tables are provided as support tables; the unpublished row-binning step used to build the final heatmap matrix is not available in current folders.",
            }
        ]
    )
    note.to_csv(OUT_DIR / "Fig8CD_source_note.csv", index=False)
    generated.append("Fig8CD_source_note.csv")

    panel_c_meta = {
        "subpanel": "Number of DEGs upregulated in OA (heatmap-style overlap)",
        "generated": [
            "Fig8CD_ctrlvsOA_diff_file_registry.csv",
            "Fig8C_up_in_OA_genes_from_diff_long.csv",
            "Fig8CD_overlap_matrix_up_in_OA.csv",
            "Fig8CD_set_size_summary_from_diff.csv",
            "Fig8CD_panel_numbers_extracted_from_heatmap_pdf.csv",
            "Fig8CD_source_note.csv",
        ],
        "raw": raw_inputs,
        "code": ["../fig8_code/fig8_source_data.py"],
        "status": "available",
        "availability_note": "Panel-level counts are directly available from summary DEG heatmap PDFs; supporting subtype diff tables are archived",
    }

    panel_d_meta = {
        "subpanel": "Number of DEGs downregulated in OA (heatmap-style overlap)",
        "generated": [
            "Fig8CD_ctrlvsOA_diff_file_registry.csv",
            "Fig8D_down_in_OA_genes_from_diff_long.csv",
            "Fig8CD_overlap_matrix_down_in_OA.csv",
            "Fig8CD_set_size_summary_from_diff.csv",
            "Fig8CD_panel_numbers_extracted_from_heatmap_pdf.csv",
            "Fig8CD_source_note.csv",
        ],
        "raw": raw_inputs,
        "code": ["../fig8_code/fig8_source_data.py"],
        "status": "available",
        "availability_note": "Panel-level counts are directly available from summary DEG heatmap PDFs; supporting subtype diff tables are archived",
    }

    return panel_c_meta, panel_d_meta


def prepare_go_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    if "pvalue" in df.columns:
        df["pvalue"] = safe_num(df["pvalue"])
        df["neg_log10_pvalue"] = df["pvalue"].apply(neg_log10)
    return df


def build_panel_e_or_f(panel: str) -> dict:
    assert panel in {"E", "F"}
    generated: list[str] = []
    raw_inputs: list[str] = []

    if panel == "E":
        up_full, down_full = ROUND_UP_FULL, ROUND_DOWN_FULL
        up_sig, down_sig = ROUND_UP_SIG, ROUND_DOWN_SIG
        terms_down_oa = PANEL_E_DOWN_IN_OA_TERMS
        terms_up_oa = PANEL_E_UP_IN_OA_TERMS
        tag = "Fig8E"
        subpanel = "GO terms of Round_Spermatids (Down/Up in OA)"
    else:
        up_full, down_full = ELONG_UP_FULL, ELONG_DOWN_FULL
        up_sig, down_sig = ELONG_UP_SIG, ELONG_DOWN_SIG
        terms_down_oa = PANEL_F_DOWN_IN_OA_TERMS
        terms_up_oa = PANEL_F_UP_IN_OA_TERMS
        tag = "Fig8F"
        subpanel = "GO terms of Elongated_Spermatids (Down/Up in OA)"

    for src in [up_full, down_full, up_sig, down_sig]:
        rel = copy_to_raw_support(src, f"{tag}_{src.name}")
        if rel:
            raw_inputs.append(rel)

    # CtrlvsOA_up => down in OA ; CtrlvsOA_down => up in OA
    down_in_oa_df = prepare_go_table(up_full)
    up_in_oa_df = prepare_go_table(down_full)

    down_sel = pick_terms(down_in_oa_df, terms_down_oa)
    down_sel.insert(0, "direction_in_OA", "Downregulated in OA")
    down_sel.to_csv(OUT_DIR / f"{tag}_GO_selected_down_in_OA.csv", index=False)
    generated.append(f"{tag}_GO_selected_down_in_OA.csv")

    up_sel = pick_terms(up_in_oa_df, terms_up_oa)
    up_sel.insert(0, "direction_in_OA", "Upregulated in OA")
    up_sel.to_csv(OUT_DIR / f"{tag}_GO_selected_up_in_OA.csv", index=False)
    generated.append(f"{tag}_GO_selected_up_in_OA.csv")

    combined = pd.concat([down_sel, up_sel], ignore_index=True)
    combined.to_csv(OUT_DIR / f"{tag}_GO_selected_combined.csv", index=False)
    generated.append(f"{tag}_GO_selected_combined.csv")

    if not down_in_oa_df.empty:
        down_in_oa_df.to_csv(OUT_DIR / f"{tag}_GO_full_down_in_OA_source.csv", index=False)
        generated.append(f"{tag}_GO_full_down_in_OA_source.csv")
    if not up_in_oa_df.empty:
        up_in_oa_df.to_csv(OUT_DIR / f"{tag}_GO_full_up_in_OA_source.csv", index=False)
        generated.append(f"{tag}_GO_full_up_in_OA_source.csv")

    found_ratio = float(combined["found"].mean()) if not combined.empty else 0.0
    status = "available" if found_ratio >= 1.0 else "partial"
    note = (
        "All displayed GO terms were matched in source enrichment tables"
        if status == "available"
        else "Some displayed GO terms could not be matched exactly; matched nearest available entries"
    )

    return {
        "subpanel": subpanel,
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../fig8_code/fig8_source_data.py"],
        "status": status,
        "availability_note": note,
    }


def build_panel_g() -> dict:
    generated: list[str] = []
    raw_inputs: list[str] = []
    panel_g_pdf = next((p for p in PANEL_G_HEATMAP_PDF_CANDIDATES if p.exists()), PANEL_G_HEATMAP_PDF_CANDIDATES[0])

    for src in [
        PANEL_GH_ARCHIVE,
        panel_g_pdf,
        PANEL_G_DRGS_XLSX,
        PANEL_G_DRGS_CSV,
        PANEL_G_ST_MATRIX_XLSX,
        PANEL_G_ST_MATRIX_CSV,
        INFERTILITY_GENE_LIST,
        AVG_EXPR_GERM,
    ]:
        rel = copy_to_raw_support(src, f"Fig8G_{src.name}")
        if rel:
            raw_inputs.append(rel)

    genes_from_pdf, group_labels_from_pdf = extract_panel_g_genes_and_groups(panel_g_pdf)

    if genes_from_pdf:
        gl_pdf = pd.DataFrame(
            [{"gene_order": i + 1, "gene": g} for i, g in enumerate(genes_from_pdf)]
        )
        gl_pdf.to_csv(OUT_DIR / "Fig8G_heatmap_pdf_gene_order.csv", index=False)
        generated.append("Fig8G_heatmap_pdf_gene_order.csv")

    if group_labels_from_pdf:
        grp_pdf = pd.DataFrame(
            [{"group_order": i + 1, "group_label": g} for i, g in enumerate(group_labels_from_pdf)]
        )
        grp_pdf.to_csv(OUT_DIR / "Fig8G_heatmap_pdf_group_labels.csv", index=False)
        generated.append("Fig8G_heatmap_pdf_group_labels.csv")

    archive_tables = read_tar_gz_csvs(PANEL_GH_ARCHIVE)
    if archive_tables:
        for name in sorted(archive_tables):
            out_name = f"Fig8G_archive_{name}"
            archive_tables[name].to_csv(OUT_DIR / out_name, index=False)
            generated.append(out_name)

    expr_source = ""
    expr_matrix = pd.DataFrame(columns=["gene", "Ctrl_mean", "OA_mean"])
    archive_genes: list[str] = []

    hm = archive_tables.get("heatmap_export.csv", pd.DataFrame())
    if not hm.empty:
        cols = {c.lower().replace(" ", ""): c for c in hm.columns}
        gene_col = cols.get("genename")
        ctrl_col = cols.get("ctrl")
        oa_col = cols.get("oa")
        if gene_col and ctrl_col and oa_col:
            m = hm[[gene_col, ctrl_col, oa_col]].copy()
            m = m.rename(columns={gene_col: "gene", ctrl_col: "Ctrl_mean", oa_col: "OA_mean"})
            m["gene"] = m["gene"].astype(str).str.strip()
            m["Ctrl_mean"] = safe_num(m["Ctrl_mean"])
            m["OA_mean"] = safe_num(m["OA_mean"])
            m = m[m["gene"] != ""].drop_duplicates(subset=["gene"]).reset_index(drop=True)
            m.insert(0, "gene_order", range(1, len(m) + 1))
            m.to_csv(OUT_DIR / "Fig8G_heatmap_export_matrix_from_archive.csv", index=False)
            generated.append("Fig8G_heatmap_export_matrix_from_archive.csv")
            expr_source = f"{PANEL_GH_ARCHIVE}::heatmap_export.csv"
            expr_matrix = m[["gene", "Ctrl_mean", "OA_mean"]].copy()
            archive_genes = m["gene"].tolist()

    # Keep legacy candidate list output, with priority: archive heatmap > panel PDF > historical infertility list.
    gene_list = pd.DataFrame(columns=["gene"])
    if archive_genes:
        gene_list = pd.DataFrame({"gene": archive_genes})
    elif genes_from_pdf:
        gene_list = pd.DataFrame({"gene": genes_from_pdf})
    elif INFERTILITY_GENE_LIST.exists():
        gl = pd.read_csv(INFERTILITY_GENE_LIST)
        first = gl.columns[0]
        gene_list = pd.DataFrame({"gene": gl[first].astype(str).str.strip()})
        gene_list = gene_list[gene_list["gene"] != ""].drop_duplicates().reset_index(drop=True)
    gene_list.to_csv(OUT_DIR / "Fig8G_candidate_gene_list.csv", index=False)
    generated.append("Fig8G_candidate_gene_list.csv")

    # If archive matrix is not available, use previous fallback expression sources.
    if expr_matrix.empty and PANEL_G_DRGS_XLSX.exists():
        df = pd.read_excel(PANEL_G_DRGS_XLSX)
        if {"Gene id", "CTRL", "OA"}.issubset(df.columns):
            expr_source = str(PANEL_G_DRGS_XLSX)
            x = df[["Gene id", "CTRL", "OA"]].copy()
            x = x.rename(columns={"Gene id": "gene", "CTRL": "Ctrl_mean", "OA": "OA_mean"})
            x["Ctrl_mean"] = safe_num(x["Ctrl_mean"])
            x["OA_mean"] = safe_num(x["OA_mean"])
            expr_matrix = x
    if expr_matrix.empty and PANEL_G_ST_MATRIX_XLSX.exists():
        df = pd.read_excel(PANEL_G_ST_MATRIX_XLSX)
        if {"gene_name", "Ctrl", "OA"}.issubset(df.columns):
            expr_source = str(PANEL_G_ST_MATRIX_XLSX)
            x = df[["gene_name", "Ctrl", "OA"]].copy()
            x = x.rename(columns={"gene_name": "gene", "Ctrl": "Ctrl_mean", "OA": "OA_mean"})
            x["Ctrl_mean"] = safe_num(x["Ctrl_mean"])
            x["OA_mean"] = safe_num(x["OA_mean"])
            expr_matrix = x
    if expr_matrix.empty and AVG_EXPR_GERM.exists():
        df = pd.read_csv(AVG_EXPR_GERM)
        gene_col = "Gene" if "Gene" in df.columns else df.columns[0]
        x = df.rename(columns={gene_col: "gene"}).copy()
        ctrl_cols = [c for c in x.columns if c.startswith("Ctrl_")]
        oa_cols = [c for c in x.columns if c.startswith("OA_")]
        for c in ctrl_cols + oa_cols:
            x[c] = safe_num(x[c])
        x["Ctrl_mean"] = x[ctrl_cols].mean(axis=1) if ctrl_cols else float("nan")
        x["OA_mean"] = x[oa_cols].mean(axis=1) if oa_cols else float("nan")
        expr_source = str(AVG_EXPR_GERM)
        expr_matrix = x[["gene", "Ctrl_mean", "OA_mean"]].copy()

    src_registry = pd.DataFrame(
        [
            {"source_path": str(PANEL_GH_ARCHIVE), "exists": PANEL_GH_ARCHIVE.exists(), "selected": expr_source.startswith(str(PANEL_GH_ARCHIVE))},
            {"source_path": str(PANEL_G_DRGS_XLSX), "exists": PANEL_G_DRGS_XLSX.exists(), "selected": expr_source == str(PANEL_G_DRGS_XLSX)},
            {"source_path": str(PANEL_G_ST_MATRIX_XLSX), "exists": PANEL_G_ST_MATRIX_XLSX.exists(), "selected": expr_source == str(PANEL_G_ST_MATRIX_XLSX)},
            {"source_path": str(AVG_EXPR_GERM), "exists": AVG_EXPR_GERM.exists(), "selected": expr_source == str(AVG_EXPR_GERM)},
        ]
    )
    src_registry.to_csv(OUT_DIR / "Fig8G_expression_source_registry.csv", index=False)
    generated.append("Fig8G_expression_source_registry.csv")

    if not gene_list.empty and not expr_matrix.empty:
        expr_matrix["gene"] = expr_matrix["gene"].astype(str)
        expr_matrix["gene_upper"] = expr_matrix["gene"].str.upper()

        target = gene_list.copy()
        target["gene"] = target["gene"].astype(str)
        target["gene_upper"] = target["gene"].str.upper()
        target["gene_order"] = range(1, len(target) + 1)

        merged = target.merge(expr_matrix, on="gene_upper", how="left", suffixes=("_target", "_expr"))
        merged["gene"] = merged["gene_target"]
        merged = merged[["gene_order", "gene", "Ctrl_mean", "OA_mean"]]
        merged["source_table"] = expr_source
        merged.to_csv(OUT_DIR / "Fig8G_candidate_ctrl_oa_expression_matrix.csv", index=False)
        generated.append("Fig8G_candidate_ctrl_oa_expression_matrix.csv")

        cov = merged[["gene_order", "gene", "Ctrl_mean", "OA_mean"]].copy()
        cov["found_in_expression_source"] = cov["Ctrl_mean"].notna() & cov["OA_mean"].notna()
        cov.to_csv(OUT_DIR / "Fig8G_gene_coverage_in_expression_source.csv", index=False)
        generated.append("Fig8G_gene_coverage_in_expression_source.csv")

        z = merged.copy()
        means = (z["Ctrl_mean"] + z["OA_mean"]) / 2.0
        stds = (((z["Ctrl_mean"] - means) ** 2 + (z["OA_mean"] - means) ** 2) / 2.0).pow(0.5)
        z["Ctrl_z"] = (z["Ctrl_mean"] - means) / stds.replace(0, pd.NA)
        z["OA_z"] = (z["OA_mean"] - means) / stds.replace(0, pd.NA)
        z = z.fillna(0.0)
        z.to_csv(OUT_DIR / "Fig8G_candidate_ctrl_oa_expression_zscore.csv", index=False)
        generated.append("Fig8G_candidate_ctrl_oa_expression_zscore.csv")

    if genes_from_pdf and archive_genes:
        pdf_u = [x.upper() for x in genes_from_pdf]
        arc_u = [x.upper() for x in archive_genes]
        comp = pd.DataFrame(
            [
                {"source": "pdf", "n_genes": len(pdf_u), "genes": ";".join(genes_from_pdf)},
                {"source": "archive_heatmap", "n_genes": len(arc_u), "genes": ";".join(archive_genes)},
                {
                    "source": "intersection",
                    "n_genes": len(set(pdf_u) & set(arc_u)),
                    "genes": ";".join([g for g in genes_from_pdf if g.upper() in set(arc_u)]),
                },
                {
                    "source": "pdf_only",
                    "n_genes": len([g for g in genes_from_pdf if g.upper() not in set(arc_u)]),
                    "genes": ";".join([g for g in genes_from_pdf if g.upper() not in set(arc_u)]),
                },
                {
                    "source": "archive_only",
                    "n_genes": len([g for g in archive_genes if g.upper() not in set(pdf_u)]),
                    "genes": ";".join([g for g in archive_genes if g.upper() not in set(pdf_u)]),
                },
            ]
        )
        comp.to_csv(OUT_DIR / "Fig8G_archive_vs_pdf_gene_comparison.csv", index=False)
        generated.append("Fig8G_archive_vs_pdf_gene_comparison.csv")

    status = "partial"
    avail_note = "gene order and group labels are extracted from panel PDF; Ctrl/OA matrix is reconstructed from DRGs all-ST.xlsx"
    if expr_source.startswith(str(PANEL_GH_ARCHIVE)) and not expr_matrix.empty:
        status = "available"
        avail_note = "direct raw exports from OA-CTRL热图.gz are available (heatmap/dotplot/violin); group labels are taken from panel PDF text"
    elif genes_from_pdf and not expr_matrix.empty:
        status = "available"
    elif not expr_matrix.empty:
        avail_note = "Ctrl/OA matrix available, but panel-PDF gene order extraction failed"

    note = pd.DataFrame(
        [
            {
                "panel": "G",
                "status": status,
                "note": "If present, OA-CTRL热图.gz is used as direct panel G source (heatmap_export/scaled_heatmap_export/dotplot_export/scaled_dotplot_export/violin_export). Group labels are extracted from heatmap_export*.pdf text. Per-gene category assignment boundaries are still not directly encoded in these exported tables.",
                "expression_source": expr_source,
                "group_label_pdf_source": str(panel_g_pdf),
                "n_genes_from_pdf": len(genes_from_pdf),
                "n_groups_from_pdf": len(group_labels_from_pdf),
                "n_genes_from_archive_heatmap": len(archive_genes),
            }
        ]
    )
    note.to_csv(OUT_DIR / "Fig8G_source_note.csv", index=False)
    generated.append("Fig8G_source_note.csv")

    return {
        "subpanel": "Heatmap of selected infertility-related genes (grouped categories)",
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../fig8_code/fig8_source_data.py"],
        "status": status,
        "availability_note": avail_note,
    }


def build_panel_h() -> dict:
    generated: list[str] = []
    raw_inputs: list[str] = []

    rel_arch = copy_to_raw_support(PANEL_GH_ARCHIVE, "Fig8H_OA-CTRL热图.gz")
    if rel_arch:
        raw_inputs.append(rel_arch)
    rel_direct_csv = copy_to_raw_support(PANEL_H_DIRECT_CSV, "Fig8H_fig8h.csv")
    if rel_direct_csv:
        raw_inputs.append(rel_direct_csv)

    has_direct_h_csv = False
    n_direct_h_rows = 0
    n_direct_h_genes_found = 0
    required_h_cols = ["cell", "gname"] + [f"{g} normalised expression value" for g in PANEL_H_GENE_LIST]

    if PANEL_H_DIRECT_CSV.exists():
        direct = pd.read_csv(PANEL_H_DIRECT_CSV)
        n_direct_h_rows = len(direct)
        col_check = pd.DataFrame(
            [{"required_column": c, "present": c in direct.columns} for c in required_h_cols]
        )
        col_check.to_csv(OUT_DIR / "Fig8H_direct_csv_required_columns_check.csv", index=False)
        generated.append("Fig8H_direct_csv_required_columns_check.csv")
        n_direct_h_genes_found = int(col_check[col_check["required_column"].str.contains("normalised expression value")]["present"].sum())

        missing_required = [c for c in required_h_cols if c not in direct.columns]
        if not missing_required:
            id_cols = [c for c in ["cell", "gname", "UMAP1", "UMAP2", "TSNE1", "TSNE2"] if c in direct.columns]
            value_cols = [f"{g} normalised expression value" for g in PANEL_H_GENE_LIST]
            d = direct[id_cols + value_cols].copy()
            d = d.rename(columns={"gname": "group"})
            for c in value_cols:
                d[c] = safe_num(d[c])

            long_all = d.melt(
                id_vars=[c for c in id_cols if c != "gname"] + (["group"] if "group" in d.columns else []),
                value_vars=value_cols,
                var_name="feature",
                value_name="normalised_expression",
            )
            long_all["gene"] = long_all["feature"].str.replace(" normalised expression value", "", regex=False)
            long_all["group"] = pd.Categorical(long_all["group"], categories=GROUP_ORDER_ALL, ordered=True)
            long_all = long_all.sort_values(["group", "gene", "cell"]).reset_index(drop=True)
            long_all.to_csv(OUT_DIR / "Fig8H_per_cell_expression_long.csv", index=False)
            generated.append("Fig8H_per_cell_expression_long.csv")

            long_ctrl_oa = long_all[long_all["group"].isin(GROUP_ORDER_CTRL_OA)].copy()
            long_ctrl_oa["group"] = pd.Categorical(long_ctrl_oa["group"], categories=GROUP_ORDER_CTRL_OA, ordered=True)
            long_ctrl_oa = long_ctrl_oa.sort_values(["group", "gene", "cell"]).reset_index(drop=True)
            long_ctrl_oa.to_csv(OUT_DIR / "Fig8H_per_cell_expression_ctrl_oa_long.csv", index=False)
            generated.append("Fig8H_per_cell_expression_ctrl_oa_long.csv")

            summary_all = (
                long_all.groupby(["group", "gene"], observed=False)["normalised_expression"]
                .agg(
                    n_cells="count",
                    mean="mean",
                    median="median",
                    q1=lambda s: s.quantile(0.25),
                    q3=lambda s: s.quantile(0.75),
                    min="min",
                    max="max",
                )
                .reset_index()
                .sort_values(["group", "gene"])
                .reset_index(drop=True)
            )
            summary_all.to_csv(OUT_DIR / "Fig8H_per_cell_expression_group_gene_summary.csv", index=False)
            generated.append("Fig8H_per_cell_expression_group_gene_summary.csv")

            summary_ctrl_oa = (
                long_ctrl_oa.groupby(["group", "gene"], observed=False)["normalised_expression"]
                .agg(
                    n_cells="count",
                    mean="mean",
                    median="median",
                    q1=lambda s: s.quantile(0.25),
                    q3=lambda s: s.quantile(0.75),
                    min="min",
                    max="max",
                )
                .reset_index()
                .sort_values(["group", "gene"])
                .reset_index(drop=True)
            )
            summary_ctrl_oa.to_csv(OUT_DIR / "Fig8H_per_cell_expression_ctrl_oa_group_gene_summary.csv", index=False)
            generated.append("Fig8H_per_cell_expression_ctrl_oa_group_gene_summary.csv")

            has_direct_h_csv = True

    registry_rows: list[dict] = []
    for gene in PANEL_H_GENE_LIST:
        p1 = PANEL_H_VIOLIN_DIR_1 / f"{gene} normalised expression value_violin_by_gname.pdf"
        p2 = PANEL_H_VIOLIN_DIR_2 / f"{gene} normalised expression value_violin_by_gname.pdf"
        p3 = PANEL_H_VIOLIN_DIR_3 / f"基因_表达_Volin_plot_{gene}.pdf"

        candidates = [
            ("significant2_violin_by_gname", p1),
            ("violin_plots2_violin_by_gname", p2),
            ("legacy_violin_plot_pdf", p3),
        ]
        chosen = False
        for rank, (tag, src) in enumerate(candidates, start=1):
            rel = copy_to_raw_support(src, f"Fig8H_{gene}_{tag}.pdf")
            if rel:
                raw_inputs.append(rel)
            exists = src.exists()
            selected_primary = exists and (not chosen)
            if selected_primary:
                chosen = True
            registry_rows.append(
                {
                    "gene": gene,
                    "source_rank": rank,
                    "source_tag": tag,
                    "source_path": str(src),
                    "exists": exists,
                    "selected_primary": selected_primary,
                    "raw_support_copy": rel,
                }
            )

    registry = pd.DataFrame(registry_rows)
    registry.to_csv(OUT_DIR / "Fig8H_raw_violin_file_registry.csv", index=False)
    generated.append("Fig8H_raw_violin_file_registry.csv")

    archive_tables = read_tar_gz_csvs(PANEL_GH_ARCHIVE)
    archive_violin = archive_tables.get("violin_export.csv", pd.DataFrame())
    archive_presence_rows: list[dict] = []
    has_exact_archive = False
    if not archive_violin.empty:
        col_map = {c.lower(): c for c in archive_violin.columns}
        for gene in PANEL_H_GENE_LIST:
            key = f"{gene.lower()} normalised expression value"
            archive_presence_rows.append(
                {
                    "gene": gene,
                    "archive_column_expected": key,
                    "found_in_archive_violin": key in col_map,
                    "archive_column_name": col_map.get(key, ""),
                }
            )
        presence = pd.DataFrame(archive_presence_rows)
        presence.to_csv(OUT_DIR / "Fig8H_archive_gene_presence.csv", index=False)
        generated.append("Fig8H_archive_gene_presence.csv")

        if bool(presence["found_in_archive_violin"].all()):
            keep_cols = ["cell", "gname"] + [col_map[f"{g.lower()} normalised expression value"] for g in PANEL_H_GENE_LIST]
            x = archive_violin[keep_cols].copy()
            x = x.rename(columns={"gname": "group"})
            long = x.melt(id_vars=["cell", "group"], var_name="feature", value_name="normalised_expression")
            long["gene"] = long["feature"].str.replace(" normalised expression value", "", regex=False)
            long = long[["cell", "group", "gene", "normalised_expression"]].copy()
            long["normalised_expression"] = safe_num(long["normalised_expression"])
            long.to_csv(OUT_DIR / "Fig8H_archive_per_cell_expression_long.csv", index=False)
            generated.append("Fig8H_archive_per_cell_expression_long.csv")
            has_exact_archive = True
    else:
        presence = pd.DataFrame(
            [{"gene": g, "archive_column_expected": f"{g} normalised expression value", "found_in_archive_violin": False, "archive_column_name": ""} for g in PANEL_H_GENE_LIST]
        )
        presence.to_csv(OUT_DIR / "Fig8H_archive_gene_presence.csv", index=False)
        generated.append("Fig8H_archive_gene_presence.csv")

    if AVG_EXPR_GERM.exists():
        expr = pd.read_csv(AVG_EXPR_GERM)
        gene_col = "Gene" if "Gene" in expr.columns else expr.columns[0]
        expr = expr.rename(columns={gene_col: "gene"}).copy()
        keep = set(g.upper() for g in PANEL_H_GENE_LIST)
        expr = expr[expr["gene"].astype(str).str.upper().isin(keep)].copy()

        sample_cols = [c for c in expr.columns if c.startswith("Ctrl_") or c.startswith("OA_")]
        for c in sample_cols:
            expr[c] = safe_num(expr[c])

        long = expr.melt(id_vars=["gene"], value_vars=sample_cols, var_name="sample", value_name="avg_expression")
        long["group"] = long["sample"].apply(lambda s: "Ctrl" if str(s).startswith("Ctrl_") else "OA")
        long = long.sort_values(["gene", "group", "sample"]).reset_index(drop=True)
        long.to_csv(OUT_DIR / "Fig8H_candidate_sample_level_expression.csv", index=False)
        generated.append("Fig8H_candidate_sample_level_expression.csv")

        summary = summarize_score(long.rename(columns={"group": "group2", "avg_expression": "score"}), "group2", "score", GROUP_ORDER_CTRL_OA)
        summary = summary.rename(columns={"group2": "group"})
        summary.to_csv(OUT_DIR / "Fig8H_candidate_group_summary.csv", index=False)
        generated.append("Fig8H_candidate_group_summary.csv")

    n_primary = int(registry["selected_primary"].sum()) if not registry.empty else 0
    n_archive_found = int(presence["found_in_archive_violin"].sum()) if not presence.empty else 0
    panel_status = "available" if n_primary >= len(PANEL_H_GENE_LIST) else "partial"
    if has_direct_h_csv or has_exact_archive:
        panel_status = "available"

    if has_direct_h_csv:
        note_txt = "Direct per-cell table fig8h.csv is available and contains all four panel-H genes (CCDC65/DNAH8/PIH1D3/CATSPER3), with group label and UMAP/TSNE coordinates."
        avail_txt = "direct per-cell input table (fig8h.csv) and raw violin PDFs are available"
    elif has_exact_archive:
        note_txt = "Primary raw violin PDFs are available and OA-CTRL热图.gz also contains exact per-cell expression columns for all four panel-H genes."
        avail_txt = "raw violin PDFs and exact per-cell archive input are available"
    else:
        note_txt = "Primary raw violin PDFs are taken from 无精症/显著 2 (with fallback to 无精症/violin_plots2 and 5/不育基因表达). OA-CTRL热图.gz was checked but does not contain all four panel-H genes. Candidate sample-level support is exported from average_expression_matrix生精.csv."
        avail_txt = "primary raw violin PDFs are available from 显著 2; per-cell plotting table is still missing"

    note = pd.DataFrame(
        [
            {
                "panel": "H",
                "status": panel_status,
                "note": note_txt,
                "n_genes_expected": len(PANEL_H_GENE_LIST),
                "n_primary_pdfs_found": n_primary,
                "n_genes_found_in_OA_CTRL_archive_violin": n_archive_found,
                "direct_h_csv_found": PANEL_H_DIRECT_CSV.exists(),
                "n_rows_in_direct_h_csv": n_direct_h_rows,
                "n_genes_found_in_direct_h_csv": n_direct_h_genes_found,
            }
        ]
    )
    note.to_csv(OUT_DIR / "Fig8H_source_note.csv", index=False)
    generated.append("Fig8H_source_note.csv")

    return {
        "subpanel": "Four marker-gene violin plots (Ctrl vs OA)",
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../fig8_code/fig8_source_data.py"],
        "status": panel_status,
        "availability_note": avail_txt,
    }


def build_panel_i() -> dict:
    generated: list[str] = []
    raw_inputs: list[str] = []
    code_refs = ["../fig8_code/fig8_source_data.py"]

    for src in [PANEL_I_INPUT, PANEL_I_SCRIPT]:
        rel = copy_to_raw_support(src, f"Fig8I_{src.name}")
        if rel:
            raw_inputs.append(rel)
            if src == PANEL_I_SCRIPT:
                code_refs.append(rel)

    if not PANEL_I_INPUT.exists():
        note = pd.DataFrame(
            [{"panel": "I", "status": "missing", "note": "2.csv not found in 体细胞表达"}]
        )
        note.to_csv(OUT_DIR / "Fig8I_source_note.csv", index=False)
        generated.append("Fig8I_source_note.csv")
        return {
            "subpanel": "Somatic-cell composition stacked bar (Ctrl vs OA)",
            "generated": generated,
            "raw": raw_inputs,
            "code": code_refs,
            "status": "missing",
            "availability_note": "somatic composition input missing",
        }

    df = pd.read_csv(PANEL_I_INPUT).rename(columns={"Cluster": "cluster", "Sample": "group", "Number": "count"})
    df = df[df["group"].isin(GROUP_ORDER_CTRL_OA)].copy()
    df["count"] = safe_num(df["count"]).fillna(0).astype(int)
    df.to_csv(OUT_DIR / "Fig8I_somatic_composition_input.csv", index=False)
    generated.append("Fig8I_somatic_composition_input.csv")

    map_cluster = {
        "TCells": "Immune cells",
        "cDCs": "Immune cells",
        "Macrophages": "Immune cells",
        "STs": "Sertoli cells",
        "Myoids": "Myoid cells",
        "LCs": "Leydig cells",
        "ECs": "Endothelial cells",
    }
    relabeled = df.copy()
    relabeled["panel_cluster"] = relabeled["cluster"].map(map_cluster).fillna(relabeled["cluster"])
    relabeled = (
        relabeled.groupby(["group", "panel_cluster"], as_index=False)
        .agg(count=("count", "sum"))
        .sort_values(["group", "panel_cluster"])
        .reset_index(drop=True)
    )

    cluster_order = ["Immune cells", "Sertoli cells", "Myoid cells", "Leydig cells", "Endothelial cells"]
    relabeled["group"] = pd.Categorical(relabeled["group"], categories=GROUP_ORDER_CTRL_OA, ordered=True)
    relabeled["panel_cluster"] = pd.Categorical(relabeled["panel_cluster"], categories=cluster_order, ordered=True)
    relabeled = relabeled.sort_values(["group", "panel_cluster"]).reset_index(drop=True)
    relabeled.to_csv(OUT_DIR / "Fig8I_somatic_composition_relabel_long.csv", index=False)
    generated.append("Fig8I_somatic_composition_relabel_long.csv")

    relabeled["percent"] = relabeled["count"] / relabeled.groupby("group", observed=False)["count"].transform("sum")
    relabeled.to_csv(OUT_DIR / "Fig8I_somatic_composition_percent_long.csv", index=False)
    generated.append("Fig8I_somatic_composition_percent_long.csv")

    wide = relabeled.pivot(index="group", columns="panel_cluster", values="percent").reset_index()
    wide.to_csv(OUT_DIR / "Fig8I_somatic_composition_percent_wide.csv", index=False)
    generated.append("Fig8I_somatic_composition_percent_wide.csv")

    return {
        "subpanel": "Somatic-cell composition stacked bar (Ctrl vs OA)",
        "generated": generated,
        "raw": raw_inputs,
        "code": code_refs,
        "status": "available",
        "availability_note": "direct input count table and plotting script are available",
    }


def build_panel_j() -> dict:
    generated: list[str] = []
    raw_inputs: list[str] = []

    rel = copy_to_raw_support(PANEL_J_ZIP, "Fig8J_OA-CTRL.zip")
    if rel:
        raw_inputs.append(rel)

    mature = read_zip_csv(PANEL_J_ZIP, PANEL_J_MEMBERS["mature_score"])
    immature = read_zip_csv(PANEL_J_ZIP, PANEL_J_MEMBERS["immature_score"])
    expr_m = read_zip_csv(PANEL_J_ZIP, PANEL_J_MEMBERS["mature_expression"])
    expr_i = read_zip_csv(PANEL_J_ZIP, PANEL_J_MEMBERS["immature_expression"])

    ok = True

    if mature.empty or immature.empty:
        ok = False
    else:
        score_col_m = next((c for c in mature.columns if "LC-mature" in c), mature.columns[1])
        score_col_i = next((c for c in immature.columns if "LC-immature" in c), immature.columns[1])

        m = mature.rename(columns={"gname": "group", score_col_m: "score"}).copy()
        m["score"] = safe_num(m["score"])
        m = m[m["group"].isin(GROUP_ORDER_ALL)].copy()
        keep_cols = [c for c in ["cell", "group", "score", "UMAP1", "UMAP2", "TSNE1", "TSNE2"] if c in m.columns]
        m = m[keep_cols].copy()
        m["group"] = pd.Categorical(m["group"], categories=GROUP_ORDER_ALL, ordered=True)
        m = m.sort_values(["group", "cell"]).reset_index(drop=True)
        m.to_csv(OUT_DIR / "Fig8J_LC_mature_score_cells.csv", index=False)
        generated.append("Fig8J_LC_mature_score_cells.csv")

        m_sum = summarize_score(m, "group", "score", GROUP_ORDER_ALL)
        m_sum.to_csv(OUT_DIR / "Fig8J_LC_mature_score_summary.csv", index=False)
        generated.append("Fig8J_LC_mature_score_summary.csv")

        i = immature.rename(columns={"gname": "group", score_col_i: "score"}).copy()
        i["score"] = safe_num(i["score"])
        i = i[i["group"].isin(GROUP_ORDER_ALL)].copy()
        keep_cols_i = [c for c in ["cell", "group", "score", "UMAP1", "UMAP2", "TSNE1", "TSNE2"] if c in i.columns]
        i = i[keep_cols_i].copy()
        i["group"] = pd.Categorical(i["group"], categories=GROUP_ORDER_ALL, ordered=True)
        i = i.sort_values(["group", "cell"]).reset_index(drop=True)
        i.to_csv(OUT_DIR / "Fig8J_LC_immature_score_cells.csv", index=False)
        generated.append("Fig8J_LC_immature_score_cells.csv")

        i_sum = summarize_score(i, "group", "score", GROUP_ORDER_ALL)
        i_sum.to_csv(OUT_DIR / "Fig8J_LC_immature_score_summary.csv", index=False)
        generated.append("Fig8J_LC_immature_score_summary.csv")

        combined = pd.concat(
            [m.assign(score_type="LC_mature_score"), i.assign(score_type="LC_immature_score")],
            ignore_index=True,
        )
        combined.to_csv(OUT_DIR / "Fig8J_score_cells_combined_long.csv", index=False)
        generated.append("Fig8J_score_cells_combined_long.csv")

    if not expr_m.empty:
        em = expr_m.rename(columns={"gname": "group"}).copy()
        em.to_csv(OUT_DIR / "Fig8J_mature_marker_expression_cells.csv", index=False)
        generated.append("Fig8J_mature_marker_expression_cells.csv")
    if not expr_i.empty:
        ei = expr_i.rename(columns={"gname": "group"}).copy()
        ei.to_csv(OUT_DIR / "Fig8J_immature_marker_expression_cells.csv", index=False)
        generated.append("Fig8J_immature_marker_expression_cells.csv")

    note = pd.DataFrame(
        [
            {
                "panel": "J",
                "status": "available" if ok else "partial",
                "note": "LC mature/immature score tables were extracted from OA-CTRL.zip (genescore-mature.csv and genescore-immature.csv).",
            }
        ]
    )
    note.to_csv(OUT_DIR / "Fig8J_source_note.csv", index=False)
    generated.append("Fig8J_source_note.csv")

    return {
        "subpanel": "LC mature score and LC immature score violins",
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../fig8_code/fig8_source_data.py"],
        "status": "available" if ok else "partial",
        "availability_note": "score tables extracted from OA-CTRL.zip" if ok else "one or more score tables missing in OA-CTRL.zip",
    }


def build_mapping_and_availability(panel_meta: dict[str, dict]) -> None:
    rows = []
    avail = []
    order = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
    for p in order:
        meta = panel_meta[p]
        rows.append(
            {
                "panel": p,
                "subpanel": meta["subpanel"],
                "generated_source_data": ";".join(meta["generated"]),
                "raw_input": ";".join(meta["raw"]),
                "code": ";".join(meta["code"]),
                "status": meta["status"],
            }
        )
        avail.append({"panel": p, "status": meta["status"], "note": meta["availability_note"]})

    pd.DataFrame(rows).to_csv(OUT_DIR / "Fig8_file_mapping.csv", index=False)
    pd.DataFrame(avail).to_csv(OUT_DIR / "Fig8_panels_data_availability.csv", index=False)


def write_readme() -> None:
    csv_files = sorted(p.name for p in OUT_DIR.glob("*.csv"))
    text = [
        "Figure 8 source data (Panels A-J)",
        "",
        "Generated files:",
        *[f"- {name}" for name in csv_files],
        "",
        "Notes:",
        "- Panels A/B/I/J: direct count/score source tables are available.",
        "- Panels C/D: panel-level counts are directly extracted from original summary heatmap PDFs; supporting subtype diff tables are archived (full unpublished row-binning matrix is not present in current folders).",
        "- Panels E/F: selected GO term bars are mapped to full GOALL enrichment tables.",
        "- Panel G: OA-CTRL热图.gz is used as primary raw source when available (heatmap/dotplot/violin exports); group labels are supplemented by heatmap_export(1).pdf text.",
        "- Panel H: if raw_support/fig8h.csv exists, it is used as direct per-cell input; OA-CTRL热图.gz is also checked for four panel-H genes and reported in Fig8H_archive_gene_presence.csv.",
        "",
        "Build script:",
        "- ../fig8_code/fig8_source_data.py",
    ]
    (OUT_DIR / "README.txt").write_text("\n".join(text) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()

    panel_a = build_panel_a()
    panel_b = build_panel_b()
    panel_c, panel_d = build_panel_cd()
    panel_e = build_panel_e_or_f("E")
    panel_f = build_panel_e_or_f("F")
    panel_g = build_panel_g()
    panel_h = build_panel_h()
    panel_i = build_panel_i()
    panel_j = build_panel_j()

    panel_meta = {
        "A": panel_a,
        "B": panel_b,
        "C": panel_c,
        "D": panel_d,
        "E": panel_e,
        "F": panel_f,
        "G": panel_g,
        "H": panel_h,
        "I": panel_i,
        "J": panel_j,
    }

    build_mapping_and_availability(panel_meta)
    write_readme()


if __name__ == "__main__":
    main()
