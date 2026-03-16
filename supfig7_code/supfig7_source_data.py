#!/usr/bin/env python3
from __future__ import annotations

import math
import re
import shutil
import subprocess
import tarfile
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parents[1]  # .../source data
ROOT = BASE.parent
OUT_DIR = BASE / "supfig7_source_data"
RAW_SUPPORT_DIR = OUT_DIR / "raw_support"

# Target panel image used in this round (legacy SupFig7 A-G layout)
TARGET_SUPFIG7_IMAGE = (
    ROOT / "原来作图的数据" / "无精症" / "20251210郭立强" / "Supplementary Figure 7.png"
)

# Panel A/B: somatic DEG overlap heatmaps (OAvsCtrl)
PANEL_A_UP_PDF = ROOT / "原来作图的数据" / "5" / "P20073103_OAvsCtrl_up_DEG_heatmap.pdf"
PANEL_B_DOWN_PDF = ROOT / "原来作图的数据" / "5" / "P20073103_OAvsCtrl_down_DEG_heatmap.pdf"
DIFF_BASE = ROOT / "原来作图的数据" / "汇总" / "44.大类亚群组间差异基因及功能富集分析" / "Diff"

SUBTYPE_ORDER = ["ECs", "LCs", "Lym", "Myoids", "STs"]
SUBTYPE_TO_FOLDER = {
    "ECs": "ECs",
    "LCs": "LCs",
    "Lym": "Lym",
    "Myoids": "Myoid",
    "STs": "STs",
}
# Numeric order in original OAvsCtrl heatmap PDFs
PDF_SUBTYPE_ORDER = ["STs", "Myoids", "Lym", "LCs", "ECs"]

# Panel C/D: GO bars
LC_GO_UP_IN_OA = ROOT / "原来作图的数据" / "5" / "LC" / "OAvsCTRL-LC" / "GO.csv"
LC_GO_DOWN_IN_OA = ROOT / "原来作图的数据" / "5" / "LC" / "CTRLvsOA-LC" / "GO.csv"
LC_GO_SUMMARY = ROOT / "原来作图的数据" / "5" / "LC" / "LC-GO.csv"

ST_GO_UP_IN_OA = ROOT / "原来作图的数据" / "5" / "ST" / "OAvsCRTL-ST" / "GO.csv"
ST_GO_DOWN_IN_OA = ROOT / "原来作图的数据" / "5" / "ST" / "CTRLvsOA-ST" / "GO.csv"
ST_GO_SUMMARY = ROOT / "原来作图的数据" / "5" / "ST" / "ST-GO.csv"

PANEL_C_DOWN_TERMS = [
    "Protein folding",
    "Response to oxidative stress",
    "Cellular oxidant detoxification",
    "RNA splicing",
    "Response to interleukin-7",
]
PANEL_C_UP_TERMS = [
    "Aerobic respiration",
    "Cellular respiration",
    "Electron transport chain",
    "Respiratory electron transport chain",
    "Oxidative phosphorylation",
]

PANEL_D_DOWN_TERMS = [
    "Supramolecular fiber organization",
    "Actin filament organization",
    "Actin cytoskeleton organization",
    "Actin filament-based process",
    "Cell-substrate adhesion",
]
PANEL_D_UP_TERMS = [
    "Cellular respiration",
    "ATP metabolic process",
    "Aerobic respiration",
    "Oxidative phosphorylation",
    "Aerobic electron transport chain",
]

# Panel E: cilia-gene dotplot
PANEL_E_PDF = ROOT / "原来作图的数据" / "无精症" / "纤毛组装相关基因OAvsCTRL气泡图.pdf"
PANEL_E_TABLE = ROOT / "原来作图的数据" / "无精症" / "Late SPC-其他vs control富集.csv"
PANEL_E_COMPARE_TAR_GZ_CANDIDATES = [
    Path("/Users/xq/Desktop/compare_by_expression.tar.gz"),
    Path("/Users/xq/Desktop/supig7e.gz"),
]
PANEL_E_GENE_ORDER = [
    "WDR63",
    "TTC25",
    "TTC29",
    "CATSPER3",
    "PIH1D3",
    "RSPH3",
    "RSPH1",
    "CFAP61",
    "CFAP44",
    "CCDC65",
    "DNAH6",
    "DNAH7",
    "DNAH8",
]
PANEL_E_GENE_PDF_DIRS = [
    ROOT / "原来作图的数据" / "无精症" / "显著 2",
    ROOT / "原来作图的数据" / "无精症" / "显著",
]

# Panel F/G: regulon RSS and heatmap
PANEL_F_OA_PDF = ROOT / "原来作图的数据" / "无精症" / "Regulon_特异性散点图_-_top_5_regulons_OA.pdf"
PANEL_F_CTRL_PDF = ROOT / "原来作图的数据" / "无精症" / "Regulon_特异性散点图_-_top_5_regulons_Ctrl.pdf"

PANEL_G_HEATMAP_CANDIDATES = [
    ROOT / "原来作图的数据" / "无精症" / "所选分组中的Regulon平均表达热图.pdf",
    ROOT / "原来作图的数据" / "4" / "所选分组中的Regulon平均表达热图.pdf",
    ROOT / "原来作图的数据" / "3" / "所选分组中的Regulon平均表达热图.pdf",
    ROOT / "原来作图的数据" / "3" / "所选分组中的Regulon平均表达热图 (1).pdf",
]

GROUP_LABEL_CANDIDATES = ["Ctrl", "OA", "AZFc_Del", "iNOA_B", "iNOA_S", "KS", "ST1", "ST2", "ST3"]


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


def safe_read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    for sep in ["\t", ","]:
        try:
            df = pd.read_csv(path, sep=sep)
            if df.shape[1] >= 2:
                return df
        except Exception:
            continue
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


def build_overlap_matrix(gene_sets: dict[str, set[str]], order: list[str]) -> pd.DataFrame:
    rows: list[dict] = []
    for a in order:
        row = {"subtype": a}
        for b in order:
            row[b] = len(gene_sets.get(a, set()) & gene_sets.get(b, set()))
        rows.append(row)
    return pd.DataFrame(rows)


def prepare_go_table(path: Path) -> pd.DataFrame:
    df = safe_read_table(path)
    if df.empty:
        return df
    if "pvalue" in df.columns:
        df["pvalue"] = safe_num(df["pvalue"])
        df["neg_log10_pvalue"] = df["pvalue"].apply(neg_log10)
    return df


def pick_terms(df: pd.DataFrame, terms: list[str]) -> pd.DataFrame:
    if df.empty or "Description" not in df.columns:
        return pd.DataFrame(
            columns=[
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
            ]
        )

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
        if "pvalue" in row.columns:
            row["pvalue"] = safe_num(row["pvalue"])
            best = row.sort_values("pvalue", ascending=True).iloc[0]
        else:
            best = row.iloc[0]

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

    return pd.DataFrame(rows).sort_values("term_order").reset_index(drop=True)


def parse_regulon_token(token: str) -> tuple[str, int | None]:
    m = re.fullmatch(r"([A-Za-z0-9_-]+)\((\d+)g\)", token)
    if not m:
        return token, None
    return m.group(1), int(m.group(2))


def extract_regulon_tokens(lines: list[str]) -> list[str]:
    tokens: list[str] = []
    for line in lines:
        for tok in re.findall(r"[A-Za-z0-9_-]+\(\d+g\)", line):
            if tok not in tokens:
                tokens.append(tok)
    return tokens


def build_panel_ab() -> tuple[dict, dict]:
    generated: list[str] = []
    raw_inputs: list[str] = []

    for src, alias in [
        (TARGET_SUPFIG7_IMAGE, "SupFig7_target_layout_legacy.png"),
        (PANEL_A_UP_PDF, "SupFig7A_P20073103_OAvsCtrl_up_DEG_heatmap.pdf"),
        (PANEL_B_DOWN_PDF, "SupFig7B_P20073103_OAvsCtrl_down_DEG_heatmap.pdf"),
    ]:
        rel = copy_to_raw_support(src, alias)
        if rel:
            raw_inputs.append(rel)

    diff_paths: dict[str, Path] = {}
    reg_rows = []
    for subtype in SUBTYPE_ORDER:
        folder = SUBTYPE_TO_FOLDER[subtype]
        fp = DIFF_BASE / folder / "P20073103_OAvsCtrl.diffexpressed.xls"
        diff_paths[subtype] = fp
        rel = copy_to_raw_support(fp, f"SupFig7AB_{subtype}_OAvsCtrl.diffexpressed.xls")
        if rel:
            raw_inputs.append(rel)
        reg_rows.append({"subtype": subtype, "file": str(fp), "exists": fp.exists(), "raw_support_copy": rel})

    pd.DataFrame(reg_rows).to_csv(OUT_DIR / "SupFig7AB_oa_vs_ctrl_diff_file_registry.csv", index=False)
    generated.append("SupFig7AB_oa_vs_ctrl_diff_file_registry.csv")

    up_sets: dict[str, set[str]] = {}
    down_sets: dict[str, set[str]] = {}
    up_rows: list[pd.DataFrame] = []
    down_rows: list[pd.DataFrame] = []

    for subtype in SUBTYPE_ORDER:
        fp = diff_paths[subtype]
        if not fp.exists():
            up_sets[subtype] = set()
            down_sets[subtype] = set()
            continue

        df = pd.read_csv(fp, sep="\t")
        for col in ["avg_logFC", "p_val", "p_val_adj", "OA", "Ctrl", "pct.1", "pct.2"]:
            if col in df.columns:
                df[col] = safe_num(df[col])

        # OAvsCtrl: avg_logFC > 0 => up in OA; < 0 => down in OA
        up = df[(df["p_val_adj"] < 0.05) & (df["avg_logFC"] > 0)].copy()
        down = df[(df["p_val_adj"] < 0.05) & (df["avg_logFC"] < 0)].copy()

        up_sets[subtype] = set(up["gene_id"].astype(str))
        down_sets[subtype] = set(down["gene_id"].astype(str))

        keep_cols = [
            c
            for c in ["gene_id", "avg_logFC", "p_val", "p_val_adj", "OA", "Ctrl", "pct.1", "pct.2"]
            if c in df.columns
        ]
        if keep_cols:
            up2 = up[keep_cols].copy()
            up2.insert(0, "subtype", subtype)
            up_rows.append(up2)

            down2 = down[keep_cols].copy()
            down2.insert(0, "subtype", subtype)
            down_rows.append(down2)

    if up_rows:
        up_df = pd.concat(up_rows, ignore_index=True).sort_values(["subtype", "gene_id"]).reset_index(drop=True)
    else:
        up_df = pd.DataFrame(columns=["subtype", "gene_id"])
    up_df.to_csv(OUT_DIR / "SupFig7A_up_in_OA_genes_from_diff_long.csv", index=False)
    generated.append("SupFig7A_up_in_OA_genes_from_diff_long.csv")

    if down_rows:
        down_df = pd.concat(down_rows, ignore_index=True).sort_values(["subtype", "gene_id"]).reset_index(drop=True)
    else:
        down_df = pd.DataFrame(columns=["subtype", "gene_id"])
    down_df.to_csv(OUT_DIR / "SupFig7B_down_in_OA_genes_from_diff_long.csv", index=False)
    generated.append("SupFig7B_down_in_OA_genes_from_diff_long.csv")

    overlap_up = build_overlap_matrix(up_sets, SUBTYPE_ORDER)
    overlap_up.to_csv(OUT_DIR / "SupFig7AB_overlap_matrix_up_in_OA.csv", index=False)
    generated.append("SupFig7AB_overlap_matrix_up_in_OA.csv")

    overlap_down = build_overlap_matrix(down_sets, SUBTYPE_ORDER)
    overlap_down.to_csv(OUT_DIR / "SupFig7AB_overlap_matrix_down_in_OA.csv", index=False)
    generated.append("SupFig7AB_overlap_matrix_down_in_OA.csv")

    summary_rows = []
    shared_up = set.intersection(*[up_sets[s] for s in SUBTYPE_ORDER]) if SUBTYPE_ORDER else set()
    shared_down = set.intersection(*[down_sets[s] for s in SUBTYPE_ORDER]) if SUBTYPE_ORDER else set()

    for subtype in SUBTYPE_ORDER:
        other_up = set().union(*[up_sets[s] for s in SUBTYPE_ORDER if s != subtype])
        other_down = set().union(*[down_sets[s] for s in SUBTYPE_ORDER if s != subtype])
        summary_rows.append(
            {
                "subtype": subtype,
                "n_up_in_OA": len(up_sets[subtype]),
                "n_up_in_OA_unique": len(up_sets[subtype] - other_up),
                "n_down_in_OA": len(down_sets[subtype]),
                "n_down_in_OA_unique": len(down_sets[subtype] - other_down),
            }
        )

    summary_rows.append(
        {
            "subtype": "Shared_all_5_subtypes",
            "n_up_in_OA": len(shared_up),
            "n_up_in_OA_unique": None,
            "n_down_in_OA": len(shared_down),
            "n_down_in_OA_unique": None,
        }
    )

    pd.DataFrame(summary_rows).to_csv(OUT_DIR / "SupFig7AB_set_size_summary_from_diff.csv", index=False)
    generated.append("SupFig7AB_set_size_summary_from_diff.csv")

    up_nums = extract_heatmap_numbers(PANEL_A_UP_PDF)
    down_nums = extract_heatmap_numbers(PANEL_B_DOWN_PDF)

    def panel_num_rows(panel: str, nums: list[int], meaning: str) -> tuple[list[dict], dict]:
        shared = nums[0] if len(nums) >= 1 else None
        spec = nums[1:6] if len(nums) >= 6 else []
        spec_map = {k: None for k in SUBTYPE_ORDER}
        for k, v in zip(PDF_SUBTYPE_ORDER, spec):
            spec_map[k] = v

        rows = []
        for subtype in SUBTYPE_ORDER:
            rows.append(
                {
                    "panel": panel,
                    "meaning": meaning,
                    "shared_count": shared,
                    "subtype": subtype,
                    "subtype_specific_count": spec_map.get(subtype),
                    "pdf_specific_order": PDF_SUBTYPE_ORDER.index(subtype) + 1,
                    "panel_x_order": SUBTYPE_ORDER.index(subtype) + 1,
                    "numbers_raw": ",".join(str(x) for x in nums),
                }
            )

        wide = {
            "panel": panel,
            "meaning": meaning,
            "shared_count": shared,
            "numbers_raw": ",".join(str(x) for x in nums),
        }
        for subtype in SUBTYPE_ORDER:
            wide[f"{subtype}_specific_count"] = spec_map.get(subtype)
        return rows, wide

    long_rows: list[dict] = []
    wide_rows: list[dict] = []

    r1, w1 = panel_num_rows("A", up_nums, "Upregulated in OA")
    r2, w2 = panel_num_rows("B", down_nums, "Downregulated in OA")
    long_rows.extend(r1)
    long_rows.extend(r2)
    wide_rows.extend([w1, w2])

    pd.DataFrame(long_rows).to_csv(OUT_DIR / "SupFig7AB_panel_numbers_extracted_from_heatmap_pdf.csv", index=False)
    generated.append("SupFig7AB_panel_numbers_extracted_from_heatmap_pdf.csv")

    pd.DataFrame(wide_rows).to_csv(OUT_DIR / "SupFig7AB_panel_numbers_wide.csv", index=False)
    generated.append("SupFig7AB_panel_numbers_wide.csv")

    note = pd.DataFrame(
        [
            {
                "panel": "A/B",
                "pdf_panelA_shared": up_nums[0] if len(up_nums) >= 1 else None,
                "pdf_panelB_shared": down_nums[0] if len(down_nums) >= 1 else None,
                "derived_shared_up_from_diff": len(shared_up),
                "derived_shared_down_from_diff": len(shared_down),
                "note": "Panel A/B displayed counts are extracted directly from original OAvsCtrl summary heatmap PDFs; subtype diff tables are provided as support and may differ due to upstream binning/filter details.",
            }
        ]
    )
    note.to_csv(OUT_DIR / "SupFig7AB_source_note.csv", index=False)
    generated.append("SupFig7AB_source_note.csv")

    panel_a_meta = {
        "panel": "A",
        "subpanel": "Number of DEGs (upregulated in OA) across somatic subtypes",
        "generated": [
            "SupFig7AB_oa_vs_ctrl_diff_file_registry.csv",
            "SupFig7A_up_in_OA_genes_from_diff_long.csv",
            "SupFig7AB_overlap_matrix_up_in_OA.csv",
            "SupFig7AB_set_size_summary_from_diff.csv",
            "SupFig7AB_panel_numbers_extracted_from_heatmap_pdf.csv",
            "SupFig7AB_panel_numbers_wide.csv",
            "SupFig7AB_source_note.csv",
        ],
        "raw": raw_inputs,
        "code": ["../supfig7_code/supfig7_source_data.py"],
        "status": "available",
        "availability_note": "Panel-level counts are directly available from original OAvsCtrl heatmap PDF; subtype diff tables are archived as support.",
    }

    panel_b_meta = {
        "panel": "B",
        "subpanel": "Number of DEGs (downregulated in OA) across somatic subtypes",
        "generated": [
            "SupFig7AB_oa_vs_ctrl_diff_file_registry.csv",
            "SupFig7B_down_in_OA_genes_from_diff_long.csv",
            "SupFig7AB_overlap_matrix_down_in_OA.csv",
            "SupFig7AB_set_size_summary_from_diff.csv",
            "SupFig7AB_panel_numbers_extracted_from_heatmap_pdf.csv",
            "SupFig7AB_panel_numbers_wide.csv",
            "SupFig7AB_source_note.csv",
        ],
        "raw": raw_inputs,
        "code": ["../supfig7_code/supfig7_source_data.py"],
        "status": "available",
        "availability_note": "Panel-level counts are directly available from original OAvsCtrl heatmap PDF; subtype diff tables are archived as support.",
    }

    return panel_a_meta, panel_b_meta


def build_panel_c_or_d(panel: str) -> dict:
    assert panel in {"C", "D"}

    generated: list[str] = []
    raw_inputs: list[str] = []

    if panel == "C":
        down_src = LC_GO_DOWN_IN_OA
        up_src = LC_GO_UP_IN_OA
        summary_src = LC_GO_SUMMARY
        down_terms = PANEL_C_DOWN_TERMS
        up_terms = PANEL_C_UP_TERMS
        tag = "SupFig7C"
        subpanel = "GO of Leydig (Down/Up in OA)"
    else:
        down_src = ST_GO_DOWN_IN_OA
        up_src = ST_GO_UP_IN_OA
        summary_src = ST_GO_SUMMARY
        down_terms = PANEL_D_DOWN_TERMS
        up_terms = PANEL_D_UP_TERMS
        tag = "SupFig7D"
        subpanel = "GO of Sertoli (Down/Up in OA)"

    for src, alias in [
        (down_src, f"{tag}_{down_src.parent.name}_GO.csv"),
        (up_src, f"{tag}_{up_src.parent.name}_GO.csv"),
        (summary_src, f"{tag}_{summary_src.name}"),
    ]:
        rel = copy_to_raw_support(src, alias)
        if rel:
            raw_inputs.append(rel)

    down_df = prepare_go_table(down_src)
    up_df = prepare_go_table(up_src)

    down_sel = pick_terms(down_df, down_terms)
    down_sel.insert(0, "direction_in_OA", "Downregulated in OA")
    down_sel.to_csv(OUT_DIR / f"{tag}_GO_selected_down_in_OA.csv", index=False)
    generated.append(f"{tag}_GO_selected_down_in_OA.csv")

    up_sel = pick_terms(up_df, up_terms)
    up_sel.insert(0, "direction_in_OA", "Upregulated in OA")
    up_sel.to_csv(OUT_DIR / f"{tag}_GO_selected_up_in_OA.csv", index=False)
    generated.append(f"{tag}_GO_selected_up_in_OA.csv")

    combined = pd.concat([down_sel, up_sel], ignore_index=True)
    combined.to_csv(OUT_DIR / f"{tag}_GO_selected_combined.csv", index=False)
    generated.append(f"{tag}_GO_selected_combined.csv")

    if not down_df.empty:
        down_df.to_csv(OUT_DIR / f"{tag}_GO_full_down_in_OA_source.csv", index=False)
        generated.append(f"{tag}_GO_full_down_in_OA_source.csv")
    if not up_df.empty:
        up_df.to_csv(OUT_DIR / f"{tag}_GO_full_up_in_OA_source.csv", index=False)
        generated.append(f"{tag}_GO_full_up_in_OA_source.csv")

    found_ratio = float(combined["found"].mean()) if not combined.empty else 0.0
    status = "available" if found_ratio >= 1.0 else "partial"
    note = (
        "All displayed GO terms were matched in source enrichment tables"
        if status == "available"
        else "Some displayed GO terms were not matched exactly in current GO table"
    )

    return {
        "panel": panel,
        "subpanel": subpanel,
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../supfig7_code/supfig7_source_data.py"],
        "status": status,
        "availability_note": note,
    }


def build_panel_e() -> dict:
    generated: list[str] = []
    raw_inputs: list[str] = []

    for src, alias in [
        (PANEL_E_PDF, "SupFig7E_cilia_gene_dotplot.pdf"),
        (PANEL_E_TABLE, "SupFig7E_Late_SPC_others_vs_control_full_table.csv"),
    ]:
        rel = copy_to_raw_support(src, alias)
        if rel:
            raw_inputs.append(rel)

    compare_tar = next((p for p in PANEL_E_COMPARE_TAR_GZ_CANDIDATES if p.exists()), None)
    if compare_tar is not None:
        rel = copy_to_raw_support(compare_tar, "SupFig7E_compare_by_expression.tar.gz")
        if rel:
            raw_inputs.append(rel)

    gene_pdf_rows = []
    for gene in PANEL_E_GENE_ORDER:
        found_path = None
        for d in PANEL_E_GENE_PDF_DIRS:
            fp = d / f"{gene} normalised expression value_violin_by_gname.pdf"
            if fp.exists():
                found_path = fp
                break

        copied = ""
        if found_path is not None:
            copied = copy_to_raw_support(found_path, f"SupFig7E_{gene}_violin_by_gname.pdf")
            if copied:
                raw_inputs.append(copied)

        gene_pdf_rows.append(
            {
                "gene": gene,
                "source_pdf": str(found_path) if found_path is not None else "",
                "exists": bool(found_path is not None),
                "raw_support_copy": copied,
            }
        )

    pd.DataFrame(gene_pdf_rows).to_csv(OUT_DIR / "SupFig7E_gene_violin_pdf_registry.csv", index=False)
    generated.append("SupFig7E_gene_violin_pdf_registry.csv")

    # Preferred source: compare_by_expression.tar.gz (Ctrl/OA dotplot/heatmap direct exports).
    if compare_tar is not None:
        tables = read_tar_gz_csvs(compare_tar)
        for name, df in sorted(tables.items()):
            out_name = f"SupFig7E_compare_{name}"
            df.to_csv(OUT_DIR / out_name, index=False)
            generated.append(out_name)

        dot = tables.get("dotplot_export.csv", pd.DataFrame())
        heat = tables.get("heatmap_export.csv", pd.DataFrame())
        violin = tables.get("violin_export.csv", pd.DataFrame())

        if not dot.empty and {"gname", "Gene name", "Normalised expression value", "Percentage"}.issubset(dot.columns):
            d = dot.copy()
            d["gname"] = d["gname"].astype(str).str.strip()
            d["gname"] = d["gname"].str.replace("control", "Ctrl", case=False, regex=False)
            d["gname"] = d["gname"].str.replace("ctrl", "Ctrl", case=False, regex=False)
            d["gname"] = d["gname"].str.replace("oa", "OA", case=False, regex=False)
            d = d[d["gname"].isin(["Ctrl", "OA"])].copy()
            d["gene"] = d["Gene name"].astype(str).str.strip()
            d["gene_upper"] = d["gene"].str.upper()
            d["normalised_expression_value"] = safe_num(d["Normalised expression value"])
            d["percentage"] = safe_num(d["Percentage"])

            wanted_rows = []
            for i, gene in enumerate(PANEL_E_GENE_ORDER, start=1):
                for grp in ["Ctrl", "OA"]:
                    wanted_rows.append(
                        {
                            "gene_order": i,
                            "gene": gene,
                            "gene_upper": gene.upper(),
                            "gname": grp,
                        }
                    )
            wanted = pd.DataFrame(wanted_rows)

            merged = wanted.merge(
                d[["gname", "gene", "gene_upper", "normalised_expression_value", "percentage"]],
                on=["gene_upper", "gname"],
                how="left",
                suffixes=("_target", ""),
            )
            merged["gene"] = merged["gene_target"]
            merged = merged.drop(columns=["gene_target"]).sort_values(["gene_order", "gname"]).reset_index(drop=True)
            merged["found_in_dotplot_export"] = merged["normalised_expression_value"].notna() & merged["percentage"].notna()

            merged[
                [
                    "gene_order",
                    "gene",
                    "gname",
                    "normalised_expression_value",
                    "percentage",
                    "found_in_dotplot_export",
                ]
            ].to_csv(OUT_DIR / "SupFig7E_cilia13_dotplot_ctrl_oa_long.csv", index=False)
            generated.append("SupFig7E_cilia13_dotplot_ctrl_oa_long.csv")

            expr_wide = (
                merged.pivot(index=["gene_order", "gene"], columns="gname", values="normalised_expression_value")
                .reset_index()
                .rename_axis(None, axis=1)
            )
            expr_wide.to_csv(OUT_DIR / "SupFig7E_cilia13_expression_ctrl_oa_wide.csv", index=False)
            generated.append("SupFig7E_cilia13_expression_ctrl_oa_wide.csv")

            pct_wide = (
                merged.pivot(index=["gene_order", "gene"], columns="gname", values="percentage")
                .reset_index()
                .rename_axis(None, axis=1)
            )
            pct_wide.to_csv(OUT_DIR / "SupFig7E_cilia13_percentage_ctrl_oa_wide.csv", index=False)
            generated.append("SupFig7E_cilia13_percentage_ctrl_oa_wide.csv")

            if not heat.empty and {"Gene name", "Ctrl", "OA"}.issubset(heat.columns):
                h = heat.copy()
                h["gene"] = h["Gene name"].astype(str).str.strip()
                h["gene_upper"] = h["gene"].str.upper()
                h["Ctrl"] = safe_num(h["Ctrl"])
                h["OA"] = safe_num(h["OA"])
                h2 = (
                    pd.DataFrame(
                        {
                            "gene_order": range(1, len(PANEL_E_GENE_ORDER) + 1),
                            "gene": PANEL_E_GENE_ORDER,
                            "gene_upper": [g.upper() for g in PANEL_E_GENE_ORDER],
                        }
                    )
                    .merge(h[["gene_upper", "Ctrl", "OA"]], on="gene_upper", how="left")
                    .drop(columns=["gene_upper"])
                )
                h2.to_csv(OUT_DIR / "SupFig7E_cilia13_heatmap_ctrl_oa_wide.csv", index=False)
                generated.append("SupFig7E_cilia13_heatmap_ctrl_oa_wide.csv")

            if not violin.empty and {"cell", "gname"}.issubset(violin.columns):
                keep_cols = ["cell", "gname"] + [
                    c for c in violin.columns if c.endswith(" normalised expression value") and c.split(" normalised expression value")[0] in PANEL_E_GENE_ORDER
                ]
                v = violin[keep_cols].copy()
                v.to_csv(OUT_DIR / "SupFig7E_cilia13_violin_cells_from_compare_export.csv", index=False)
                generated.append("SupFig7E_cilia13_violin_cells_from_compare_export.csv")

            n_found_genes = int(merged.groupby("gene")["found_in_dotplot_export"].max().sum())
            n_found_pairs = int(merged["found_in_dotplot_export"].sum())
            missing_genes = merged.groupby("gene")["found_in_dotplot_export"].max()
            missing_list = [g for g, ok in missing_genes.items() if not ok]

            check = pd.DataFrame(
                [
                    {
                        "source": str(compare_tar),
                        "n_target_genes": len(PANEL_E_GENE_ORDER),
                        "n_found_genes": n_found_genes,
                        "n_found_gene_group_pairs": n_found_pairs,
                        "n_expected_gene_group_pairs": len(PANEL_E_GENE_ORDER) * 2,
                        "missing_genes": ";".join(missing_list),
                    }
                ]
            )
            check.to_csv(OUT_DIR / "SupFig7E_gene_coverage_check.csv", index=False)
            generated.append("SupFig7E_gene_coverage_check.csv")

            fully_available = n_found_genes == len(PANEL_E_GENE_ORDER) and n_found_pairs == len(PANEL_E_GENE_ORDER) * 2
            status = "available" if fully_available else "partial"
            availability_note = (
                "Direct Ctrl/OA expression and percentage tables are available from compare_by_expression.tar.gz."
                if fully_available
                else "compare_by_expression.tar.gz is available, but some gene/group values are missing."
            )

            note = pd.DataFrame(
                [
                    {
                        "panel": "E",
                        "status": status,
                        "source": str(compare_tar),
                        "note": availability_note,
                    }
                ]
            )
            note.to_csv(OUT_DIR / "SupFig7E_source_note.csv", index=False)
            generated.append("SupFig7E_source_note.csv")

            return {
                "panel": "E",
                "subpanel": "Cilia-assembly genes dotplot (Ctrl vs OA)",
                "generated": generated,
                "raw": raw_inputs,
                "code": ["../supfig7_code/supfig7_source_data.py"],
                "status": status,
                "availability_note": availability_note,
            }

    if not PANEL_E_TABLE.exists():
        note = pd.DataFrame(
            [
                {
                    "panel": "E",
                    "status": "missing",
                    "note": "Late SPC-其他vs control富集.csv not found in current workspace.",
                }
            ]
        )
        note.to_csv(OUT_DIR / "SupFig7E_source_note.csv", index=False)
        generated.append("SupFig7E_source_note.csv")
        return {
            "panel": "E",
            "subpanel": "Cilia-assembly genes dotplot (Ctrl vs OA)",
            "generated": generated,
            "raw": raw_inputs,
            "code": ["../supfig7_code/supfig7_source_data.py"],
            "status": "missing",
            "availability_note": "Panel E input table missing",
        }

    df = pd.read_csv(PANEL_E_TABLE)
    gene_col = "Gene id" if "Gene id" in df.columns else next((c for c in df.columns if "gene" in c.lower()), df.columns[0])

    x = df.copy()
    x["gene_upper"] = x[gene_col].astype(str).str.upper()

    wanted = pd.DataFrame(
        {
            "gene_order": range(1, len(PANEL_E_GENE_ORDER) + 1),
            "gene": PANEL_E_GENE_ORDER,
            "gene_upper": [g.upper() for g in PANEL_E_GENE_ORDER],
        }
    )

    merged = wanted.merge(x, on="gene_upper", how="left")
    merged["found_in_source"] = merged[gene_col].notna()

    rename_map = {
        "Average expression (normalized)": "average_expression_normalized",
        "Log2FC": "log2FC",
        "Pct1": "pct1",
        "Pct2": "pct2",
        "P_val": "p_value",
        "Adjusted p value": "adjusted_p_value",
        "Z score": "z_score",
    }
    for old, new in rename_map.items():
        if old in merged.columns:
            merged[new] = merged[old]

    out_cols = [
        "gene_order",
        "gene",
        "found_in_source",
        "average_expression_normalized",
        "log2FC",
        "pct1",
        "pct2",
        "z_score",
        "p_value",
        "adjusted_p_value",
    ]
    for c in out_cols:
        if c not in merged.columns:
            merged[c] = None

    merged[out_cols].to_csv(OUT_DIR / "SupFig7E_cilia13_table_from_LateSPC_vs_control.csv", index=False)
    generated.append("SupFig7E_cilia13_table_from_LateSPC_vs_control.csv")

    pct_long_rows = []
    for _, row in merged.iterrows():
        pct_long_rows.append(
            {
                "gene_order": row["gene_order"],
                "gene": row["gene"],
                "group_inferred": "Ctrl",
                "percent": row.get("pct1", None),
                "inference": "from Pct1",
            }
        )
        pct_long_rows.append(
            {
                "gene_order": row["gene_order"],
                "gene": row["gene"],
                "group_inferred": "OA",
                "percent": row.get("pct2", None),
                "inference": "from Pct2",
            }
        )

    pd.DataFrame(pct_long_rows).to_csv(OUT_DIR / "SupFig7E_cilia13_percent_by_group_inferred.csv", index=False)
    generated.append("SupFig7E_cilia13_percent_by_group_inferred.csv")

    check = pd.DataFrame(
        [
            {
                "n_target_genes": len(PANEL_E_GENE_ORDER),
                "n_found_in_table": int(merged["found_in_source"].sum()),
                "missing_genes": ";".join(merged.loc[~merged["found_in_source"], "gene"].astype(str).tolist()),
                "note": "Fallback source (Late SPC-其他vs control富集.csv): Pct1/Pct2 are available; explicit two-group mean expression columns (Ctrl/OA) are not directly present.",
            }
        ]
    )
    check.to_csv(OUT_DIR / "SupFig7E_gene_coverage_check.csv", index=False)
    generated.append("SupFig7E_gene_coverage_check.csv")

    note = pd.DataFrame(
        [
            {
                "panel": "E",
                "status": "partial",
                "note": "Fallback source used: gene list is present in Late SPC-其他vs control富集.csv, but direct Ctrl/OA mean-expression inputs are missing.",
            }
        ]
    )
    note.to_csv(OUT_DIR / "SupFig7E_source_note.csv", index=False)
    generated.append("SupFig7E_source_note.csv")

    return {
        "panel": "E",
        "subpanel": "Cilia-assembly genes dotplot (Ctrl vs OA)",
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../supfig7_code/supfig7_source_data.py"],
        "status": "partial",
        "availability_note": "Gene-level source table is available, but direct Ctrl/OA mean-expression inputs for bubble color are not explicitly provided.",
    }


def build_panel_f() -> tuple[dict, list[str]]:
    generated: list[str] = []
    raw_inputs: list[str] = []

    for src, alias in [
        (PANEL_F_OA_PDF, "SupFig7F_top5_regulons_OA.pdf"),
        (PANEL_F_CTRL_PDF, "SupFig7F_top5_regulons_Ctrl.pdf"),
    ]:
        rel = copy_to_raw_support(src, alias)
        if rel:
            raw_inputs.append(rel)

    oa_lines = pdftotext_lines(PANEL_F_OA_PDF)
    ctrl_lines = pdftotext_lines(PANEL_F_CTRL_PDF)

    oa_tokens = extract_regulon_tokens(oa_lines)
    ctrl_tokens = extract_regulon_tokens(ctrl_lines)

    def token_df(tokens: list[str], group: str) -> pd.DataFrame:
        rows = []
        for i, tok in enumerate(tokens, start=1):
            regulon, n_target = parse_regulon_token(tok)
            rows.append(
                {
                    "group": group,
                    "rank_in_pdf_text": i,
                    "token": tok,
                    "regulon": regulon,
                    "n_target_genes": n_target,
                }
            )
        return pd.DataFrame(rows)

    oa_df = token_df(oa_tokens, "OA")
    ctrl_df = token_df(ctrl_tokens, "Ctrl")

    oa_df.to_csv(OUT_DIR / "SupFig7F_top5_regulons_OA_from_pdf_text.csv", index=False)
    generated.append("SupFig7F_top5_regulons_OA_from_pdf_text.csv")

    ctrl_df.to_csv(OUT_DIR / "SupFig7F_top5_regulons_Ctrl_from_pdf_text.csv", index=False)
    generated.append("SupFig7F_top5_regulons_Ctrl_from_pdf_text.csv")

    combined = pd.concat([oa_df, ctrl_df], ignore_index=True)
    combined.to_csv(OUT_DIR / "SupFig7F_top5_regulons_combined_from_pdf_text.csv", index=False)
    generated.append("SupFig7F_top5_regulons_combined_from_pdf_text.csv")

    oa_set = set(oa_tokens)
    ctrl_set = set(ctrl_tokens)
    union = list(dict.fromkeys(oa_tokens + ctrl_tokens))

    union_rows = []
    for tok in union:
        regulon, n_target = parse_regulon_token(tok)
        union_rows.append(
            {
                "token": tok,
                "regulon": regulon,
                "n_target_genes": n_target,
                "in_OA_top5": tok in oa_set,
                "in_Ctrl_top5": tok in ctrl_set,
                "in_both": (tok in oa_set) and (tok in ctrl_set),
            }
        )

    union_df = pd.DataFrame(union_rows)
    union_df.to_csv(OUT_DIR / "SupFig7F_union_regulons_from_OA_Ctrl_top5.csv", index=False)
    generated.append("SupFig7F_union_regulons_from_OA_Ctrl_top5.csv")

    tick_rows = []
    for label, lines in [("OA", oa_lines), ("Ctrl", ctrl_lines)]:
        for line in lines:
            if re.fullmatch(r"\d+\.\d+", line):
                tick_rows.append({"group": label, "rss_axis_tick": float(line)})

    if tick_rows:
        pd.DataFrame(tick_rows).drop_duplicates().to_csv(
            OUT_DIR / "SupFig7F_rss_axis_ticks_from_pdf_text.csv", index=False
        )
        generated.append("SupFig7F_rss_axis_ticks_from_pdf_text.csv")

    note = pd.DataFrame(
        [
            {
                "panel": "F",
                "n_top_tokens_oa": len(oa_tokens),
                "n_top_tokens_ctrl": len(ctrl_tokens),
                "n_union": len(union),
                "note": "Top regulon labels are extracted from OA/Ctrl RSS PDFs. Full regulon-by-RSS numeric vector is not available in current workspace.",
            }
        ]
    )
    note.to_csv(OUT_DIR / "SupFig7F_source_note.csv", index=False)
    generated.append("SupFig7F_source_note.csv")

    return {
        "panel": "F",
        "subpanel": "Regulon specificity score curve (OA top5 highlighted)",
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../supfig7_code/supfig7_source_data.py"],
        "status": "partial",
        "availability_note": "Top regulon labels are available from PDF text; full RSS numeric vectors are not located.",
    }, union


def build_panel_g(union_tokens_from_f: list[str]) -> dict:
    generated: list[str] = []
    raw_inputs: list[str] = []

    chosen_pdf = next((p for p in PANEL_G_HEATMAP_CANDIDATES if p.exists()), PANEL_G_HEATMAP_CANDIDATES[0])

    for src in PANEL_G_HEATMAP_CANDIDATES:
        rel = copy_to_raw_support(src, f"SupFig7G_candidate_{src.parent.name}_{src.name}")
        if rel:
            raw_inputs.append(rel)

    lines = pdftotext_lines(chosen_pdf)

    group_hits: list[str] = []
    for g in GROUP_LABEL_CANDIDATES:
        for line in lines:
            if re.fullmatch(re.escape(g), line):
                if g not in group_hits:
                    group_hits.append(g)

    # Also scan inline tokens in each line.
    flat_tokens: list[str] = []
    for line in lines:
        parts = re.split(r"\s+", line)
        for tok in parts:
            tok = tok.strip()
            if tok:
                flat_tokens.append(tok)

    for g in GROUP_LABEL_CANDIDATES:
        if g in flat_tokens and g not in group_hits:
            group_hits.append(g)

    skip = set(GROUP_LABEL_CANDIDATES + ["Matrix", "2", "1", "0", "-1", "-2"])
    regulon_hits: list[str] = []
    for tok in flat_tokens:
        if tok in skip:
            continue
        if re.fullmatch(r"[A-Za-z0-9_-]+", tok):
            if tok not in regulon_hits:
                regulon_hits.append(tok)

    if regulon_hits:
        pd.DataFrame(
            [{"order": i + 1, "label": x} for i, x in enumerate(regulon_hits)]
        ).to_csv(OUT_DIR / "SupFig7G_regulon_labels_from_pdf_text.csv", index=False)
        generated.append("SupFig7G_regulon_labels_from_pdf_text.csv")

    if group_hits:
        pd.DataFrame(
            [{"order": i + 1, "group": x} for i, x in enumerate(group_hits)]
        ).to_csv(OUT_DIR / "SupFig7G_group_labels_from_pdf_text.csv", index=False)
        generated.append("SupFig7G_group_labels_from_pdf_text.csv")

    if union_tokens_from_f:
        rows = []
        for i, tok in enumerate(union_tokens_from_f, start=1):
            reg, n_target = parse_regulon_token(tok)
            rows.append(
                {
                    "order": i,
                    "token": tok,
                    "regulon": reg,
                    "n_target_genes": n_target,
                }
            )
        exp_df = pd.DataFrame(rows)
        exp_df.to_csv(OUT_DIR / "SupFig7G_expected_regulon_union_from_panelF.csv", index=False)
        generated.append("SupFig7G_expected_regulon_union_from_panelF.csv")

        template_rows = []
        for reg in exp_df["regulon"].tolist():
            for grp in ["Ctrl", "OA"]:
                template_rows.append(
                    {
                        "regulon": reg,
                        "group": grp,
                        "average_expression_value": None,
                        "value_source": "to_fill_if_raw_matrix_available",
                    }
                )
        pd.DataFrame(template_rows).to_csv(OUT_DIR / "SupFig7G_template_regulon_matrix_ctrl_oa.csv", index=False)
        generated.append("SupFig7G_template_regulon_matrix_ctrl_oa.csv")

    note = pd.DataFrame(
        [
            {
                "panel": "G",
                "selected_pdf": str(chosen_pdf),
                "n_regulon_labels_detected_from_text": len(regulon_hits),
                "n_group_labels_detected_from_text": len(group_hits),
                "note": "Current workspace contains candidate regulon heatmap PDFs, but direct 8x2 numeric matrix used in final panel G is not found as csv/xls; extracted labels and template are provided.",
            }
        ]
    )
    note.to_csv(OUT_DIR / "SupFig7G_source_note.csv", index=False)
    generated.append("SupFig7G_source_note.csv")

    return {
        "panel": "G",
        "subpanel": "Regulon average-expression heatmap (Ctrl vs OA)",
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../supfig7_code/supfig7_source_data.py"],
        "status": "partial",
        "availability_note": "Candidate heatmap exports are located, but exact final 8x2 numeric matrix is not available as a direct table.",
    }


def write_mapping(rows: list[dict]) -> None:
    mapping_rows = []
    availability_rows = []

    for row in rows:
        mapping_rows.append(
            {
                "panel": row["panel"],
                "subpanel": row["subpanel"],
                "generated_source_data": ";".join(row["generated"]),
                "raw_input": ";".join(row["raw"]),
                "code": ";".join(row["code"]),
                "status": row["status"],
            }
        )
        availability_rows.append(
            {
                "panel": row["panel"],
                "status": row["status"],
                "note": row["availability_note"],
            }
        )

    pd.DataFrame(mapping_rows).to_csv(OUT_DIR / "SupFig7_file_mapping.csv", index=False)
    pd.DataFrame(availability_rows).to_csv(OUT_DIR / "SupFig7_panels_data_availability.csv", index=False)


def write_readme(rows: list[dict]) -> None:
    lines = [
        "Supplementary Figure 7 source data (legacy A-G layout)",
        "",
        "Generated files:",
    ]

    files = []
    for row in rows:
        files.extend(row["generated"])
    files.extend(["SupFig7_file_mapping.csv", "SupFig7_panels_data_availability.csv"])

    for name in sorted(set(files)):
        lines.append(f"- {name}")

    lines.extend(
        [
            "",
            "Notes:",
            "- Panels A/B: displayed DEG counts are directly extracted from original OAvsCtrl summary heatmap PDFs; subtype diff tables are included as support.",
            "- Panels C/D: selected GO terms are matched directly to LC/ST GO.csv enrichment tables.",
            "- Panel E: priority source is compare_by_expression.tar.gz (Ctrl/OA dotplot/heatmap/violin exports). If missing, fallback is Late SPC-其他vs control富集.csv.",
            "- Panel F: top regulon labels are extracted from OA/Ctrl RSS PDFs; full RSS vectors are not located.",
            "- Panel G: candidate regulon heatmap PDFs are archived; exact final 8x2 numeric matrix is still missing in current workspace.",
            "",
            "Build script:",
            "- ../supfig7_code/supfig7_source_data.py",
        ]
    )

    (OUT_DIR / "README.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()

    panel_a, panel_b = build_panel_ab()
    panel_c = build_panel_c_or_d("C")
    panel_d = build_panel_c_or_d("D")
    panel_e = build_panel_e()
    panel_f, union_tokens = build_panel_f()
    panel_g = build_panel_g(union_tokens)

    rows = [panel_a, panel_b, panel_c, panel_d, panel_e, panel_f, panel_g]

    write_mapping(rows)
    write_readme(rows)

    print("Done: generated SupFig7 source-data package in", OUT_DIR)


if __name__ == "__main__":
    main()
