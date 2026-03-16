#!/usr/bin/env python3
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parents[1]  # .../source data
ROOT = BASE.parent
OUT_DIR = BASE / "supfig6_source_data"
RAW_SUPPORT_DIR = OUT_DIR / "raw_support"

SUPFIG6_AI = ROOT / "Supplementary Figure 6.ai"
SUPFIG6_PDF = ROOT / "Supplementary Figure 6.pdf"
SUPFIG6_PNG = ROOT / "Supplementary Figure 6.png"

PANEL_A_DIR = ROOT / "原来作图的数据" / "汇总" / "61.生殖亚群的scMetabolism分析"
PANEL_A_MATRIX = PANEL_A_DIR / "P20073103.scMeta_heatmap.xls"
PANEL_A_TOP_ORDER_TXT = PANEL_A_DIR / "topn" / "P20073103.scMeta_heatmap_plot.txt"
PANEL_A_TOP_HEATMAP_PDF = PANEL_A_DIR / "topn" / "P20073103_top.scMeta_heatmap.pdf"
PANEL_A_TOP_HEATMAP_PNG = PANEL_A_DIR / "topn" / "P20073103_top.scMeta_heatmap.png"
PANEL_A_README = PANEL_A_DIR / "scMetabolism_score_README.pdf"

PANEL_B_DIR = ROOT / "原来作图的数据" / "汇总" / "15.生殖细胞的scFEA分析"
PANEL_B_PDF = PANEL_B_DIR / "P20073103_shengzhi.Metabolism_Flux_heatmap.pdf"
PANEL_B_PDF_OLD = PANEL_B_DIR / "P20073103_shengzhi.Metabolism_Flux_heatmap-25.06.17.pdf"
PANEL_B_PNG = PANEL_B_DIR / "P20073103_shengzhi.Metabolism_Flux_heatmap.png"
PANEL_B_README = PANEL_B_DIR / "scFEA_readme.pdf"

PANEL_C_DIR = ROOT / "原来作图的数据" / "汇总" / "17.生殖细胞的scFEA组间比较分析"
PANEL_C_PDF = PANEL_C_DIR / "P20073103.group_Metabolism_Flux_top10_heatmap.pdf"
PANEL_C_PDF_OLD = PANEL_C_DIR / "P20073103.group_Metabolism_Flux_top10_heatmap的副本25.06.17.pdf"
PANEL_C_PNG = PANEL_C_DIR / "P20073103.group_Metabolism_Flux_top10_heatmap.png"
PANEL_C_ALL_PDF = PANEL_C_DIR / "P20073103.group_Metabolism_Flux_heatmap_all.pdf"
PANEL_C_ALL_PNG = PANEL_C_DIR / "P20073103.group_Metabolism_Flux_heatmap_all.png"

CELLTYPE_ORDER = [
    "SSCs",
    "Early_primary_SPCs",
    "SPGs",
    "Late_primary_SPCs",
    "Round_Spermatids",
    "Elongated_Spermatids",
    "Sperm",
]
GROUP_ORDER = ["AZFc_Del", "Ctrl", "iNOA_B", "iNOA_S", "KS", "OA"]


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


def safe_read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    for sep in ["\t", ","]:
        try:
            df = pd.read_csv(path, sep=sep)
            if not df.empty and len(df.columns) >= 2:
                return df
        except Exception:
            continue
    return pd.DataFrame()


def normalize_pathway(text: str) -> str:
    s = str(text).strip()
    s = re.sub(r"\s+", " ", s)
    if s.startswith("Glycosylphosphatidylinositol") and "GPI" not in s:
        s = s.replace("Glycosylphosphatidylinositol", "Glycosylphosphatidylinositol (GPI)")
    s = s.replace("(GPI) -anchor", "(GPI)-anchor")
    return s


def build_panel_a() -> dict:
    generated: list[str] = []
    raw_inputs: list[str] = []

    for src, alias in [
        (SUPFIG6_PDF, "SupFig6_main.pdf"),
        (SUPFIG6_PNG, "SupFig6_main.png"),
        (SUPFIG6_AI, "SupFig6_main.ai"),
        (PANEL_A_MATRIX, "SupFig6A_P20073103.scMeta_heatmap.xls"),
        (PANEL_A_TOP_ORDER_TXT, "SupFig6A_P20073103.scMeta_heatmap_plot.txt"),
        (PANEL_A_TOP_HEATMAP_PDF, "SupFig6A_P20073103_top.scMeta_heatmap.pdf"),
        (PANEL_A_TOP_HEATMAP_PNG, "SupFig6A_P20073103_top.scMeta_heatmap.png"),
        (PANEL_A_README, "SupFig6A_scMetabolism_score_README.pdf"),
    ]:
        rel = copy_to_raw_support(src, alias)
        if rel:
            raw_inputs.append(rel)

    full_df = safe_read_table(PANEL_A_MATRIX)
    if full_df.empty:
        note = pd.DataFrame(
            [
                {
                    "panel": "A",
                    "status": "missing",
                    "note": "P20073103.scMeta_heatmap.xls not found; cannot export panel A numeric matrix.",
                }
            ]
        )
        note.to_csv(OUT_DIR / "SupFig6A_source_note.csv", index=False)
        generated.append("SupFig6A_source_note.csv")
        return {
            "panel": "A",
            "subpanel": "scMetabolism pathway heatmap (celltype x pathway)",
            "generated": generated,
            "raw": raw_inputs,
            "code": ["../supfig6_code/supfig6_source_data.py"],
            "status": "partial",
            "availability_note": "Panel A matrix missing in current workspace",
        }

    id_col = "id" if "id" in full_df.columns else full_df.columns[0]
    value_cols = [c for c in full_df.columns if c != id_col]

    full_df = full_df[[id_col] + value_cols].copy()
    full_df = full_df.rename(columns={id_col: "pathway"})
    full_df["pathway"] = full_df["pathway"].astype(str).map(normalize_pathway)

    for col in value_cols:
        full_df[col] = pd.to_numeric(full_df[col], errors="coerce")

    full_df.to_csv(OUT_DIR / "SupFig6A_scmeta_heatmap_full.csv", index=False)
    generated.append("SupFig6A_scmeta_heatmap_full.csv")

    top_order: list[str] = []
    top_txt_df = safe_read_table(PANEL_A_TOP_ORDER_TXT)
    if not top_txt_df.empty:
        top_col = top_txt_df.columns[0]
        for value in top_txt_df[top_col].tolist():
            pathway = normalize_pathway(value)
            if pathway and pathway not in top_order:
                top_order.append(pathway)

    # Preferred: the topn text file already contains the exact panel-A values.
    top30_df = pd.DataFrame()
    if not top_txt_df.empty:
        top_id_col = top_txt_df.columns[0]
        top_value_cols: list[str] = []
        for col in top_txt_df.columns[1:]:
            numeric = pd.to_numeric(top_txt_df[col], errors="coerce")
            if int(numeric.notna().sum()) >= max(3, int(0.5 * len(top_txt_df))):
                top_value_cols.append(col)

        if top_value_cols:
            top30_df = top_txt_df[[top_id_col] + top_value_cols].copy()
            top30_df = top30_df.rename(columns={top_id_col: "pathway"})
            top30_df["pathway"] = top30_df["pathway"].astype(str).map(normalize_pathway)
            for col in top_value_cols:
                top30_df[col] = pd.to_numeric(top30_df[col], errors="coerce")
            if top_order:
                top30_df["pathway"] = pd.Categorical(top30_df["pathway"], categories=top_order, ordered=True)
                top30_df = top30_df.sort_values("pathway").reset_index(drop=True)
                top30_df["pathway"] = top30_df["pathway"].astype(str)

    # Fallback: subset full matrix by top order.
    if top30_df.empty:
        top30_df = full_df.copy()
        if top_order:
            top30_df = top30_df[top30_df["pathway"].isin(top_order)].copy()
            top30_df["pathway"] = pd.Categorical(top30_df["pathway"], categories=top_order, ordered=True)
            top30_df = top30_df.sort_values("pathway").reset_index(drop=True)
            top30_df["pathway"] = top30_df["pathway"].astype(str)
        else:
            top30_df = top30_df.head(30).copy()

    top30_df.to_csv(OUT_DIR / "SupFig6A_scmeta_heatmap_top30.csv", index=False)
    generated.append("SupFig6A_scmeta_heatmap_top30.csv")

    top30_long = top30_df.melt(id_vars=["pathway"], var_name="celltype", value_name="expression_zscore")
    pathway_rank = {p: i + 1 for i, p in enumerate(top30_df["pathway"].tolist())}
    top30_long["pathway_rank"] = top30_long["pathway"].map(pathway_rank)
    celltype_rank = {c: i + 1 for i, c in enumerate(CELLTYPE_ORDER)}
    top30_long["celltype_rank"] = top30_long["celltype"].map(celltype_rank)
    top30_long = top30_long.sort_values(["pathway_rank", "celltype_rank", "celltype"]).reset_index(drop=True)
    top30_long.to_csv(OUT_DIR / "SupFig6A_scmeta_heatmap_top30_long.csv", index=False)
    generated.append("SupFig6A_scmeta_heatmap_top30_long.csv")

    pd.DataFrame(
        {
            "pathway_rank": range(1, len(top30_df) + 1),
            "pathway": top30_df["pathway"].tolist(),
        }
    ).to_csv(OUT_DIR / "SupFig6A_pathway_order.csv", index=False)
    generated.append("SupFig6A_pathway_order.csv")

    summary = (
        top30_long.groupby("celltype", as_index=False)
        .agg(
            n_pathways=("pathway", "count"),
            mean_expression=("expression_zscore", "mean"),
            median_expression=("expression_zscore", "median"),
            min_expression=("expression_zscore", "min"),
            max_expression=("expression_zscore", "max"),
        )
        .copy()
    )
    summary["celltype"] = pd.Categorical(summary["celltype"], categories=CELLTYPE_ORDER, ordered=True)
    summary = summary.sort_values("celltype").reset_index(drop=True)
    summary.to_csv(OUT_DIR / "SupFig6A_celltype_summary.csv", index=False)
    generated.append("SupFig6A_celltype_summary.csv")

    note = pd.DataFrame(
        [
            {
                "panel": "A",
                "status": "available",
                "note": "Top pathway values are sourced from P20073103.scMeta_heatmap_plot.txt; full heatmap matrix is from P20073103.scMeta_heatmap.xls.",
            }
        ]
    )
    note.to_csv(OUT_DIR / "SupFig6A_source_note.csv", index=False)
    generated.append("SupFig6A_source_note.csv")

    return {
        "panel": "A",
        "subpanel": "scMetabolism pathway heatmap (celltype x pathway)",
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../supfig6_code/supfig6_source_data.py", "../supfig6_code/supfig6_panelA_heatmap.py"],
        "status": "available",
        "availability_note": "Numeric matrix and top30 order are available",
    }


def pdftotext_lines(pdf_path: Path) -> list[str]:
    if not pdf_path.exists():
        return []
    try:
        result = subprocess.run(
            ["pdftotext", str(pdf_path), "-"],
            capture_output=True,
            text=True,
            check=True,
        )
    except Exception:
        return []

    lines: list[str] = []
    for raw in result.stdout.splitlines():
        text = re.sub(r"\s+", " ", raw.strip())
        if text:
            lines.append(text)
    return lines


def extract_reaction_labels(lines: list[str]) -> list[str]:
    out: list[str] = []
    for line in lines:
        if ">" not in line:
            continue
        if re.fullmatch(r"[-0-9. ]+", line):
            continue
        label = line.replace("−", "-")
        label = re.sub(r"\s*->\s*", " -> ", label)
        label = re.sub(r"\s+", " ", label).strip()
        # pdftotext may drop the leading "(Glc)3 (Glc" fragment for this long reaction.
        if label.startswith("NAc)2 (Man)9 (PP-Dol)1+Protein asparagine ->"):
            label = "(Glc)3 (GlcNAc)2 (Man)9 (PP-Dol)1+Protein asparagine -> (Glc)3 (GlcNAc)2 (Man)9 (Asn)1"
        if label and label not in out:
            out.append(label)
    return out


def extract_ordered_labels(lines: list[str], allowed: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for line in lines:
        for item in allowed:
            if line == item and item not in seen:
                out.append(item)
                seen.add(item)
    return out


def build_panel_b_or_c(panel: str) -> dict:
    assert panel in {"B", "C"}

    generated: list[str] = []
    raw_inputs: list[str] = []

    if panel == "B":
        src_files = [
            (PANEL_B_PDF, "SupFig6B_P20073103_shengzhi.Metabolism_Flux_heatmap.pdf"),
            (PANEL_B_PDF_OLD, "SupFig6B_P20073103_shengzhi.Metabolism_Flux_heatmap-25.06.17.pdf"),
            (PANEL_B_PNG, "SupFig6B_P20073103_shengzhi.Metabolism_Flux_heatmap.png"),
            (PANEL_B_README, "SupFig6B_scFEA_readme.pdf"),
        ]
        pdf_for_extract = PANEL_B_PDF
        subpanel = "scFEA flux heatmap across germ cell subtypes"
    else:
        src_files = [
            (PANEL_C_PDF, "SupFig6C_P20073103.group_Metabolism_Flux_top10_heatmap.pdf"),
            (PANEL_C_PDF_OLD, "SupFig6C_P20073103.group_Metabolism_Flux_top10_heatmap-25.06.17-copy.pdf"),
            (PANEL_C_PNG, "SupFig6C_P20073103.group_Metabolism_Flux_top10_heatmap.png"),
            (PANEL_C_ALL_PDF, "SupFig6C_P20073103.group_Metabolism_Flux_heatmap_all.pdf"),
            (PANEL_C_ALL_PNG, "SupFig6C_P20073103.group_Metabolism_Flux_heatmap_all.png"),
        ]
        pdf_for_extract = PANEL_C_PDF
        subpanel = "scFEA group-comparison flux heatmap (celltype + group)"

    for src, alias in src_files:
        rel = copy_to_raw_support(src, alias)
        if rel:
            raw_inputs.append(rel)

    lines = pdftotext_lines(pdf_for_extract)
    reaction_labels = extract_reaction_labels(lines)
    celltypes = extract_ordered_labels(lines, CELLTYPE_ORDER)
    groups = extract_ordered_labels(lines, GROUP_ORDER)

    reaction_df = pd.DataFrame(
        {
            "reaction_rank": range(1, len(reaction_labels) + 1),
            "reaction_label": reaction_labels,
        }
    )
    reaction_file = f"SupFig6{panel}_reaction_axis_labels.csv"
    reaction_df.to_csv(OUT_DIR / reaction_file, index=False)
    generated.append(reaction_file)

    celltype_df = pd.DataFrame(
        {
            "celltype_rank": range(1, len(celltypes) + 1),
            "celltype": celltypes,
        }
    )
    celltype_file = f"SupFig6{panel}_celltype_axis_labels.csv"
    celltype_df.to_csv(OUT_DIR / celltype_file, index=False)
    generated.append(celltype_file)

    if panel == "C":
        group_df = pd.DataFrame(
            {
                "group_rank": range(1, len(groups) + 1),
                "group": groups,
            }
        )
        group_file = "SupFig6C_group_axis_labels.csv"
        group_df.to_csv(OUT_DIR / group_file, index=False)
        generated.append(group_file)

    note_text = (
        f"Panel {panel} numeric flux matrix (reaction x samples) was not found in current workspace; "
        "only final heatmap exports and axis labels extracted from PDF are available."
    )
    note_df = pd.DataFrame(
        [
            {
                "panel": panel,
                "status": "partial",
                "n_reactions_from_pdf": len(reaction_labels),
                "n_celltypes_from_pdf": len(celltypes),
                "n_groups_from_pdf": len(groups),
                "note": note_text,
            }
        ]
    )
    note_file = f"SupFig6{panel}_source_note.csv"
    note_df.to_csv(OUT_DIR / note_file, index=False)
    generated.append(note_file)

    return {
        "panel": panel,
        "subpanel": subpanel,
        "generated": generated,
        "raw": raw_inputs,
        "code": ["../supfig6_code/supfig6_source_data.py"],
        "status": "partial",
        "availability_note": "Only heatmap exports and axis labels are available; raw numeric matrix not located",
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

    pd.DataFrame(mapping_rows).to_csv(OUT_DIR / "SupFig6_file_mapping.csv", index=False)
    pd.DataFrame(availability_rows).to_csv(OUT_DIR / "SupFig6_panels_data_availability.csv", index=False)


def write_readme(rows: list[dict]) -> None:
    lines = [
        "Supplementary Figure 6 source data (Panels A-C)",
        "",
        "Generated files:",
    ]

    files = []
    for row in rows:
        files.extend(row["generated"])
    files.extend(["SupFig6_file_mapping.csv", "SupFig6_panels_data_availability.csv"])

    for name in sorted(set(files)):
        lines.append(f"- {name}")

    lines.extend(
        [
            "",
            "Notes:",
            "- Panel A: direct numeric matrix and top30 pathway order are available.",
            "- Panels B/C: current workspace has final heatmap PDFs/PNGs, but no raw flux numeric matrix table.",
            "- For full reproducibility of B/C heatmaps, add original reaction-by-sample flux matrix files (csv/xls/tsv).",
            "",
            "Build script:",
            "- ../supfig6_code/supfig6_source_data.py",
        ]
    )

    (OUT_DIR / "README.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()

    rows = [
        build_panel_a(),
        build_panel_b_or_c("B"),
        build_panel_b_or_c("C"),
    ]

    write_mapping(rows)
    write_readme(rows)

    print("Done: generated SupFig6 source-data package in", OUT_DIR)


if __name__ == "__main__":
    main()
