from __future__ import annotations

import tarfile
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
RAW_C = BASE / "fig6c.gz"
RAW_D = BASE / "fig6d.gz"
OUT_DIR = BASE / "fig6_source_data"
RAW_SUPPORT_C = OUT_DIR / "raw_support_fig6c"
RAW_SUPPORT_D = OUT_DIR / "raw_support_fig6d"

GENE_ORDER = ["GAS1", "CDON", "BOC", "HHIP", "DHH"]
GROUP_ORDER_NO_OA = ["Ctrl", "iNOA_B", "iNOA_S", "AZFc_Del", "KS"]
GROUP_ORDER_WITH_OA = ["Ctrl", "iNOA_B", "iNOA_S", "OA", "AZFc_Del", "KS"]


def extract_tar_gz(src: Path, dst: Path) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    with tarfile.open(src, "r:gz") as tf:
        tf.extractall(dst)


def order_dotplot(df: pd.DataFrame, groups: list[str]) -> pd.DataFrame:
    out = df.copy()
    out = out[out["gname"].isin(groups) & out["Gene name"].isin(GENE_ORDER)]
    out["gname"] = pd.Categorical(out["gname"], categories=groups, ordered=True)
    out["Gene name"] = pd.Categorical(out["Gene name"], categories=GENE_ORDER, ordered=True)
    return out.sort_values(["Gene name", "gname"]).reset_index(drop=True)


def order_heatmap(df: pd.DataFrame, groups: list[str]) -> pd.DataFrame:
    keep_cols = ["Gene name"] + [g for g in groups if g in df.columns]
    out = df.copy()
    out = out[out["Gene name"].isin(GENE_ORDER)][keep_cols]
    out["Gene name"] = pd.Categorical(out["Gene name"], categories=GENE_ORDER, ordered=True)
    return out.sort_values("Gene name").reset_index(drop=True)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    extract_tar_gz(RAW_C, RAW_SUPPORT_C)
    extract_tar_gz(RAW_D, RAW_SUPPORT_D)

    # Panel C (LC_a)
    c_dot = pd.read_csv(RAW_SUPPORT_C / "dotplot_export.csv")
    c_heat = pd.read_csv(RAW_SUPPORT_C / "heatmap_export.csv")
    c_violin = pd.read_csv(RAW_SUPPORT_C / "violin_export.csv")

    c_dot_ordered = order_dotplot(c_dot, GROUP_ORDER_NO_OA)
    c_heat_ordered = order_heatmap(c_heat, GROUP_ORDER_NO_OA)

    c_dot_ordered.to_csv(OUT_DIR / "Fig6C_dotplot_LCa.csv", index=False)
    c_heat_ordered.to_csv(OUT_DIR / "Fig6C_heatmap_LCa.csv", index=False)
    c_violin.to_csv(OUT_DIR / "Fig6C_violin_cells_LCa.csv", index=False)

    # Panel D (ST)
    d_dot = pd.read_csv(RAW_SUPPORT_D / "dotplot_export.csv")
    d_heat = pd.read_csv(RAW_SUPPORT_D / "heatmap_export.csv")
    d_violin = pd.read_csv(RAW_SUPPORT_D / "violin_export.csv")

    d_dot_full = order_dotplot(d_dot, GROUP_ORDER_WITH_OA)
    d_dot_no_oa = order_dotplot(d_dot, GROUP_ORDER_NO_OA)
    d_heat_full = order_heatmap(d_heat, GROUP_ORDER_WITH_OA)
    d_heat_no_oa = order_heatmap(d_heat, GROUP_ORDER_NO_OA)

    d_dot_full.to_csv(OUT_DIR / "Fig6D_dotplot_ST_full_with_OA.csv", index=False)
    d_dot_no_oa.to_csv(OUT_DIR / "Fig6D_dotplot_ST_no_OA.csv", index=False)
    d_heat_full.to_csv(OUT_DIR / "Fig6D_heatmap_ST_full_with_OA.csv", index=False)
    d_heat_no_oa.to_csv(OUT_DIR / "Fig6D_heatmap_ST_no_OA.csv", index=False)
    d_violin.to_csv(OUT_DIR / "Fig6D_violin_cells_ST.csv", index=False)

    # quick metadata note
    meta = pd.DataFrame(
        [
            {
                "panel": "C",
                "raw_archive": str(RAW_C.relative_to(BASE)),
                "raw_dotplot": str((RAW_SUPPORT_C / "dotplot_export.csv").relative_to(BASE)),
                "raw_heatmap": str((RAW_SUPPORT_C / "heatmap_export.csv").relative_to(BASE)),
                "groups_used": ",".join(GROUP_ORDER_NO_OA),
                "genes": ",".join(GENE_ORDER),
            },
            {
                "panel": "D",
                "raw_archive": str(RAW_D.relative_to(BASE)),
                "raw_dotplot": str((RAW_SUPPORT_D / "dotplot_export.csv").relative_to(BASE)),
                "raw_heatmap": str((RAW_SUPPORT_D / "heatmap_export.csv").relative_to(BASE)),
                "groups_used": ",".join(GROUP_ORDER_NO_OA),
                "groups_raw_has_oa": True,
                "genes": ",".join(GENE_ORDER),
            },
        ]
    )
    meta.to_csv(OUT_DIR / "Fig6CD_source_note.csv", index=False)


if __name__ == "__main__":
    main()
