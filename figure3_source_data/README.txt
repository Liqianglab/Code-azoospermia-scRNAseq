Figure 3 source data (Panels A–F)

This folder contains normalized source-data tables built from local raw-support files in:
- raw_support/

Generated tables
- Fig3A_somatic_umap_cells.csv
- Fig3A_somatic_cell_counts_by_group.csv
- Fig3A_somatic_cell_ratio_by_group.csv
- Fig3B_radar_matrix_group_by_celltype.csv
- Fig3B_radar_matrix_celltype_by_group.csv
- Fig3B_radar_matrix_long.csv
- Fig3C_venn_gene_sets_long.csv
- Fig3C_venn_group_gene_counts.csv
- Fig3C_venn_LC_wide.csv
- Fig3C_venn_ST_wide.csv
- Fig3C_venn_EC_wide.csv
- Fig3C_venn_Myoid_wide.csv
- Fig3C_venn_Immune_wide.csv
- Fig3C_bubble_pathways_long.csv
- Fig3C_bubble_top10_pathways_by_group.csv
- Fig3D_cytokine_score_cells.csv
- Fig3D_cytokine_score_group_summary.csv
- Fig3E_sasp_score_cells.csv
- Fig3E_sasp_score_group_summary.csv
- Fig3F_heatmap_matrix_wide.csv
- Fig3F_heatmap_matrix_long.csv

Notes
- Panel A uses `cell_annotation_from_clusters.csv` (with UMAP + cell_type + group) and merges `Lym` + `Myeloid` as `Immune cells`.
- Panel B uses `demo (2).xlsx` matrix directly.
- Panel C Venn uses 5 separate `Demo_Venn_5.xlsx` files (LC/ST/EC/Myoid/Immune); bubble uses `table_PKD1_down_*.csv`.
- Panel D and Panel E use per-cell pathway score files and preserve original `gname` groups.
- Panel F uses `sasp-gene.zip` -> `heatmap_export/heatmap_export.xlsx`.

Build script
- ../figure3_code/figure3_source_data.py

Run
- python3 ../figure3_code/figure3_source_data.py
