Fig1 Source Data generated from available exports.

Created files:
- Fig1B_UMAP_by_group.csv (from expression.csv)
- Fig1D_marker_feature_long.csv (from expression.csv)
- cell_annotation_from_clusters.csv (from expression (1).csv + cluster->cell type mapping)
- Fig1C_UMAP_by_celltype.csv (from cell annotations, if UMAP available)
- Fig1G_dotplot_full_genes.csv (cell type dotplot from dot.csv, if available)
- Fig1G_dotplot_available_genes.csv (fallback: cell type dotplot for available genes only)
- Fig1E_group_composition_* (from view_cellular_composition_of_samples.zip)
- Fig1F_sample_somatic_germ_counts.csv (from view_cellular_composition_of_samples (1)/count.csv, with age/count override from count.xlsx if present)
- Fig1F_countxlsx_mapping.csv (sample-level mapping pulled from count.xlsx, if present)
- compare_dotplot_export.csv, compare_scaled_dotplot_export.csv (sample-level; not cell-type-level)
- cluster_to_celltype_template.csv (manual fill if needed)

Missing data (not available in exports):
- Donor ages remain NA only if count.xlsx is unavailable.

Notes:
- Cluster-to-cell-type mapping comes from screenshot iShot_2026-02-12_10.56.46.png.
- Fig1G_dotplot_full_genes.csv uses dot.csv (marker list export). GAPDH is excluded as a housekeeping gene.
- Fig1G_dotplot_available_genes.csv only includes genes present in expression.csv.
