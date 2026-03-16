Supplementary Figure 3 source data (Panels A-B)

Files:
- SupFig3A_Late_primary_SPCs_GSEA_all.csv: all Hallmark rows from 4 comparisons (OA/AZFc_Del/KS/iNOA_B vs Ctrl).
- SupFig3B_Round_Spermatids_GSEA_all.csv: all Hallmark rows from 3 comparisons (OA/AZFc_Del/KS vs Ctrl).
- SupFig3A_plot_matrix.csv: panel A bubble-matrix values after pathway selection (max_pathways=24).
- SupFig3B_plot_matrix.csv: panel B bubble-matrix values after pathway selection (max_pathways=28).
- raw_support/*.xls: original GSEA enrichment tables used to build the above CSV files.

Selection logic follows:
- pathway ranking by minimum adjusted p-value, then by maximum absolute NES.
- point size = -log10(FDR), capped at 20.