Supplementary Figure 4 source data (Panels A-F)

Generated files:
- SupFig4A_ST_cells_pseudotime.csv
- SupFig4A_stage_markers_from_diffgene.csv
- SupFig4B_regulon_stage_matrix_target.csv
- SupFig4B_regulon_stage_matrix_detected.csv
- SupFig4C_gene_set_candidates.csv
- SupFig4C_glycolysis_score_cells.csv
- SupFig4C_glycolysis_score_summary.csv
- SupFig4D_ST1_disease_vs_ctrl_gsea_all.csv
- SupFig4D_ST1_disease_vs_ctrl_top10.csv
- SupFig4E_gene_set_candidates.csv
- SupFig4E_inflammatory_score_cells.csv
- SupFig4E_inflammatory_score_summary.csv
- SupFig4F_ST3_disease_vs_ctrl_gsea_all.csv
- SupFig4F_ST3_disease_vs_ctrl_top10.csv
- SupFig4C_E_missing_cell_score_raw.csv
- SupFig4_file_mapping.csv
- SupFig4_panels_data_availability.csv
- SourceData_SupFig4.xlsx

Raw support copied from panel_raw:
- 21 files in ./raw_support/

Notes:
- Panel B target regulons follow Supplementary Figure 4 label set and are joined to ST stage via monocle cell IDs.
- Panel C glycolysis score table has been located and exported at cell level.
- Panel E inflammatory/immune score table has been located and exported at cell level.
- Panel D/F now use original ST disease-vs-control KEGG GSEA tables (ST1-diseasevsctrl.csv / ST3-diseasevsctrl.csv).
- Legacy ST1vsST3/ ST3vsST1 summary CSVs are removed to avoid version confusion.
- Plotting scripts for reconstructed/proxy views are in ../supfig4_code/.