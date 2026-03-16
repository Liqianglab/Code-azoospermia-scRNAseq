Figure 6 source data (Panels A-I)

Generated files:
- Fig6A_network_edges_manual.csv
- Fig6AB_hedgehog_interaction_matrix_long.csv
- Fig6AB_hedgehog_interaction_summary.csv
- Fig6C_dotplot_LCa.csv
- Fig6C_heatmap_LCa.csv
- Fig6C_violin_cells_LCa.csv
- Fig6D_dotplot_ST_full_with_OA.csv
- Fig6D_dotplot_ST_no_OA.csv
- Fig6D_heatmap_ST_full_with_OA.csv
- Fig6D_heatmap_ST_no_OA.csv
- Fig6D_violin_cells_ST.csv
- Fig6CD_source_note.csv
- Fig6E_volcano_group1_all.csv
- Fig6E_volcano_group1_threshold_summary.csv
- Fig6E_volcano_candidate_counts.csv
- Fig6E_note.csv
- Fig6FG_GO_all_long.csv
- Fig6F_GO_terms_selected.csv
- Fig6G_GO_terms_selected.csv
- Fig6H_image_path_inventory.csv
- Fig6I_note.csv
- Fig6_file_mapping.csv
- Fig6_panels_data_availability.csv

Notes:
- Panels C/D now use original archives: ../fig6c.gz and ../fig6d.gz.
- Panel E uses group1_result DEG table when available; threshold summary is in Fig6E_volcano_group1_threshold_summary.csv.
- If group1_result is unavailable, fallback candidate LCs diff tables are summarized in Fig6E_volcano_candidate_counts.csv.
- Panels F/G selected term bars are extracted from LCa_ctrl_vs_disease and LCa_disease_vs_ctrl GO.csv files.

Build scripts:
- ../fig6_code/fig6_cd_prepare_source_data.py
- ../fig6_code/fig6_cd_dotplot.py
- ../fig6_code/fig6_source_data.py
