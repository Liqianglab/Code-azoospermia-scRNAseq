Figure 4 source data (Panels A-I)

Generated tables:
- Fig4A_ST_cells_trajectory.csv
- Fig4A_trajectory_line_points.csv
- Fig4A_stage_candidate_mapping.csv
- Fig4A_ST_counts_by_group.csv
- Fig4A_ST_counts_by_state.csv
- Fig4A_ST_pseudotime_summary.csv
- Fig4B_heatmap_matrix_wide.csv
- Fig4B_heatmap_matrix_long.csv
- Fig4B_heatmap_gene_order.csv
- Fig4C_stage_marker_cells.csv
- Fig4C_stage_marker_summary.csv
- Fig4D_stage_counts_by_group.csv
- Fig4D_stage_ratio_by_group.csv
- Fig4D_stage_ratio_all_group.csv
- Fig4D_stage_counts_selected_samples.csv
- Fig4D_stage_ratio_selected_samples.csv
- Fig4E_group_marker_cells.csv
- Fig4E_group_marker_summary.csv
- Fig4F_regulon_activity_matrix_wide.csv
- Fig4F_regulon_activity_matrix_long.csv
- Fig4F_regulon_candidates_from_scaled_top.csv
- Fig4F_note.csv
- Fig4G_BTB_score_cells.csv
- Fig4G_BTB_group_score_summary.csv
- Fig4G_BTB_stage_score_summary.csv
- Fig4G_note.csv
- Fig4HI_image_panels_note.csv
- Fig4_file_mapping.csv
- Fig4_panels_data_availability.csv

Notes:
- ST stage labels in this export use candidate mapping ST1->Stage_a, ST2->Stage_b, ST3->Stage_c (see Fig4A_stage_candidate_mapping.csv).
- Panel F uses raw pyscenic regulon matrices extracted from ST细胞提取-scenic1_regulon.tar.gz, and excludes OA to match Figure 4.
- Panel G uses raw cell-level BTB integrity score table from ST细胞提取 - BTB_Integrity_Score_expression.csv, and excludes OA to match Figure 4.

Build script:
- ../figure4_code/figure4_source_data.py