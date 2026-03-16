SupFigure4 code files

1) supfig4_source_data.py
- Build SupFigure4 source-data CSV files (Panels A-F) from files in:
  ../supfig4_panel_raw/panel_*/raw_data/
- For panel D/F, also pull the original ST disease-vs-control GSEA files from:
  ../../../睾丸文章数据/ST/疾病的GSEA/
- Copies all panel raw files into:
  ../supfig4_source_data/raw_support/
- Removes stale legacy outputs (ST1vsST3/ST3vsST1 tables) in:
  ../supfig4_source_data/
- Packs all panel-level source tables into:
  ../supfig4_source_data/SourceData_SupFig4.xlsx

2) supfig4_panelB_heatmap_plot.py
- Plot Panel B regulon-by-stage heatmap SVG from:
  ../supfig4_source_data/SupFig4B_regulon_stage_matrix_target.csv
- Output:
  ./supfig4_plots/SupFig4B_regulon_stage_heatmap.svg

3) supfig4_panelDF_bubble_proxy_plot.py
- Plot proxy bubble charts for Panel D/F from:
  ../supfig4_source_data/SupFig4D_ST1_disease_vs_ctrl_top10.csv
  ../supfig4_source_data/SupFig4F_ST3_disease_vs_ctrl_top10.csv
- Output:
  ./supfig4_plots/SupFig4D_ST1_disease_vs_ctrl_proxy_bubble.svg
  ./supfig4_plots/SupFig4F_ST3_disease_vs_ctrl_proxy_bubble.svg

4) supfig4_panelDF_bubble_original.R
- Original bubble-plot script copied from ST disease GSEA folder.
- Supports custom input/output:
  Rscript supfig4_panelDF_bubble_original.R <input_csv> <output_pdf> <plot_title>

Run order:
1. python3 supfig4_source_data.py
2. python3 supfig4_panelB_heatmap_plot.py
3. python3 supfig4_panelDF_bubble_proxy_plot.py
