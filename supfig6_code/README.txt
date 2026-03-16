SupFigure6 code files

1) supfig6_source_data.py
- Build SupFigure6 source-data package for Panels A-C.
- Output directory:
  ../supfig6_source_data/
- Generates:
  - panel-level source CSV tables
  - raw_support file copies
  - SupFig6_file_mapping.csv
  - SupFig6_panels_data_availability.csv

2) supfig6_panelA_heatmap.py
- Replot Panel A heatmap from:
  ../supfig6_source_data/SupFig6A_scmeta_heatmap_top30.csv
- Output:
  ./supfig6_plots/SupFig6A_heatmap_top30.png
  ./supfig6_plots/SupFig6A_heatmap_top30.svg

Run order:
1. python3 supfig6_source_data.py
2. python3 supfig6_panelA_heatmap.py
