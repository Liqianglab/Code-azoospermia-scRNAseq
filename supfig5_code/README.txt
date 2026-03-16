SupFigure5 code files

1) supfig5_source_data.py
- Build SupFigure5 source-data CSV files from:
  ../supfig5_source_data/raw_support/
- Output CSV tables to:
  ../supfig5_source_data/
- Panel A behavior:
  prefer raw score tables (three Gene set enrichment CSV files),
  fallback to legacy proxy reconstruction only when raw tables are missing.

2) supfig5_panelAB_plot.py
- Plot panel A/B from generated CSV tables.
- Output:
  ./supfig5_plots/SupFig5A_violin.svg
  ./supfig5_plots/SupFig5A_violin_proxy.svg
  ./supfig5_plots/SupFig5B_GO_updown_top15.svg
  ./supfig5_plots/SupFig5AB_layout.svg
  ./supfig5_plots/SupFig5AB_proxy_layout.svg

3) raw_original_scripts/
- Original script recovered from LC-GO top15 bundle:
  SupFig5B_GO_top15_original.R

Run order:
1. python3 supfig5_source_data.py
2. python3 supfig5_panelAB_plot.py
