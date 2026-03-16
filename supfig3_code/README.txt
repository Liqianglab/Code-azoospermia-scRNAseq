SupFigure3 code files

1) supfig3_source_data.py
- Build SupFig3 source-data CSV files from original GSEA enrichment tables.
- Also copies the exact input `.xls` tables into `../supfig3_source_data/raw_support/`.

2) supfig3_plot_svg.py
- Render Panel A / Panel B bubble charts from the prepared matrix CSV files.
- Output SVG files to `./supfig3_plots/`.

Run order:
1. python3 supfig3_source_data.py
2. python3 supfig3_plot_svg.py
