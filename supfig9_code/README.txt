Supplementary Figure 9 code

Files:
- supfig9_source_data.py
  Builds panel-wise source-data tables, copies raw support files, writes file mapping,
  panel availability, and panel-source_code-confidence checklist for Supplementary Figure 9 (A-K).

- supfig9_plot_panels.py
  Replots check figures for Panels A-K from generated source data.
  Outputs to: ./supfig9_plots/

Run order:
1) python3 source data/supfig9_code/supfig9_source_data.py
2) python3 source data/supfig9_code/supfig9_plot_panels.py
