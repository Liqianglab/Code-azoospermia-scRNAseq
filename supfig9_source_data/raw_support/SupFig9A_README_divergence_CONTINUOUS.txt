External validation divergence plots (CONTINUOUS version)

This bundle contains:
1) Final plots (2:1 aspect, unified colors)
   - ExternalValidation_GSE235321_developmental_divergence_entropy_FINAL_CONTINUOUS.pdf/png
   - ExternalValidation_GSE149512_developmental_divergence_entropy_FINAL_CONTINUOUS.pdf/png

2) Intermediate data
   - GSE235321_externalValidation_cellLevel_pseudotime_entropy_moduleScores.csv
     (cell-level pseudotime2 + sliding-window mixing entropy + stage labels)
   - GSE149512_crossPlatformIntegrated_cellLevel_pseudotime_entropy.csv
     (cell-level pseudotime_scaled + mixing entropy after platform-corrected PCs)

   - Smoothed curves:
     * GSE235321_developmental_divergence_entropy_smoothedCurve_FINAL_CONTINUOUS.csv
     * GSE149512_developmental_divergence_entropy_smoothedCurve_FINAL_CONTINUOUS.csv

3) Code
   - external_validation_divergence_plot_unified_CONTINUOUS.py

Run:
  python external_validation_divergence_plot_unified_CONTINUOUS.py

Font note:
  Uses sans-serif with fallback [Arial, Liberation Sans, DejaVu Sans]. If Arial is absent, swap in Illustrator.
