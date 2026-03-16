Bundle: GeneScores 8-panels (cell-level p-values)

Files:
- ExternalValidation_GSE149512_GeneScores_8panels_cellPvalues.pdf/.png
- ExternalValidation_GSE235321_GeneScores_8panels_cellPvalues.pdf/.png
- ExternalValidation_GeneScores_8panels_cellPvalues_BOTH.pdf   (2 pages)
- GSE149512_cellLevel_geneScores8_withCellPvalues_data.csv      (intermediate cell-level table)
- GSE235321_cellLevel_geneScores8_withCellPvalues_data.csv      (intermediate cell-level table)
- GSE149512_geneScores8_cellLevel_pvalues.csv                   (cell-level MWU p-values)
- GSE235321_geneScores8_cellLevel_pvalues.csv                   (cell-level MWU p-values)
- external_validation_geneScores8_cellPvalues.py                (reproducible code)
- signature_manifest_externalValidation.csv                     (gene set manifest; from earlier steps)

Notes:
- P-values are computed using Mann–Whitney U (two-sided) at the cell level.
  This can inflate significance due to pseudo-replication; donor medians are overlaid visually.
- If you need donor-level p-values (more appropriate), you can adapt the script to test donor medians.
