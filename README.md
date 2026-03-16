# scAzoo

Code and curated source-data tables for a single-cell RNA-seq study of azoospermia and related male infertility subtypes.

This repository is organized by manuscript figure. Each figure has a code directory and a matching source-data directory that stores exported CSV tables, figure notes, and archived support files when available. The figure set covers atlas construction, developmental trajectory analysis, somatic niche remodeling, metabolic programs, and external-validation analyses.

## Manuscript

This repository is intended as the companion code and source-data archive for the associated azoospermia scRNA-seq manuscript. The formal paper title, author list, DOI, and primary dataset accession numbers can be added here after submission or publication.

## Repository At A Glance

- Scope: 8 main figures and 9 supplementary figures.
- Structure: figure-oriented code and source-data packages, not one monolithic end-to-end pipeline.
- Data included: curated source-data CSVs, plotting scripts, source-data builders, panel-wise raw assets, and recovered original scripts where available.
- Groups appearing across figures: `Ctrl`, `OA`, `AZFc_Del`, `iNOA_B`, `iNOA_S`, `KS`; some external-validation panels also use `AZFa` and `iNOA`.

## Repository Structure

| Path | Description |
| --- | --- |
| `figure1_code/` ... `figure8_code/` | Scripts for main-figure source-data generation or figure reconstruction |
| `figure1_source_data/` ... `figure8_source_data/` | Curated CSV tables, figure-specific notes, and `raw_support/` inputs |
| `supfig1_code/` ... `supfig9_code/` | Scripts for supplementary figures |
| `supfig1_source_data/` ... `supfig9_source_data/` | Curated supplementary source-data tables and notes |
| `figure4_panel_raw/` | Panel-wise raw package for Figure 4, including recovered plotting snippets when found |
| `*_plots/` | Reconstructed panel or figure outputs produced by selected scripts |
| `raw_support/` | Archived intermediate/raw inputs copied into the public package when available |

## How To Navigate This Repository

If you want to inspect a specific figure:

1. Open the matching `figureX_source_data/` or `supfigX_source_data/` directory.
2. Read the local `README.txt` when present.
3. Use the CSV files in that directory as the canonical public source-data tables.
4. Check the sibling `figureX_code/` or `supfigX_code/` directory for reconstruction scripts.

Recommended entry points:

- Data behind a panel: the corresponding `figureX_source_data/` or `supfigX_source_data/` directory.
- Reconstruction code: the corresponding `figureX_code/` or `supfigX_code/` directory.
- Raw supporting inputs: `raw_support/` plus any `*_file_mapping.csv` or `*_panels_data_availability.csv`.
- Figure 4 raw panel assets: `figure4_panel_raw/`.

## Typical Usage

Many plotting scripts read directly from the curated `*_source_data` directories already included in this repository. Representative examples:

```bash
Rscript figure2_code/fig2_plot.R
python3 figure3_code/figure3_source_data.py
python3 figure3_code/fig3_plot_panels.py
python3 figure4_code/figure4_source_data.py
python3 figure4_code/fig4_plot_panels.py
Rscript supfig2_code/supfig2_plot.R
```

Outputs are typically written to a sibling `*_plots/` directory inside the corresponding `*_code` folder.

## Data Availability

- This repository already contains curated source-data tables for the manuscript figures and supplementary figures.
- It does not currently serve as a raw sequencing-data repository and does not include every author-local intermediate file referenced during figure reconstruction.
- Public accession links for the primary study dataset can be added here once they are finalized.
- Supplementary Figure 9 references bundled external-validation assets related to `GSE149512` and `GSE235321`.

## Software Requirements

No environment lockfile is currently included. Based on the scripts in this repository, the main software dependencies are:

- Python: `pandas`, `numpy`, `matplotlib`, `scipy`, `statsmodels`, `openpyxl`
- R: `Seurat`, `ggplot2`, `dplyr`, `tidyr`, `readr`, `stringr`, `cowplot`, `patchwork`, `scales`

Example setup:

```bash
python3 -m pip install pandas numpy matplotlib scipy statsmodels openpyxl
```

```r
install.packages(c(
  "Seurat", "ggplot2", "dplyr", "tidyr", "readr",
  "stringr", "cowplot", "patchwork", "scales"
))
```

## Reproducibility Notes

- This repository is best understood as a figure-oriented source-data and code archive rather than a fully frozen pipeline from raw sequencing files to final manuscript figures.
- The safest public starting point is the existing CSV data in each `*_source_data/` directory.
- Several reconstruction utilities still reference author-local paths such as `原来作图的数据/...` or `/Users/...`. These scripts document how source tables were assembled, but some of the raw inputs are not fully included in this repository.
- Some Python utilities still use legacy output names such as `fig1_source_data/`, `fig6_source_data/`, `fig7_source_data/`, or `fig8_source_data/`, while the current repository stores directories such as `figure1_source_data/`, `figure6_source_data/`, `figure7_source_data/`, and `figure8_source_data/`. Check and adjust paths before re-running builders.
- `figure4_panel_raw/` is a special case that stores per-panel raw materials and recovered original plotting code where available.
- Supplementary Figure 9 includes external-validation assets built around bundled files related to `GSE149512` and `GSE235321`.

## Manuscript And Citation

Manuscript title, author list, DOI, accession links, and release citation can be added here after submission or publication.

If you use this repository, please cite the associated manuscript and the archived repository release once available.
