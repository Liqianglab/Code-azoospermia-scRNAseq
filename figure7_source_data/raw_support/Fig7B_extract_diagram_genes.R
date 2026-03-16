#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

base_dir <- "/Users/xq/Desktop/睾丸文章总/最新内容汇总/figure汇总/figure7"
expr_path <- file.path(base_dir, "late阶段所有基因表达.csv")
out_path <- file.path(base_dir, "late阶段_代谢通路图基因表达_提取.csv")

# Genes shown in the metabolic diagram (glycolysis, TCA, OXPHOS)
genes_tbl <- tibble::tribble(
  ~pathway, ~subcategory, ~gene, ~diagram_order,
  # Glycolysis
  "Glycolysis", "Glucose transport", "SLC2A1", 1,
  "Glycolysis", "Glucose transport", "SLC2A4", 2,
  "Glycolysis", "Upper glycolysis", "HK1", 3,
  "Glycolysis", "Upper glycolysis", "HK2", 4,
  "Glycolysis", "Upper glycolysis", "HK3", 5,
  "Glycolysis", "Upper glycolysis", "GPI", 6,
  "Glycolysis", "Rate-limiting", "PFKM", 7,
  "Glycolysis", "Rate-limiting", "PFKP", 8,
  "Glycolysis", "Middle glycolysis", "ALDOA", 9,
  "Glycolysis", "Middle glycolysis", "ALDOC", 10,
  "Glycolysis", "Middle glycolysis", "TPI1", 11,
  "Glycolysis", "Lower glycolysis", "PGK1", 12,
  "Glycolysis", "Lower glycolysis", "PGAM1", 13,
  "Glycolysis", "Lower glycolysis", "ENO1", 14,
  "Glycolysis", "Lower glycolysis", "ENO2", 15,
  "Glycolysis", "Pyruvate kinase", "PKLR", 16,
  "Glycolysis", "Pyruvate kinase", "PKM", 17,
  "Glycolysis", "Lactate dehydrogenase", "LDHA", 18,
  "Glycolysis", "Lactate dehydrogenase", "LDHB", 19,
  # Link to TCA (PDH)
  "TCA cycle", "PDH complex", "DLAT", 20,
  "TCA cycle", "PDH complex", "DLD", 21,
  "TCA cycle", "PDH complex", "PDHB", 22,
  # TCA cycle
  "TCA cycle", "Entry/Anaplerosis", "CS", 23,
  "TCA cycle", "Entry/Anaplerosis", "PC", 24,
  "TCA cycle", "Isocitrate dehydrogenase", "IDH1", 25,
  "TCA cycle", "Isocitrate dehydrogenase", "IDH2", 26,
  "TCA cycle", "Isocitrate dehydrogenase", "IDH3A", 27,
  "TCA cycle", "α-KGDH complex", "DLST", 28,
  "TCA cycle", "α-KGDH complex", "OGDH", 29,
  "TCA cycle", "Succinyl-CoA synthetase", "SUCLG1", 30,
  "TCA cycle", "Succinate dehydrogenase", "SDHB", 31,
  "TCA cycle", "Fumarase", "FH", 32,
  # Oxidative phosphorylation
  "Oxidative phosphorylation", "CoQ", "COQ2", 33,
  "Oxidative phosphorylation", "CoQ", "COQ4", 34,
  "Oxidative phosphorylation", "CoQ", "COQ6", 35,
  "Oxidative phosphorylation", "Complex I", "MT-ND1", 36,
  "Oxidative phosphorylation", "Complex I", "MT-ND3", 37,
  "Oxidative phosphorylation", "Complex I", "NDUFV1", 38,
  "Oxidative phosphorylation", "Complex I", "NDUFA2", 39,
  "Oxidative phosphorylation", "Complex I", "NDUFB1", 40,
  "Oxidative phosphorylation", "Complex I", "NDUFC1", 41,
  "Oxidative phosphorylation", "Complex II", "SDHAF1", 42,
  "Oxidative phosphorylation", "Complex II", "SDHAF2", 43,
  "Oxidative phosphorylation", "Complex II", "SDHA", 44,
  "Oxidative phosphorylation", "Complex II", "SDHC", 45,
  "Oxidative phosphorylation", "Complex II", "SDHD", 46,
  "Oxidative phosphorylation", "Complex IV", "MT-CO1", 47,
  "Oxidative phosphorylation", "Complex IV", "MT-CO2", 48,
  "Oxidative phosphorylation", "Complex IV", "MT-CO3", 49,
  "Oxidative phosphorylation", "Complex IV", "COX4I1", 50,
  "Oxidative phosphorylation", "Complex IV", "COX6A1", 51,
  "Oxidative phosphorylation", "Complex IV", "COX6B2", 52,
  "Oxidative phosphorylation", "Complex V", "ATP5F1B", 53,
  "Oxidative phosphorylation", "Complex V", "ATP5F1D", 54,
  "Oxidative phosphorylation", "Complex V", "ATP5MC1", 55,
  "Oxidative phosphorylation", "Complex V", "ATP5ME", 56,
  "Oxidative phosphorylation", "Complex V", "ATP5MF", 57,
  "Oxidative phosphorylation", "Complex V", "ATP5PF", 58
)

expr <- read_csv(expr_path, show_col_types = FALSE) %>%
  rename(gene = Gene, ctrl = CTRL, disease = Disease)

eps <- 1e-6
out <- genes_tbl %>%
  left_join(expr, by = "gene") %>%
  mutate(
    log2fc = log2((disease + eps) / (ctrl + eps))
  ) %>%
  arrange(diagram_order)

missing <- out %>% filter(is.na(ctrl) | is.na(disease)) %>% pull(gene)
if (length(missing) > 0) {
  warning("Missing genes in expression table: ", paste(missing, collapse = ", "))
}

write_csv(out, out_path)
message("Wrote: ", out_path)

