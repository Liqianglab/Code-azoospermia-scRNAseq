#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ---- resolve script root ----
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
root_dir <- if (length(script_path) > 0) dirname(dirname(script_path)) else "."

# ---- user inputs ----
# Update these paths to your raw 10X directories
samples <- tibble::tribble(
  ~sample_id,   ~group,      ~path,
  "Ctrl_1",     "Ctrl",      "path/to/Ctrl_1/filtered_feature_bc_matrix",
  "Ctrl_2",     "Ctrl",      "path/to/Ctrl_2/filtered_feature_bc_matrix",
  "Ctrl_3",     "Ctrl",      "path/to/Ctrl_3/filtered_feature_bc_matrix",
  "OA_2",       "OA",        "path/to/OA_2/filtered_feature_bc_matrix",
  "OA_3",       "OA",        "path/to/OA_3/filtered_feature_bc_matrix",
  "OA_4",       "OA",        "path/to/OA_4/filtered_feature_bc_matrix",
  "AZFc_Del_1","AZFc_Del",  "path/to/AZFc_Del_1/filtered_feature_bc_matrix",
  "AZFc_Del_2","AZFc_Del",  "path/to/AZFc_Del_2/filtered_feature_bc_matrix",
  "AZFc_Del_3","AZFc_Del",  "path/to/AZFc_Del_3/filtered_feature_bc_matrix",
  "AZFc_Del_4","AZFc_Del",  "path/to/AZFc_Del_4/filtered_feature_bc_matrix",
  "KS_1",       "KS",        "path/to/KS_1/filtered_feature_bc_matrix",
  "KS_2",       "KS",        "path/to/KS_2/filtered_feature_bc_matrix",
  "iNOA_B",     "iNOA_B",    "path/to/iNOA_B/filtered_feature_bc_matrix",
  "iNOA_S1",    "iNOA_S",    "path/to/iNOA_S1/filtered_feature_bc_matrix",
  "iNOA_S2",    "iNOA_S",    "path/to/iNOA_S2/filtered_feature_bc_matrix",
  "iNOA_S3",    "iNOA_S",    "path/to/iNOA_S3/filtered_feature_bc_matrix"
)

# QC thresholds (adjust to match your pipeline)
qc_min_features <- 200
qc_max_features <- 7500
qc_max_mt <- 20

# Marker list for Sup Fig 1A
markers <- c(
  "BEX1","WT1","PECAM1","CDH5","INSL3","IGF1",
  "CD163","CD14","FAM129A","ACTA2","UTF1","PIWIL4",
  "MKI67","AURKA","STRA8","DMC1","BRCA1","PSMC3IP",
  "PRM2","SPEM1"
)

# Cluster -> cell type mapping (update if cluster IDs change)
cluster_to_celltype <- c(
  `1`="SSCs", `2`="STs", `3`="LCs", `4`="Spermatids", `5`="SPCs",
  `6`="SPCs", `7`="Myeloid", `8`="STs", `9`="SPGs", `10`="Spermatids",
  `11`="Myeloid", `12`="SSCs", `13`="Spermatids", `14`="Myeloid", `15`="Spermatids",
  `16`="LCs", `17`="SPCs", `18`="ECs", `19`="Spermatids", `20`="Spermatids",
  `21`="Spermatids", `22`="Spermatids", `23`="SPCs", `24`="SPCs", `25`="STs",
  `26`="SPCs", `27`="Myeloid", `28`="Spermatids", `29`="LCs", `30`="Spermatids",
  `31`="LCs", `32`="SPCs", `33`="Myeloid", `34`="Spermatids", `35`="STs",
  `36`="Lym", `37`="Myeloid"
)

out_dir <- file.path(root_dir, "supfig1_source_data_from_raw")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- build Seurat objects ----
objs <- lapply(seq_len(nrow(samples)), function(i) {
  s <- samples[i,]
  mat <- Read10X(s$path)
  obj <- CreateSeuratObject(counts = mat, project = s$sample_id)
  obj$sample_id <- s$sample_id
  obj$group <- s$group
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_RNA > qc_min_features & nFeature_RNA < qc_max_features & percent.mt < qc_max_mt)
  obj
})

# ---- integration (SCT) ----
objs <- lapply(objs, SCTransform, verbose = FALSE)
features <- SelectIntegrationFeatures(objs, nfeatures = 3000)
objs <- PrepSCTIntegration(objs, anchor.features = features)
anchors <- FindIntegrationAnchors(objs, normalization.method = "SCT", anchor.features = features)
obj <- IntegrateData(anchors, normalization.method = "SCT")

# ---- dimensionality reduction & clustering ----
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)

# ---- assign cell types ----
obj$cell_type <- unname(cluster_to_celltype[as.character(obj$seurat_clusters)])

# ---- Sup Fig 1A (feature plot source data) ----
markers_present <- intersect(markers, rownames(obj))
feat <- FetchData(obj, vars = c("UMAP_1", "UMAP_2", markers_present))
feat$cell_id <- rownames(feat)
meta <- obj@meta.data %>%
  select(sample_id, group, cell_type) %>%
  mutate(cell_id = rownames(obj@meta.data))

sup1a <- feat %>%
  left_join(meta, by = "cell_id") %>%
  pivot_longer(cols = all_of(markers_present), names_to = "gene", values_to = "expression")
write_csv(sup1a, file.path(out_dir, "SupFig1A_featureplot_long.csv"))

# ---- Sup Fig 1B (heatmap source data) ----
# Option 1: use a predefined gene list (e.g., from average_expression_export.csv)
heatmap_genes <- NULL
heatmap_candidates <- c(
  file.path(root_dir, "Average expression Heatmap", "average_expression_export.csv"),
  file.path(root_dir, "figure1", "Average expression Heatmap", "average_expression_export.csv")
)
heatmap_path <- heatmap_candidates[file.exists(heatmap_candidates)][1]
if (!is.na(heatmap_path) && file.exists(heatmap_path)) {
  tmp <- read_csv(heatmap_path, show_col_types = FALSE)
  heatmap_genes <- setdiff(names(tmp), "cluster")
}

# Option 2: derive markers if no predefined list is available
if (is.null(heatmap_genes)) {
  markers_all <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  heatmap_genes <- markers_all %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10) %>%
    pull(gene) %>%
    unique()
}

avg_expr <- AverageExpression(obj, group.by = "cell_type", features = heatmap_genes, return.seurat = FALSE)$RNA
avg_expr <- as.data.frame(avg_expr)
avg_expr$gene <- rownames(avg_expr)

avg_long <- avg_expr %>%
  pivot_longer(cols = -gene, names_to = "cell_type", values_to = "value")
write_csv(avg_long, file.path(out_dir, "SupFig1B_heatmap_avgexpr_long.csv"))

# Scaled heatmap (z-score by gene)
scaled <- t(scale(t(as.matrix(avg_expr[, setdiff(colnames(avg_expr), "gene")]))))
scaled <- as.data.frame(scaled)
scaled$gene <- avg_expr$gene
scaled_long <- scaled %>%
  pivot_longer(cols = -gene, names_to = "cell_type", values_to = "value")
write_csv(scaled_long, file.path(out_dir, "SupFig1B_heatmap_scaled_long.csv"))

# ---- Sup Fig 1B (GO enrichment: optional) ----
# Requires clusterProfiler + org.Hs.eg.db
# Uncomment if you want to run GO directly from raw data
# library(clusterProfiler)
# library(org.Hs.eg.db)
# markers_all <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# go_results <- list()
# for (ct in unique(obj$cell_type)) {
#   genes <- markers_all %>% filter(cluster == ct) %>% pull(gene) %>% unique()
#   eg <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
#   if (nrow(eg) == 0) next
#   ego <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
#                   pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
#   res <- as.data.frame(ego)
#   res$group <- ct
#   go_results[[ct]] <- res
# }
# if (length(go_results) > 0) {
#   go_df <- bind_rows(go_results)
#   write_csv(go_df, file.path(out_dir, "SupFig1B_GO_all_terms.csv"))
# }

saveRDS(obj, file.path(out_dir, "supfig1_seurat.rds"))
message("Done. Outputs in ", out_dir)
