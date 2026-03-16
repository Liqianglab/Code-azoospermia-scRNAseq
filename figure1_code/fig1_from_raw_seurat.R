#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
})

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

# Marker list for Fig1D/G
markers <- c(
  "VIM","SOX9","VWF","DLK1","PTPRC","NOTCH3",
  "DDX4","FGFR3","KIT","SYCP3","TEX29","PRM3",
  "CPEB1","MAGEB2","DMRT1","SOHLH1","SOHLH2","BRCA1","DMC1",
  "MLH3","PSMC3IP","WT1","CLDN11","AMH","S100A13","ENO1","BEX1",
  "CST9L","HOPX","SPAG6","ACTL7B","ATG4D","ZPBP","DNAH6","SPATA18",
  "HOOK1","TNP1","NR2F2","CFD","CLDN5","CDH5","PECAM1","LYZ",
  "CD14","MRC1","C1QC","AIF1","CD3D","CD79A","ACTA2","MYH11"
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

out_dir <- file.path("source data", "fig1_source_data_from_raw")
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

# ---- Fig1B (UMAP by group) ----
emb <- Embeddings(obj, "umap") %>% as.data.frame()
emb$cell_id <- rownames(emb)
fig1b <- emb %>%
  left_join(obj@meta.data %>%
              select(sample_id, group) %>%
              mutate(cell_id = rownames(obj@meta.data)), by = "cell_id") %>%
  select(cell_id, sample_id, group, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
write_csv(fig1b, file.path(out_dir, "Fig1B_UMAP_by_group.csv"))

# ---- Fig1C (UMAP by cell type) ----
fig1c <- emb %>%
  left_join(obj@meta.data %>%
              select(sample_id, group, cell_type) %>%
              mutate(cell_id = rownames(obj@meta.data)), by = "cell_id") %>%
  select(cell_id, cell_type, sample_id, group, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
write_csv(fig1c, file.path(out_dir, "Fig1C_UMAP_by_celltype.csv"))

# ---- Fig1D (marker feature data, long) ----
feat <- FetchData(obj, vars = c("UMAP_1", "UMAP_2", markers))
feat$cell_id <- rownames(feat)
meta <- obj@meta.data %>%
  select(sample_id, group) %>%
  mutate(cell_id = rownames(obj@meta.data))
fig1d <- feat %>%
  left_join(meta, by = "cell_id") %>%
  pivot_longer(cols = all_of(markers), names_to = "gene", values_to = "expression")
write_csv(fig1d, file.path(out_dir, "Fig1D_marker_feature_long.csv"))

# ---- Fig1G (dotplot stats by cell type) ----
expr_mat <- GetAssayData(obj, slot = "data")
expr_mat <- expr_mat[intersect(markers, rownames(expr_mat)), , drop = FALSE]
meta_ct <- obj@meta.data$cell_type

avg_expr <- apply(expr_mat, 1, function(g) tapply(g, meta_ct, mean))
avg_expr <- as.data.frame(avg_expr)
avg_expr$cell_type <- rownames(avg_expr)

pct_expr <- apply(expr_mat > 0, 1, function(g) tapply(g, meta_ct, mean))
pct_expr <- as.data.frame(pct_expr)
pct_expr$cell_type <- rownames(pct_expr)

fig1g <- avg_expr %>%
  pivot_longer(-cell_type, names_to = "gene", values_to = "avg_expr") %>%
  left_join(pct_expr %>% pivot_longer(-cell_type, names_to = "gene", values_to = "pct_expr"),
            by = c("cell_type", "gene"))
write_csv(fig1g, file.path(out_dir, "Fig1G_dotplot_full_genes.csv"))

# ---- Fig1E/F (composition) ----
comp <- obj@meta.data %>%
  count(group, cell_type, name = "n_cells") %>%
  group_by(group) %>%
  mutate(percent = n_cells / sum(n_cells)) %>%
  ungroup()
write_csv(comp, file.path(out_dir, "Fig1E_group_composition_long.csv"))

somatic <- c("ECs","LCs","Lym","Myeloid","Myoid")
germ <- c("SPCs","SPGs","SSCs","STs","Spermatids")

fig1f <- obj@meta.data %>%
  mutate(cell_group = case_when(
    cell_type %in% somatic ~ "somatic",
    cell_type %in% germ ~ "germ",
    TRUE ~ "other"
  )) %>%
  count(sample_id, group, cell_group, name = "n_cells") %>%
  filter(cell_group %in% c("somatic","germ")) %>%
  pivot_wider(names_from = cell_group, values_from = n_cells) %>%
  mutate(age = NA_real_)
write_csv(fig1f, file.path(out_dir, "Fig1F_sample_somatic_germ_counts.csv"))

# Save Seurat object for plotting
saveRDS(obj, file.path(out_dir, "fig1_seurat.rds"))

message("Done. Outputs in ", out_dir)
