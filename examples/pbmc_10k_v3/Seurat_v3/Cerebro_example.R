##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
.libPaths(
  c(
    "/usr/local/lib/R/site-library",
    "/usr/lib/R/site-library",
    "/usr/lib/R/library"
  )
)

setwd("/data")

set.seed(1234567)

##----------------------------------------------------------------------------##
## Load libraries.
##----------------------------------------------------------------------------##
library("Seurat")
library("cerebroPrepare")
library("tidyverse")

##----------------------------------------------------------------------------##
## Load count matrix.
##----------------------------------------------------------------------------##
path_to_data <- "./raw_data"

feature_matrix <- Matrix::readMM(paste0(path_to_data, "/matrix.mtx.gz"))
feature_matrix <- as.matrix(feature_matrix)
feature_matrix <- as.data.frame(feature_matrix)

colnames(feature_matrix) <- read_tsv(paste0(path_to_data, "/barcodes.tsv.gz"), col_names = FALSE) %>%
  select(1) %>%
  t() %>%
  as.vector()

gene_names <- read_tsv(paste0(path_to_data, "/features.tsv.gz"), col_names = FALSE) %>%
  select(2) %>%
  t() %>%
  as.vector()

feature_matrix <- feature_matrix %>% 
  mutate(gene = gene_names) %>%
  select("gene", everything()) %>%
  group_by(gene) %>%
  summarise_all(sum)

genes <- feature_matrix$gene

feature_matrix <- select(feature_matrix, -c("gene"))
feature_matrix <- as.data.frame(feature_matrix)
rownames(feature_matrix) <- genes

##----------------------------------------------------------------------------##
## Basic Seurat analysis.
##----------------------------------------------------------------------------##
seurat <- CreateSeuratObject(
  project = "PBMC_10k_v3",
  counts = feature_matrix,
  min.cells = 10
)

seurat <- subset(seurat, subset = nCount_RNA > 100 & nFeature_RNA > 50)

seurat <- NormalizeData(seurat)

seurat <- FindVariableFeatures(seurat)

seurat <- ScaleData(seurat, vars.to.regress = "nCount_RNA")

seurat <- CellCycleScoring(
  seurat,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)

seurat <- RunPCA(
  seurat,
  npcs = 30,
  features = seurat@assays$RNA@var.features
)

seurat <- FindNeighbors(seurat)
seurat <- FindClusters(seurat, resolution = 0.5)

seurat <- BuildClusterTree(
  seurat,
  dims = 1:30,
  reorder = TRUE,
  reorder.numeric = TRUE
)

seurat[["cluster"]] <- factor(
  as.character(seurat@meta.data$tree.ident),
  levels = sort(unique(seurat@meta.data$tree.ident))
)

seurat@meta.data$RNA_snn_res.0.7 <- NULL
seurat@meta.data$tree.ident <- NULL

##----------------------------------------------------------------------------##
## Dimensional reductions.
##----------------------------------------------------------------------------##
seurat <- RunTSNE(
  seurat,
  reduction.name = "tSNE",
  dims = 1:30,
  dim.embed = 2,
  perplexity = 30,
  seed.use = 100
)

seurat <- RunTSNE(
  seurat,
  reduction.name = "tSNE_3D",
  dims = 1:30,
  dim.embed = 3,
  perplexity = 30,
  seed.use = 100
)

seurat <- RunUMAP(
  seurat,
  reduction.name = "UMAP",
  dims = 1:30,
  n.components = 2,
  seed.use = 100
)

seurat <- RunUMAP(
  seurat,
  reduction.name = "UMAP_3D",
  dims = 1:30,
  n.components = 3,
  seed.use = 100
)

##----------------------------------------------------------------------------##
## Randomly assign cells to three samples.
##----------------------------------------------------------------------------##
sample_dump <- c(
  rep("PBMC_5k_1", 10000),
  rep("PBMC_5k_2", 10000),
  rep("PBMC_5k_3", 10000)
)

seurat[["sample"]] <- sample(sample_dump, nrow(seurat@meta.data))

seurat@meta.data$sample <- factor(
  seurat@meta.data$sample,
  levels = c("PBMC_5k_1","PBMC_5k_2","PBMC_5k_3")
)

##----------------------------------------------------------------------------##
## Store meta data and parameters of analysis.
##----------------------------------------------------------------------------##
seurat@misc$experiment <- list(
  experiment_name = "PBMC_10k",
  organism = "hg",
  date_of_analysis = Sys.Date()
)

seurat@misc$parameters <- list(
  gene_nomenclature = "gene_name",
  discard_genes_expressed_in_fewer_cells_than = 10,
  keep_mitochondrial_genes = TRUE,
  variables_to_regress_out = "nUMI",
  number_PCs = 30,
  tSNE_perplexity = 30,
  cluster_resolution = 0.7
)

seurat@misc$parameters$filtering <- list(
  UMI_min = 100,
  UMI_max = Inf,
  genes_min = 50,
  genes_max = Inf
)

##----------------------------------------------------------------------------##
## Perform further analysis with cerebroPrepare and export file for Cerebro.
##----------------------------------------------------------------------------##
seurat <- cerebroPrepare::addPercentMtRibo(
  seurat,
  organism = "hg",
  gene_nomenclature = "name"
)

seurat <- cerebroPrepare::getMostExpressedGenes(
  seurat,
  column_sample = "sample",
  column_cluster = "cluster"
)

seurat <- cerebroPrepare::getMarkerGenes(
  seurat,
  organism = "hg",
  column_sample = "sample",
  column_cluster = "cluster"
)

seurat <- cerebroPrepare::getEnrichedPathways(
  seurat,
  column_sample = "sample",
  column_cluster = "cluster",
  adj_p_cutoff = 0.01,
  max_terms = 100
)

cerebroPrepare::exportFromSeurat(
  seurat,
  experiment_name = "PBMC_10k",
  file = paste0("Seurat_v3/cerebro_PBMC_10k_", Sys.Date(), ".crb"),
  organism = "hg",
  column_nUMI = "nCount_RNA",
  column_nGene = "nFeature_RNA",
  column_cell_cycle_seurat = "Phase"
)

##----------------------------------------------------------------------------##
## Save Seurat object.
##----------------------------------------------------------------------------##
saveRDS(seurat, "Seurat_v3/seurat.rds")

