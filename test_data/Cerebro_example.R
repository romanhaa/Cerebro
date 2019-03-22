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

setwd("/data/pbmc_10k_v3")

set.seed(1234567)

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
library("cerebroPrepare")
library("magrittr")
library("Seurat")
library("tidyverse")

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
path_to_data <- "."

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

feature_matrix %<>% 
  mutate(gene = gene_names) %>%
  select("gene", everything()) %>%
  group_by(gene) %>%
  summarise_all(sum)

genes <- feature_matrix$gene

feature_matrix <- select(feature_matrix, -c("gene"))
feature_matrix <- as.data.frame(feature_matrix)
rownames(feature_matrix) <- genes

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
seurat <- CreateSeuratObject(
  project = "PBMC_10k_v3",
  raw.data = feature_matrix,
  min.cells = 10
)

seurat <- FilterCells(
  seurat,
  subset.names = c("nUMI", "nGene"),
  low.thresholds = c(100, 50),
  high.thresholds = c(Inf, Inf)
)

seurat <- NormalizeData(seurat)

seurat <- FindVariableGenes(seurat)
length(seurat@var.genes)
# 500

seurat <- ScaleData(seurat, vars.to.regress = "nUMI")

seurat <- CellCycleScoring(
  seurat,
  g2m.genes = cc.genes$g2m.genes,
  s.genes = cc.genes$s.genes
)

seurat <- RunPCA(
  seurat,
  pcs.compute = 30,
  pc.genes = seurat@var.genes
)

seurat <- FindClusters(
  seurat,
  reduction.type = "pca",
  dims.use = 1:30,
  resolution = 0.7,
  print.output = FALSE,
  save.SNN = TRUE
)

seurat <- BuildClusterTree(
  seurat,
  pcs.use = 1:30,
  do.plot = FALSE,
  do.reorder = TRUE,
  reorder.numeric = TRUE
)

seurat <- AddMetaData(
  seurat,
  metadata = factor(seurat@ident, ordered = FALSE),
  col.name = "cluster"
)

seurat@meta.data$res.0.7 <- NULL

seurat <- RunTSNE(
  seurat,
  reduction.name = "tSNE",
  reduction.key = "tSNE",
  dims.use = 1:30,
  do.fast = TRUE,
  seed.use = 100
)

seurat <- RunTSNE(
  seurat,
  reduction.name = "tSNE_3D",
  reduction.key = "tSNE",
  dims.use = 1:30,
  dim.embed = 3,
  do.fast = TRUE,
  seed.use = 100
)

seurat <- RunUMAP(
  seurat,
  reduction.name = "UMAP",
  dims.use = 1:30,
  do.fast = TRUE,
  seed.use = 100
)

seurat <- RunUMAP(
  seurat,
  reduction.name = "UMAP_3D",
  dims.use = 1:30,
  max.dim = 3,
  do.fast = TRUE,
  seed.use = 100
)

seurat <- RunDiffusion(
  seurat,
  reduction.name = "DM",
  dims.use = 1:30
)

seurat <- RunDiffusion(
  seurat,
  reduction.name = "DM_3D",
  dims.use = 1:30,
  max.dim = 3
)

seurat <- RunPHATE(
  seurat,
  reduction.name = "PHATE",
  npca = 30,
  max.dim = 2,
  n.jobs = -1
)

seurat <- RunPHATE(
  seurat,
  reduction.name = "PHATE_3D",
  npca = 30,
  max.dim = 3,
  n.jobs = -1
)

##----------------------------------------------------------------------------##
## Randomly assign cells to three samples.
##----------------------------------------------------------------------------##
sample_dump <- c(
  rep("PBMC_5k_1", 10000),
  rep("PBMC_5k_2", 10000),
  rep("PBMC_5k_3", 10000)
)

sample <- data.frame(
  "sample" = sample(sample_dump, nrow(seurat@meta.data)),
  row.names = seurat@cell.names
)

seurat <- AddMetaData(
  seurat,
  metadata = sample,
  col.name = "sample"
)

seurat@meta.data$sample <- factor(
  seurat@meta.data$sample,
  levels = c("PBMC_5k_1","PBMC_5k_2","PBMC_5k_3")
)

##----------------------------------------------------------------------------##
##
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
  column_cluster = "cluster",
  only.pos = TRUE,
  min.pct = 0.7,
  thresh.use = 0.25,
  test.use = "t"
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
  file = paste0("cerebro_PBMC_10k_", Sys.Date(), ".crb"),
  organism = "hg",
  column_cell_cycle_seurat = "Phase"
)

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
saveRDS(seurat, "seurat.rds")

