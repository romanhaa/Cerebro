##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
.libPaths(
  c(
    '/usr/local/lib/R/site-library',
    '/usr/lib/R/site-library',
    '/usr/lib/R/library'
  )
)

options(width = 100)

setwd('/data')

set.seed(1234567)

##----------------------------------------------------------------------------##
## Load libraries.
##----------------------------------------------------------------------------##
library('Seurat', lib.loc = '/other_R_packages')
library('cerebroPrepare')
library('tidyverse')
library('monocle')

##----------------------------------------------------------------------------##
## Load count matrix.
##----------------------------------------------------------------------------##
path_to_data <- './raw_data'

feature_matrix <- Matrix::readMM(paste0(path_to_data, '/matrix.mtx.gz'))
feature_matrix <- as.matrix(feature_matrix)
feature_matrix <- as.data.frame(feature_matrix)

colnames(feature_matrix) <- read_tsv(paste0(path_to_data, '/barcodes.tsv.gz'), col_names = FALSE) %>%
  select(1) %>%
  t() %>%
  as.vector()

gene_names <- read_tsv(paste0(path_to_data, '/features.tsv.gz'), col_names = FALSE) %>%
  select(2) %>%
  t() %>%
  as.vector()

feature_matrix <- feature_matrix %>%
  mutate(gene = gene_names) %>%
  select('gene', everything()) %>%
  group_by(gene) %>%
  summarise_all(sum)

genes <- feature_matrix$gene

feature_matrix <- select(feature_matrix, -c('gene'))
feature_matrix <- as.data.frame(feature_matrix)
rownames(feature_matrix) <- genes

##----------------------------------------------------------------------------##
## Basic Seurat analysis.
##----------------------------------------------------------------------------##
seurat <- CreateSeuratObject(
  project = 'PBMC_10k_v3',
  raw.data = feature_matrix,
  min.cells = 10
)

seurat <- FilterCells(
  seurat,
  subset.names = c('nUMI', 'nGene'),
  low.thresholds = c(100, 50),
  high.thresholds = c(Inf, Inf)
)

seurat <- NormalizeData(seurat)

seurat <- FindVariableGenes(seurat)

seurat <- ScaleData(seurat, vars.to.regress = 'nUMI')

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
  reduction.type = 'pca',
  dims.use = 1:30,
  resolution = 0.5,
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
  col.name = 'cluster'
)

seurat@meta.data$res.0.7 <- NULL

##----------------------------------------------------------------------------##
## Dimensional reductions.
##----------------------------------------------------------------------------##
seurat <- RunTSNE(
  seurat,
  reduction.name = 'tSNE',
  reduction.key = 'tSNE_',
  dims.use = 1:30,
  perplexity = 30,
  do.fast = TRUE,
  seed.use = 100
)

seurat <- RunTSNE(
  seurat,
  reduction.name = 'tSNE_3D',
  reduction.key = 'tSNE3D_',
  dims.use = 1:30,
  dim.embed = 3,
  perplexity = 30,
  do.fast = TRUE,
  seed.use = 100
)

seurat <- RunUMAP(
  seurat,
  reduction.name = 'UMAP',
  reduction.key = 'UMAP_',
  dims.use = 1:30,
  do.fast = TRUE,
  seed.use = 100
)

seurat <- RunUMAP(
  seurat,
  reduction.name = 'UMAP_3D',
  reduction.key = 'UMAP3D_',
  dims.use = 1:30,
  max.dim = 3,
  do.fast = TRUE,
  seed.use = 100
)

##----------------------------------------------------------------------------##
## Randomly assign cells to three samples.
##----------------------------------------------------------------------------##
sample_dump <- c(
  rep('PBMC_5k_1', 10000),
  rep('PBMC_5k_2', 10000),
  rep('PBMC_5k_3', 10000)
)

sample <- data.frame(
  'sample' = sample(sample_dump, nrow(seurat@meta.data)),
  row.names = seurat@cell.names
)

seurat <- AddMetaData(
  seurat,
  metadata = sample,
  col.name = 'sample'
)

seurat@meta.data$sample <- factor(
  seurat@meta.data$sample,
  levels = c('PBMC_5k_1','PBMC_5k_2','PBMC_5k_3')
)

##----------------------------------------------------------------------------##
## Store meta data and parameters of analysis.
##----------------------------------------------------------------------------##
seurat@misc$experiment <- list(
  experiment_name = 'PBMC_10k',
  organism = 'hg',
  date_of_analysis = Sys.Date()
)

seurat@misc$parameters <- list(
  gene_nomenclature = 'gene_name',
  discard_genes_expressed_in_fewer_cells_than = 10,
  keep_mitochondrial_genes = TRUE,
  variables_to_regress_out = 'nUMI',
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

seurat@misc$technical_info <- list(
  'R' = capture.output(devtools::session_info())
)

##----------------------------------------------------------------------------##
## Perform further analysis with cerebroPrepare and export file for Cerebro.
##----------------------------------------------------------------------------##
seurat <- cerebroPrepare::addPercentMtRibo(
  seurat,
  organism = 'hg',
  gene_nomenclature = 'name'
)

seurat <- cerebroPrepare::getMostExpressedGenes(
  seurat,
  column_sample = 'sample',
  column_cluster = 'cluster'
)

seurat <- cerebroPrepare::getMarkerGenes(
  seurat,
  organism = 'hg',
  column_sample = 'sample',
  column_cluster = 'cluster'
)

seurat <- cerebroPrepare::getEnrichedPathways(
  seurat,
  column_sample = 'sample',
  column_cluster = 'cluster',
  adj_p_cutoff = 0.01,
  max_terms = 100
)

##----------------------------------------------------------------------------##
## Trajectory analysis with Monocle using all cells.
##----------------------------------------------------------------------------##
monocle_all_cells <- newCellDataSet(
  seurat@data,
  phenoData = new('AnnotatedDataFrame', data = seurat@meta.data),
  featureData = new('AnnotatedDataFrame', data = data.frame(gene_short_name = rownames(seurat@data), row.names = rownames(seurat@data)))
)

monocle_all_cells <- estimateSizeFactors(monocle_all_cells)
monocle_all_cells <- estimateDispersions(monocle_all_cells)
monocle_all_cells <- setOrderingFilter(monocle_all_cells, seurat@var.genes)
monocle_all_cells <- reduceDimension(monocle_all_cells, max_components = 2, method = 'DDRTree')
monocle_all_cells <- orderCells(monocle_all_cells)

seurat <- cerebroPrepare::extractMonocleTrajectory(monocle_all_cells, seurat, 'all_cells')

##----------------------------------------------------------------------------##
## Trajectory analysis with Monocle using only cells in G1 phase.
##----------------------------------------------------------------------------##
G1_cells <- which(seurat@meta.data$Phase == 'G1')

monocle_subset_of_cells <- newCellDataSet(
  seurat@data[,G1_cells],
  phenoData = new('AnnotatedDataFrame', data = seurat@meta.data[G1_cells,]),
  featureData = new('AnnotatedDataFrame', data = data.frame(gene_short_name = rownames(seurat@data), row.names = rownames(seurat@data)))
)
monocle_subset_of_cells <- estimateSizeFactors(monocle_subset_of_cells)
monocle_subset_of_cells <- estimateDispersions(monocle_subset_of_cells)
monocle_subset_of_cells <- setOrderingFilter(monocle_subset_of_cells, seurat@var.genes)
monocle_subset_of_cells <- reduceDimension(monocle_subset_of_cells, max_components = 2, method = 'DDRTree')
monocle_subset_of_cells <- orderCells(monocle_subset_of_cells)

seurat <- cerebroPrepare::extractMonocleTrajectory(monocle_subset_of_cells, seurat, 'subset_of_cells')

##----------------------------------------------------------------------------##
## Export data for Cerebro.
##----------------------------------------------------------------------------##
cerebroPrepare::exportFromSeurat(
  seurat,
  experiment_name = 'PBMC_10k',
  file = paste0('Seurat_v2/cerebro_PBMC_10k_', Sys.Date(), '.crb'),
  organism = 'hg',
  column_cell_cycle_seurat = 'Phase'
)

##----------------------------------------------------------------------------##
## Save Seurat object.
##----------------------------------------------------------------------------##
saveRDS(seurat, 'Seurat_v2/seurat.rds')

