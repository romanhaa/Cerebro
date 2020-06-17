# How to generate the pre-loaded example data set for the Cerebro user interface

## Download raw data

Download raw data from Cerebro GitHub repository: <https://github.com/romanhaa/Cerebro/blob/master/examples/pbmc_10k_v3/raw_data/filtered_feature_bc_matrix.h5>

## Preparation

```r
library('dplyr')
library('Seurat')
library('monocle')
library('cerebroApp')

options(width = 100)
set.seed(1234567)
```

## Load and randomly downsample data to 501 cells

```r
feature_matrix <- Read10X_h5('filtered_feature_bc_matrix.h5')
feature_matrix <- feature_matrix[ , sample(1:ncol(feature_matrix), 501) ]
```

## Basic Seurat workflow

```r
seurat <- CreateSeuratObject(
  project = 'pbmc_10k_v3',
  counts = feature_matrix,
  min.cells = 10
)
seurat <- subset(seurat, subset = nCount_RNA > 100 & nFeature_RNA > 50)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat, vars.to.regress = 'nCount_RNA')
seurat <- RunPCA(seurat, npcs = 30, features = seurat@assays$RNA@var.features)
seurat <- FindNeighbors(seurat)
seurat <- FindClusters(seurat, resolution = 0.5)
seurat <- BuildClusterTree(
  seurat,
  dims = 1:30,
  reorder = TRUE,
  reorder.numeric = TRUE
)
seurat[['cluster']] <- factor(
  as.character(seurat@meta.data$tree.ident),
  levels = sort(unique(seurat@meta.data$tree.ident))
)
seurat@meta.data$seurat_clusters <- NULL
seurat@meta.data$RNA_snn_res.0.5 <- NULL
seurat@meta.data$tree.ident <- NULL
```

## Randomly assign cells to samples

```r
sample_info <- rep(NA, 501)
sample_info[1:166] <- 'pbmc_10k_v3_rep1'
sample_info[167:334] <- 'pbmc_10k_v3_rep2'
sample_info[335:501] <- 'pbmc_10k_v3_rep3'
sample_info <- factor(sample_info, levels = c('pbmc_10k_v3_rep1','pbmc_10k_v3_rep2','pbmc_10k_v3_rep3'))
seurat@meta.data$sample <- sample_info
```

## Cell cycle analysis

```r
seurat <- CellCycleScoring(
seurat,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)

seurat@misc$gene_lists$G2M_phase_genes <- cc.genes$g2m.genes
seurat@misc$gene_lists$S_phase_genes <- cc.genes$s.genes
```

## Dimensional reduction

```r
seurat <- RunTSNE(
  seurat,
  reduction.name = 'tSNE',
  reduction.key = 'tSNE_',
  dims = 1:30,
  dim.embed = 2,
  perplexity = 30,
  seed.use = 100
)

seurat <- RunTSNE(
  seurat,
  reduction.name = 'tSNE_3D',
  reduction.key = 'tSNE3D_',
  dims = 1:30,
  dim.embed = 3,
  perplexity = 30,
  seed.use = 100
)

seurat <- RunUMAP(
  seurat,
  reduction.name = 'UMAP',
  reduction.key = 'UMAP_',
  dims = 1:30,
  n.components = 2,
  seed.use = 100
)

seurat <- RunUMAP(
  seurat,
  reduction.name = 'UMAP_3D',
  reduction.key = 'UMAP3D_',
  dims = 1:30,
  n.components = 3,
  seed.use = 100
)
```

## Meta data

```r
seurat@misc$experiment <- list(
  experiment_name = 'pbmc_10k_v3',
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
  cluster_resolution = 0.5
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
```

## cerebroApp functions (optional but recommended)

```r
seurat <- cerebroApp::addPercentMtRibo(
  seurat,
  organism = 'hg',
  gene_nomenclature = 'name'
)

seurat <- cerebroApp::getMostExpressedGenes(
  seurat,
  column_sample = 'sample',
  column_cluster = 'cluster'
)

seurat <- cerebroApp::getMarkerGenes(
  seurat,
  organism = 'hg',
  column_sample = 'sample',
  column_cluster = 'cluster'
)

seurat <- cerebroApp::getEnrichedPathways(
  seurat,
  column_sample = 'sample',
  column_cluster = 'cluster',
  adj_p_cutoff = 0.01,
  max_terms = 100
)

seurat <- cerebroApp::performGeneSetEnrichmentAnalysis(
  seurat,
  GMT_file = 'c2.all.v7.0.symbols.gmt',
  column_sample = 'sample',
  column_cluster = 'cluster',
  thresh_p_val = 0.05,
  thresh_q_val = 0.1,
  parallel.sz = 1,
  verbose = FALSE
)
```

## Trajectory analysis with Monocle

### All cells

```r
monocle_all_cells <- newCellDataSet(
  seurat@assays$RNA@data,
  phenoData = new('AnnotatedDataFrame', data = seurat@meta.data),
  featureData = new('AnnotatedDataFrame', data = data.frame(
    gene_short_name = rownames(seurat@assays$RNA@data),
    row.names = rownames(seurat@assays$RNA@data))
  )
)

monocle_all_cells <- estimateSizeFactors(monocle_all_cells)
monocle_all_cells <- estimateDispersions(monocle_all_cells)
monocle_all_cells <- setOrderingFilter(monocle_all_cells, seurat@assays$RNA@var.features)
monocle_all_cells <- reduceDimension(monocle_all_cells, max_components = 2, method = 'DDRTree')
monocle_all_cells <- orderCells(monocle_all_cells)

seurat <- cerebroApp::extractMonocleTrajectory(monocle_all_cells, seurat, 'all_cells')
```

### Cells in G1 phase

```r
G1_cells <- which(seurat@meta.data$Phase == 'G1')

monocle_subset_of_cells <- newCellDataSet(
  seurat@assays$RNA@data[,G1_cells],
  phenoData = new('AnnotatedDataFrame', data = seurat@meta.data[G1_cells,]),
  featureData = new('AnnotatedDataFrame', data = data.frame(
    gene_short_name = rownames(seurat@assays$RNA@data),
    row.names = rownames(seurat@assays$RNA@data))
  )
)

monocle_subset_of_cells <- estimateSizeFactors(monocle_subset_of_cells)
monocle_subset_of_cells <- estimateDispersions(monocle_subset_of_cells)
monocle_subset_of_cells <- setOrderingFilter(monocle_subset_of_cells, seurat@assays$RNA@var.features)
monocle_subset_of_cells <- reduceDimension(monocle_subset_of_cells, max_components = 2, method = 'DDRTree')
monocle_subset_of_cells <- orderCells(monocle_subset_of_cells)

seurat <- cerebroApp::extractMonocleTrajectory(monocle_subset_of_cells, seurat, 'subset_of_cells')
```

## Randomly downsample genes to 1000

```r
seurat@assays$RNA@data <- seurat@assays$RNA@data[ sample(1:nrow(seurat@assays$RNA@data), 1000) , ]
```

## Export to Cerebro format

```r
cerebroApp::exportFromSeurat(
  seurat,
  experiment_name = 'pbmc_10k_v3',
  file = 'example.crb',
  organism = 'hg',
  column_nUMI = 'nCount_RNA',
  column_nGene = 'nFeature_RNA',
  column_cell_cycle_seurat = 'Phase'
)
```
