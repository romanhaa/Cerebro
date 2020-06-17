# How to generate the Seurat PBMC example data set used in the examples parts of cerebroApp function

```r
library(Seurat)
library(cerebroApp)

pbmc_counts <- read.table(
  file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
  as.is = TRUE
)

pbmc <- CreateSeuratObject(counts = pbmc_counts)

pbmc@meta.data$sample <- factor('A', levels = 'A')

pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
pbmc <- FindNeighbors(object = pbmc)
pbmc <- FindClusters(object = pbmc)

pbmc <- RunUMAP(
  pbmc,
  reduction.name = 'UMAP',
  reduction.key = 'UMAP_',
  dims = 1:30,
  n.components = 2,
  seed.use = 100,
  verbose = FALSE
)

pbmc <- getMarkerGenes(
  pbmc,
  organism = 'hg',
  column_sample = 'sample',
  column_cluster = 'seurat_clusters'
)

saveRDS(pbmc, 'inst/extdata/v1.2/seurat_pbmc.rds')
```
