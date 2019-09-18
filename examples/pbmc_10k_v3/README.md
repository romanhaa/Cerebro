# Example data

As an example data set, we used a public scRNA-seq data set containing about 10k human PBMCs from a healthy donor (`pbmc_10k_v3`), generated using the v3 chemistry and processed using Cell Ranger 3.0.0.
The data set is available through the [10x Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3) and is licensed under the [Creative Commons Attribution 4.0 license](https://creativecommons.org/licenses/by/4.0/).

Briefly, in the example we followed the basic [Seurat](https://satijalab.org/seurat/) workflow (both Seurat v2 and Seurat v3 are supported):

* Load the feature matrix and create a Seurat object (`CreateSeuratObject()`).
* Filter cells based on the number of transcripts and expressed genes (`FilterCells()`).
* Normalize the transcript counts using the `LogNormalize` method and scaled each cell to contain 10,000 transcripts (`NormalizeData()`).
* Identify variable genes (`FindVariableGenes()`).
* Scale the expression matrix and regressing out the number of transcripts (`ScaleData()`).
* Perform cell cycle analysis (`CellCycleScoring()`).
* Perform principal component analysis (`RunPCA`).
* Identify clusters and build a cluster tree (`FindClusters()`, `BuildClusterTree()`).
* Perform 2D and 3D dimensional reduction: tSNE, UMAP (`RunTSNE`, `RunUMAP`).

Then, we add some meta data, randomly assign each cell to one of three samples to simulate a dataset with multiple samples and subsequently use cerebroPrepare to:

* Calculate the percent of mitochondrial and ribosomal gene expression (`addPercentMtRibo()`).
* Get the most expressed genes in each sample and cluster (`getMostExpressedGenes()`).
* Get marker genes for each sample and cluster (`getMarkerGenes()`).
* Perform pathway enrichment analysis using the marker genes of each sample and cluster (`getEnrichedPathways()`).
* Perform gene set enrichment analysis for each sample and cluster (`getEnrichedPathways()`).

Next, we calculate trajectories of (1) all cells and (2) a subset of cells (those in G1 phase) using Monocle v2 and the variable features identified by Seurat.
We extract these trajectories from the respective Monocle objects and add them to our Seurat object through the `extractMonocleTrajectory()` function.

Lastly, from the Seurat object we export a Cerebro file (`.crb` extension) that can be loaded into Cerebro (`exportFromSeurat()`).

To test Cerebro, download the `.crb` file from either [Seurat v2](Seurat_v2) or [Seurat v3](Seurat_v3) and load it into Cerebro.

## How to reproduce

The example data sets were generated in a container using [Singularity](https://singularity.lbl.gov/) (here I used Singularity 2.6.0).
The workflows for Seurat v2 and Seurat v3 are conceptually identical with some differences due to changes in the Seurat package.

### Seurat v2

First, we clone the Cerebro repository (or manually download it) because it contains the raw data of our example data set.
One (optional) step of our analysis will require us to provide some gene sets in a GMT file.
We manually download the `c2.all.v7.0.symbols.gmt` file from [MSigDB](http://software.broadinstitute.org/gsea/downloads.jsp#msigdb) and put it in our current working directory.
Then, we pull the Docker image from the Docker Hub, convert it to Singularity, and start an R session inside.

```sh
git clone https://github.com/romanhaa/Cerebro
cd Cerebro/examples/pbmc_10k_v3
# download GMT file (if you want) and place it inside this folder
singularity build <path_to>/cerebro-example_2019-09-18.simg docker://romanhaa/cerebro-example:2019-09-18
singularity exec --bind ./:/data <path_to>/cerebro-example_2019-09-18.simg R
```

Then, we set the console width to `100`, change the working directory, and set the seed.

```r
options(width = 100)
setwd('/data')
set.seed(1234567)
```

To do our analysis, we load the following libraries.

```r
library('dplyr')
library('Seurat', lib.loc = '/other_R_packages')
library('monocle')
library('cerebroPrepare')
```

Next, we load the raw data and add gene names / cell barcodes.

```r
path_to_data <- './raw_data'

feature_matrix <- Matrix::readMM(paste0(path_to_data, '/matrix.mtx.gz'))
feature_matrix <- as.matrix(feature_matrix)
feature_matrix <- as.data.frame(feature_matrix)

colnames(feature_matrix) <- readr::read_tsv(paste0(path_to_data, '/barcodes.tsv.gz'), col_names = FALSE) %>%
  dplyr::select(1) %>%
  t() %>%
  as.vector()

gene_names <- readr::read_tsv(paste0(path_to_data, '/features.tsv.gz'), col_names = FALSE) %>%
  dplyr::select(2) %>%
  t() %>%
  as.vector()

feature_matrix <- feature_matrix %>%
  dplyr::mutate(gene = gene_names) %>%
  dplyr::select('gene', everything()) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise_all(sum)

genes <- feature_matrix$gene

feature_matrix <- dplyr::select(feature_matrix, -c('gene'))
feature_matrix <- as.data.frame(feature_matrix)
rownames(feature_matrix) <- genes
```

With the transcript count loaded, we create a Seurat object, remove cells with less than `100` transcripts or fewer than `50` expressed genes.
Then, we follow the standard Seurat workflow, including normalization, identifying highly variably genes, scaling expression values and regressing out the number of transcripts per cell, perform principal component analysis (PCA), find neighbors and clusters.
Furthermore, we build a cluster tree that represents the similarity between clusters and create a dedicated `cluster` column in the meta data.

```r
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
seurat@meta.data$res.0.5 <- NULL
seurat@meta.data$tree.ident <- NULL
```

We also perform cell cycle analysis using the `CellCycleScoring` built into Seurat.

```r
seurat <- CellCycleScoring(
  seurat,
  g2m.genes = cc.genes$g2m.genes,
  s.genes = cc.genes$s.genes
)
```

Next, we generate 4 dimensional reductions: tSNE, tSNE (3D), UMAP, UMAP (3D).

```r
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
```

This example data set consists of a single sample.
To highlight the functionality of Cerebro when working with a multi-sample data set, we the cells of clusters 1-5 to `sample_A`, those in clusters 6-10 to `sample_B`, and those of clusters 11-16 to `sample_C`.

```r
meta_sample <- seurat@meta.data$cluster %>% as.character()
meta_sample[which(meta_sample %in% c('1','2','3','4','5'))] <- 'sample_A'
meta_sample[which(meta_sample %in% c('6','7','8','9','10'))] <- 'sample_B'
meta_sample[which(meta_sample %in% c('11','12','13','14','15','16'))] <- 'sample_C'
seurat@meta.data$sample <- factor(meta_sample, levels = c('sample_A','sample_B','sample_C'))
```

In order to later be able to understand how we did the analysis, we add some meta data to the `misc` slot of the Seurat object.

```r
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

Using the functions provided by cerebroPrepare, we check the percentage of mitochondrial and ribosomal genes and, for every sample and cluster, we...

* get the 100 most expressed genes,
* identify marker genes (with the `FindAllMarkers` of Seurat),
* get enriched pathways in marker lists (using Enrichr),
* and perform gene set enrichment analysis (using GSVA).

```r
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

seurat <- cerebroPrepare::performGeneSetEnrichmentAnalysis(
  seurat,
  gene_sets = 'c2.all.v7.0.symbols.gmt',
  column_sample = 'sample',
  column_cluster = 'cluster',
  thresh_p_val = 0.05,
  thresh_q_val = 0.1,
  verbose = FALSE
)
```

Next, we perform trajectory analysis of all cells with Monocle using the previously identified highly variable genes.
We extract the trajectory from the generated Monocle object with the `extractMonocleTrajectory()` function of cerebroPrepare and attach it to our Seurat object.

```r
monocle_all_cells <- newCellDataSet(
  seurat@data,
  phenoData = new('AnnotatedDataFrame', data = seurat@meta.data),
  featureData = new('AnnotatedDataFrame', data = data.frame(
    gene_short_name = rownames(seurat@data),
    row.names = rownames(seurat@data))
  )
)

monocle_all_cells <- estimateSizeFactors(monocle_all_cells)
monocle_all_cells <- estimateDispersions(monocle_all_cells)
monocle_all_cells <- setOrderingFilter(monocle_all_cells, seurat@var.genes)
monocle_all_cells <- reduceDimension(monocle_all_cells, max_components = 2, method = 'DDRTree')
monocle_all_cells <- orderCells(monocle_all_cells)

seurat <- cerebroPrepare::extractMonocleTrajectory(monocle_all_cells, seurat, 'all_cells')
```

Then, we do the same procedure again, however this time only with a subset of cells (those which are in G1 phase of the cell cycle).

```r
G1_cells <- which(seurat@meta.data$Phase == 'G1')

monocle_subset_of_cells <- newCellDataSet(
  seurat@data[,G1_cells],
  phenoData = new('AnnotatedDataFrame', data = seurat@meta.data[G1_cells,]),
  featureData = new('AnnotatedDataFrame', data = data.frame(
    gene_short_name = rownames(seurat@data),
    row.names = rownames(seurat@data))
  )
)

monocle_subset_of_cells <- estimateSizeFactors(monocle_subset_of_cells)
monocle_subset_of_cells <- estimateDispersions(monocle_subset_of_cells)
monocle_subset_of_cells <- setOrderingFilter(monocle_subset_of_cells, seurat@var.genes)
monocle_subset_of_cells <- reduceDimension(monocle_subset_of_cells, max_components = 2, method = 'DDRTree')
monocle_subset_of_cells <- orderCells(monocle_subset_of_cells)

seurat <- cerebroPrepare::extractMonocleTrajectory(monocle_subset_of_cells, seurat, 'subset_of_cells')
```

Finally, we use the `exportFromSeurat()` function of cerebroPrepare to export our Seurat object to a `.crb` file which can be loaded into Cerebro.

```r
cerebroPrepare::exportFromSeurat(
  seurat,
  experiment_name = 'PBMC_10k',
  file = paste0('Seurat_v2/cerebro_PBMC_10k_', Sys.Date(), '.crb'),
  organism = 'hg',
  column_cell_cycle_seurat = 'Phase'
)
```

Very last step: Save the Seurat object.

```r
saveRDS(seurat, 'Seurat_v2/seurat.rds')
```

### Seurat v3

First, we clone the Cerebro repository (or manually download it) because it contains the raw data of our example data set.
One (optional) step of our analysis will require us to provide some gene sets in a GMT file.
We manually download the `c2.all.v7.0.symbols.gmt` file from [MSigDB](http://software.broadinstitute.org/gsea/downloads.jsp#msigdb) and put it in our current working directory.
Then, we pull the Docker image from the Docker Hub, convert it to Singularity, and start an R session inside.

```sh
git clone https://github.com/romanhaa/Cerebro
cd Cerebro/examples/pbmc_10k_v3
# download GMT file (if you want) and place it inside this folder
singularity build <path_to>/cerebro-example_2019-09-18.simg docker://romanhaa/cerebro-example:2019-09-18
singularity exec --bind ./:/data <path_to>/cerebro-example_2019-09-18.simg R
```

Then, we set the console width to `100`, change the working directory, and set the seed.

```r
options(width = 100)
setwd('/data')
set.seed(1234567)
```

To do our analysis, we load the following libraries.

```r
library('dplyr')
library('Seurat')
library('monocle')
library('cerebroPrepare')
```

Next, we load the raw data and add gene names / cell barcodes.

```r
path_to_data <- './raw_data'

feature_matrix <- Matrix::readMM(paste0(path_to_data, '/matrix.mtx.gz'))
feature_matrix <- as.matrix(feature_matrix)
feature_matrix <- as.data.frame(feature_matrix)

colnames(feature_matrix) <- readr::read_tsv(paste0(path_to_data, '/barcodes.tsv.gz'), col_names = FALSE) %>%
  dplyr::select(1) %>%
  t() %>%
  as.vector()

gene_names <- readr::read_tsv(paste0(path_to_data, '/features.tsv.gz'), col_names = FALSE) %>%
  dplyr::select(2) %>%
  t() %>%
  as.vector()

feature_matrix <- feature_matrix %>%
  dplyr::mutate(gene = gene_names) %>%
  dplyr::select('gene', everything()) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise_all(sum)

genes <- feature_matrix$gene

feature_matrix <- dplyr::select(feature_matrix, -c('gene'))
feature_matrix <- as.data.frame(feature_matrix)
rownames(feature_matrix) <- genes
```

With the transcript count loaded, we create a Seurat object, remove cells with less than `100` transcripts or fewer than `50` expressed genes.
Then, we follow the standard Seurat workflow, including normalization, identifying highly variably genes, scaling expression values and regressing out the number of transcripts per cell, perform principal component analysis (PCA), find neighbors and clusters.
Furthermore, we build a cluster tree that represents the similarity between clusters and create a dedicated `cluster` column in the meta data.

```r
seurat <- CreateSeuratObject(
  project = 'PBMC_10k_v3',
  counts = feature_matrix,
  min.cells = 10
)
seurat <- subset(seurat, subset = nCount_RNA > 100 & nFeature_RNA > 50)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat, vars.to.regress = 'nCount_RNA')
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
seurat[['cluster']] <- factor(
  as.character(seurat@meta.data$tree.ident),
  levels = sort(unique(seurat@meta.data$tree.ident))
)
seurat@meta.data$seurat_clusters <- NULL
seurat@meta.data$RNA_snn_res.0.5 <- NULL
seurat@meta.data$tree.ident <- NULL
```

We also perform cell cycle analysis using the `CellCycleScoring` built into Seurat.

```r
seurat <- CellCycleScoring(
  seurat,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)
```

Next, we generate 4 dimensional reductions: tSNE, tSNE (3D), UMAP, UMAP (3D).

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

This example data set consists of a single sample.
To highlight the functionality of Cerebro when working with a multi-sample data set, we the cells of clusters 1-6 to `sample_A`, those in clusters 7-12 to `sample_B`, and those of clusters 13-18 to `sample_C`.

```r
meta_sample <- seurat@meta.data$cluster %>% as.character()
meta_sample[which(meta_sample %in% c('1','2','3','4','5','6'))] <- 'sample_A'
meta_sample[which(meta_sample %in% c('7','8','9','10','11','12'))] <- 'sample_B'
meta_sample[which(meta_sample %in% c('13','14','15','16','17','18'))] <- 'sample_C'
seurat@meta.data$sample <- factor(meta_sample, levels = c('sample_A','sample_B','sample_C'))
```

In order to later be able to understand how we did the analysis, we add some meta data to the `misc` slot of the Seurat object.

```r
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

Using the functions provided by cerebroPrepare, we check the percentage of mitochondrial and ribosomal genes and, for every sample and cluster, we...

* get the 100 most expressed genes,
* identify marker genes (with the `FindAllMarkers` of Seurat),
* get enriched pathways in marker lists (using Enrichr),
* and perform gene set enrichment analysis (using GSVA).

```r
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

seurat <- cerebroPrepare::performGeneSetEnrichmentAnalysis(
  seurat,
  gene_sets = 'c2.all.v7.0.symbols.gmt',
  column_sample = 'sample',
  column_cluster = 'cluster',
  thresh_p_val = 0.05,
  thresh_q_val = 0.1,
  verbose = FALSE
)
```

Next, we perform trajectory analysis of all cells with Monocle using the previously identified highly variable genes.
We extract the trajectory from the generated Monocle object with the `extractMonocleTrajectory()` function of cerebroPrepare and attach it to our Seurat object.

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

seurat <- cerebroPrepare::extractMonocleTrajectory(monocle_all_cells, seurat, 'all_cells')
```

Then, we do the same procedure again, however this time only with a subset of cells (those which are in G1 phase of the cell cycle).

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

seurat <- cerebroPrepare::extractMonocleTrajectory(monocle_subset_of_cells, seurat, 'subset_of_cells')
```

Finally, we use the `exportFromSeurat()` function of cerebroPrepare to export our Seurat object to a `.crb` file which can be loaded into Cerebro.

```r
cerebroPrepare::exportFromSeurat(
  seurat,
  experiment_name = 'PBMC_10k',
  file = paste0('Seurat_v3/cerebro_PBMC_10k_', Sys.Date(), '.crb'),
  organism = 'hg',
  column_nUMI = 'nCount_RNA',
  column_nGene = 'nFeature_RNA',
  column_cell_cycle_seurat = 'Phase'
)
```

Very last step: Save the Seurat object.

```r
saveRDS(seurat, 'Seurat_v3/seurat.rds')
```
