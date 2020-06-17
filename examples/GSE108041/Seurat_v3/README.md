# Seurat v3 workflow for `GSE108041` data set

Here, we analyze a `GSE108041` data set ("Extreme heterogeneity of influenza virus infection in single cells", Russell *et al.*, eLIFE (2018), [DOI](https://doi.org/10.7554/eLife.32303), [GEO submission](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108041)) using [Seurat](https://satijalab.org/seurat/) framework, following the basic Seurat workflow.

## Preparation

Before starting, we clone the Cerebro repository (or manually download it) because it contains the raw data of our example data set.
One (optional) step of our analysis will require us to provide some gene sets in a `GMT` file.
We manually download the `c2.all.v7.0.symbols.gmt` file from [MSigDB](http://software.broadinstitute.org/gsea/downloads.jsp#msigdb) and put it in our current working directory.
Then, we pull the Docker image from the Docker Hub, convert it to Singularity, and start an R session inside.

```sh
git clone https://github.com/romanhaa/Cerebro
cd Cerebro/examples/GSE108041
# download GMT file (if you want) and place it inside this folder
singularity build <path_to>/cerebro_v1.1.simg docker://romanhaa/cerebro:v1.1
singularity exec --bind ./:/data <path_to>/cerebro_v1.1.simg R
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
library('cerebroApp')
```

## Load transcript counts

For each of the four samples we load the transcript count matrix (`.h5` format), add a tag for the sample of origin to the cellular barcodes, and then merge the transcript counts together.

```r
sample_uninfected <- Read10X_h5('raw_data/GSM2888370_Uninfected.h5') %>%
  as.matrix() %>%
  as.data.frame()
colnames(sample_uninfected) <- colnames(sample_uninfected) %>%
  gsub(., pattern = '-1', replacement = '-uninfected')
sample_uninfected <- sample_uninfected %>%
  dplyr::mutate(gene = rownames(.))

sample_6h <- Read10X_h5('raw_data/GSM2888371_6h.h5') %>%
  as.matrix() %>%
  as.data.frame()
colnames(sample_6h) <- colnames(sample_6h) %>%
  gsub(., pattern = '-1', replacement = '-6h')
sample_6h <- sample_6h %>%
  dplyr::mutate(gene = rownames(.))

sample_8h <- Read10X_h5('raw_data/GSM2888372_8h.h5') %>%
  as.matrix() %>%
  as.data.frame()
colnames(sample_8h) <- colnames(sample_8h) %>%
  gsub(., pattern = '-1', replacement = '-8h')
sample_8h <- sample_8h %>%
  dplyr::mutate(gene = rownames(.))

sample_10h <- Read10X_h5('raw_data/GSM2888374_10h.h5') %>%
  as.matrix() %>%
  as.data.frame()
colnames(sample_10h) <- colnames(sample_10h) %>%
  gsub(., pattern = '-1', replacement = '-10h')
sample_10h <- sample_10h %>%
  dplyr::mutate(gene = rownames(.))

feature_matrix <- dplyr::full_join(sample_uninfected, sample_6h, by = 'gene') %>%
  dplyr::full_join(sample_8h, , by = 'gene') %>%
  dplyr::full_join(sample_10h, , by = 'gene')
rownames(feature_matrix) <- feature_matrix$gene
feature_matrix <- dplyr::select(feature_matrix, -gene)
```

## Pre-processing with Seurat

With the merged transcript count matrix ready, we create a Seurat object, add sample info to meta data, and remove cells with less fewer `100` transcripts or `50` expressed genes.
Then, we follow the standard Seurat workflow, including...

* normalization,
* identifying highly variably genes,
* scaling expression values and regressing out the number of transcripts per cell,
* perform principal component analysis (PCA),
* find neighbors and clusters.

Furthermore, we build a cluster tree that represents the similarity between clusters and create a dedicated `cluster` column in the meta data.

```r
seurat <- CreateSeuratObject(
  project = 'GSE108041',
  counts = feature_matrix,
  min.cells = 10
)
seurat@meta.data$sample <- rownames(seurat@meta.data) %>%
  strsplit('-') %>%
  vapply(`[`, 2, FUN.VALUE = character(1)) %>%
  factor(levels = c('uninfected','6h','8h','10h'))
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

## Cell cycle analysis

We also perform cell cycle analysis using the `CellCycleScoring` built into Seurat.
The S and G2M phase-specific gene lists are stored in the Seurat object so we have access to these lists in Cerebro.

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

## Meta data

In order to later be able to understand how we did the analysis, we add some meta data to the `misc` slot of the Seurat object.

```r
seurat@misc$experiment <- list(
  experiment_name = 'GSE108041',
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

Using the functions provided by cerebroApp, we check the percentage of mitochondrial and ribosomal genes and, for every sample and cluster, we...

* get the 100 most expressed genes,
* identify marker genes (with the `FindAllMarkers` of Seurat),
* get enriched pathways in marker lists (using Enrichr),
* and perform gene set enrichment analysis (using GSVA).

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

## Trajectory analysis with Monocle (optional)

### All cells

Next, we perform trajectory analysis of all cells with Monocle using the previously identified highly variable genes.
We extract the trajectory from the generated Monocle object with the `extractMonocleTrajectory()` function of cerebroApp and attach it to our Seurat object.

```r
monocle_all_cells <- newCellDataSet(
  seurat@assays$RNA@data,
  phenoData = new('AnnotatedDataFrame', data = seurat@meta.data),
  featureData = new('AnnotatedDataFrame', data = data.frame(
    gene_short_name = rownames(seurat@assays$RNA@counts),
    row.names = rownames(seurat@assays$RNA@counts))
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

Then, we do the same procedure again, however this time only with a subset of cells (those which are in G1 phase of the cell cycle).

```r
G1_cells <- which(seurat@meta.data$Phase == 'G1')

monocle_subset_of_cells <- newCellDataSet(
  seurat@assays$RNA@data[,G1_cells],
  phenoData = new('AnnotatedDataFrame', data = seurat@meta.data[G1_cells,]),
  featureData = new('AnnotatedDataFrame', data = data.frame(
    gene_short_name = rownames(seurat@assays$RNA@counts),
    row.names = rownames(seurat@assays$RNA@counts))
  )
)

monocle_subset_of_cells <- estimateSizeFactors(monocle_subset_of_cells)
monocle_subset_of_cells <- estimateDispersions(monocle_subset_of_cells)
monocle_subset_of_cells <- setOrderingFilter(monocle_subset_of_cells, seurat@assays$RNA@var.features)
monocle_subset_of_cells <- reduceDimension(monocle_subset_of_cells, max_components = 2, method = 'DDRTree')
monocle_subset_of_cells <- orderCells(monocle_subset_of_cells)

seurat <- cerebroApp::extractMonocleTrajectory(monocle_subset_of_cells, seurat, 'subset_of_cells')
```

## Export to Cerebro format

Finally, we use the `exportFromSeurat()` function of cerebroApp to export our Seurat object to a `.crb` file which can be loaded into Cerebro.

```r
cerebroApp::exportFromSeurat(
  seurat,
  experiment_name = 'GSE108041',
  file = paste0('Seurat_v3/cerebro_GSE108041_', Sys.Date(), '.crb'),
  organism = 'hg',
  column_nUMI = 'nCount_RNA',
  column_nGene = 'nFeature_RNA',
  column_cell_cycle_seurat = 'Phase'
)
```

## Save Seurat object

Very last step: Save the Seurat object.

```r
saveRDS(seurat, 'Seurat_v3/seurat.rds')
```
