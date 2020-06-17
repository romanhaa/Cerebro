# Seurat v2 workflow for `pbmc_10k_v3` data set

Here, we analyze the `pbmc_10k_v3` data set using [Seurat](https://satijalab.org/seurat/) framework, following the basic [Seurat](https://satijalab.org/seurat/) workflow.

## Preparation

Before starting, we clone the Cerebro repository (or manually download it) because it contains the raw data of our example data set.
One (optional) step of our analysis will require us to provide some gene sets in a `GMT` file.
We manually download the `c2.all.v7.0.symbols.gmt` file from [MSigDB](http://software.broadinstitute.org/gsea/downloads.jsp#msigdb) and put it in our current working directory.
Then, we pull the Docker image from the Docker Hub, convert it to Singularity, and start an R session inside.

```sh
git clone https://github.com/romanhaa/Cerebro
cd Cerebro/examples/pbmc_10k_v3
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
library('Seurat', lib.loc = '/other_R_packages')
library('monocle')
library('cerebroApp')
```

## Load transcript counts

Unfortunately, the `Read10X_h5()` function of Seurat v2 has problems with the `.h5` file downloaded from the 10x Genomics website so instead we load it manually, convert it to a sparse matrix, and merge transcripts from genes with the same name.

```r
h5_data <- hdf5r::H5File$new('raw_data/filtered_feature_bc_matrix.h5', mode = 'r')

feature_matrix <- Matrix::sparseMatrix(
  i = h5_data[['matrix/indices']][],
  p = h5_data[['matrix/indptr']][],
  x = h5_data[['matrix/data']][],
  dimnames = list(
    h5_data[['matrix/features/name']][],
    h5_data[['matrix/barcodes']][]
  ),
  dims = h5_data[['matrix/shape']][],
  index1 = FALSE
)

genes <- rownames(feature_matrix)

feature_matrix <- feature_matrix %>%
  as.matrix() %>%
  as.data.frame() %>%
  dplyr::mutate(gene = genes) %>%
  dplyr::select(gene, dplyr::everything()) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise_all(sum) %>%
  dplyr::ungroup()

genes <- feature_matrix$gene

feature_matrix <- feature_matrix %>%
  dplyr::select(-gene) %>%
  as.matrix() %>%
  as('sparseMatrix')

rownames(feature_matrix) <- genes
```

## Pre-processing with Seurat

With the transcript count loaded, we create a Seurat object and remove cells with less than `100` transcripts or fewer than `50` expressed genes.
Then, we follow the standard Seurat workflow, including...

* normalization,
* identifying highly variably genes,
* scaling expression values and regressing out the number of transcripts per cell,
* perform principal component analysis (PCA),
* find neighbors and clusters.

Furthermore, we build a cluster tree that represents the similarity between clusters and create a dedicated `cluster` column in the meta data.

```r
seurat <- CreateSeuratObject(
  project = 'pbmc_10k_v3',
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
seurat <- RunPCA(seurat, pcs.compute = 30, pc.genes = seurat@var.genes)
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

## Cell cycle analysis

We also perform cell cycle analysis using the `CellCycleScoring` built into Seurat.
The S and G2M phase-specific gene lists are stored in the Seurat object so we have access to these lists in Cerebro.

```r
seurat <- CellCycleScoring(
  seurat,
  g2m.genes = cc.genes$g2m.genes,
  s.genes = cc.genes$s.genes
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

## Meta data

This example data set consists of a single sample so we just add that name to the meta data.
Moreover, in order to later be able to understand how we did the analysis, we add some meta data to the `misc` slot of the Seurat object.

```r
seurat@meta.data$sample <- factor('pbmc_10k_v3', levels = 'pbmc_10k_v3')

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
  seurat@data,
  phenoData = new('AnnotatedDataFrame', data = seurat@meta.data),
  featureData = new('AnnotatedDataFrame', data = data.frame(
    gene_short_name = rownames(seurat@counts),
    row.names = rownames(seurat@counts))
  )
)

monocle_all_cells <- estimateSizeFactors(monocle_all_cells)
monocle_all_cells <- estimateDispersions(monocle_all_cells)
monocle_all_cells <- setOrderingFilter(monocle_all_cells, seurat@var.genes)
monocle_all_cells <- reduceDimension(monocle_all_cells, max_components = 2, method = 'DDRTree')
monocle_all_cells <- orderCells(monocle_all_cells)

seurat <- cerebroApp::extractMonocleTrajectory(monocle_all_cells, seurat, 'all_cells')
```

### Cells in G1 phase

Then, we do the same procedure again, however this time only with a subset of cells (those which are in G1 phase of the cell cycle).

```r
G1_cells <- which(seurat@meta.data$Phase == 'G1')

monocle_subset_of_cells <- newCellDataSet(
  seurat@data[,G1_cells],
  phenoData = new('AnnotatedDataFrame', data = seurat@meta.data[G1_cells,]),
  featureData = new('AnnotatedDataFrame', data = data.frame(
    gene_short_name = rownames(seurat@counts),
    row.names = rownames(seurat@counts))
  )
)

monocle_subset_of_cells <- estimateSizeFactors(monocle_subset_of_cells)
monocle_subset_of_cells <- estimateDispersions(monocle_subset_of_cells)
monocle_subset_of_cells <- setOrderingFilter(monocle_subset_of_cells, seurat@var.genes)
monocle_subset_of_cells <- reduceDimension(monocle_subset_of_cells, max_components = 2, method = 'DDRTree')
monocle_subset_of_cells <- orderCells(monocle_subset_of_cells)

seurat <- cerebroApp::extractMonocleTrajectory(monocle_subset_of_cells, seurat, 'subset_of_cells')
```

## Export to Cerebro format

Finally, we use the `exportFromSeurat()` function of cerebroApp to export our Seurat object to a `.crb` file which can be loaded into Cerebro.

```r
cerebroApp::exportFromSeurat(
  seurat,
  experiment_name = 'pbmc_10k_v3',
  file = paste0('Seurat_v2/cerebro_pbmc_10k_v3_', Sys.Date(), '.crb'),
  organism = 'hg',
  column_cell_cycle_seurat = 'Phase'
)
```

## Save Seurat object

Very last step: Save the Seurat object.

```r
saveRDS(seurat, 'Seurat_v2/seurat.rds')
```
