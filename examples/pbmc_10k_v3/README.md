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

Next, we calculate trajectories of (1) all cells and (2) a subset of cells (those in G1 phase) using Monocle v2 and the variable features identified by Seurat.
We extract these trajectories from the respective Monocle objects and add them to our Seurat object through the `extractMonocleTrajectory()` function.

Lastly, from the Seurat object we export a Cerebro file (`.crb` extension) that can be loaded into Cerebro (`exportFromSeurat()`).

To test Cerebro, download the `.crb` file from either [Seurat v2](Seurat_v2) or [Seurat v3](Seurat_v3) and load it into Cerebro.

## How to reproduce

I suggest to run the example R script in a container using [Singularity](https://singularity.lbl.gov/) (here I used Singularity 2.6.0).

```sh
git clone https://github.com/romanhaa/Cerebro
cd Cerebro/examples/pbmc_10k_v3
singularity build cerebro-example_2019-08-20.simg docker://romanhaa/cerebro-example:2019-08-20
singularity exec --bind ./:/data cerebro-example.simg Rscript /data/Seurat_v2/Cerebro_example.R
singularity exec --bind ./:/data cerebro-example.simg Rscript /data/Seurat_v3/Cerebro_example.R
```

This will reproduce the `.crb` files used as an example.
