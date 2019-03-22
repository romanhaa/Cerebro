# Example data

As an example/test data set, we used a public scRNA-seq data set containing about 10k human PBMCs from a healthy donor (`pbmc_10k_v3`), generated using the v3 chemistry and processed using Cell Ranger 3.0.0.
The dataset is available throught the 10x Genomics website (<https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3>) and is licensed under the Creative Commons Attribution license (<https://creativecommons.org/licenses/by/4.0/>).

Briefly, in the example we followed the basic [Seurat](https://satijalab.org/seurat/) workflow:

* Load the feature matrix and create a Seurat object (`CreateSeuratObject()`).
* Filter cells based on the number of transcripts and expressed genes (`FilterCells()`).
* Normalize the transcript counts using the `LogNormalize` method and scaled each cell to contain 10,000 transcripts (`NormalizeData()`).
* Identify variable genes (`FindVariableGenes()`).
* Scale the expression matrix and regressing out the number of transcripts (`ScaleData()`).
* Perform cell cycle analysis (`CellCycleScoring()`).
* Perform principal component analysis (`RunPCA`).
* Identify clusters and build a cluster tree (`FindClusters()`, `BuildClusterTree()`).
* Perform 2D and 3D dimensional reduction: tSNE, UMAP, diffusion map, PHATE (`RunTSNE`, `RunUMAP`, `RunDiffusion`, `RunPHATE`).

Then, we randomly assign each cell to one of three samples to simulate a dataset with multiple samples and subsequently use cerebroPrepare to:

* Calculate the percent of mitochondrial and ribosomal gene expression (`addPercentMtRibo()`).
* Get the most expressed genes in each sample and cluster (`getMostExpressedGenes()`).
* Get marker genes for each sample and cluster (`getMarkerGenes()`).
* Perform pathway enrichment analysis using the marker genes of each sample and cluster (`getEnrichedPathways()`).
* Export a `.crb` file that can be loaded into Cerebro (`exportFromSeurat()`).

To test Cerebro, download the `.crb` file and load it into Cerebro.

## How to reproduce

I suggest to run the example R script in a container using [Singularity](https://singularity.lbl.gov/) (here I used Singularity 2.6.0).

```
git clone romanhaa/Cerebro
cd Cerebro/test_data
singularity build cerebro-manuscript.simg docker://romanhaa/cerebro-manuscript:2019-03-01
singularity exec --bind ./:/data cerebroManuscript.simg Rscript /data/Cerebro_example.R
```

This will reproduce the `.crb` file used as an example.
