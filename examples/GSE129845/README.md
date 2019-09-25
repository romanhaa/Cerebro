# `GSE129845` data set

This data set comes from the publication "Single-Cell Transcriptomic Map of the Human and Mouse Bladders" by Yu *et al.*, J Am Soc Nephrol (2019) ([DOI](https://doi.org/10.1681/ASN.2019040335), [GEO submission](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129845)).
It contains ~13,000 cells from 3 patients.

To test Cerebro, download the `.crb` file from [Seurat v3](Seurat_v3) and load it into Cerebro.

## Workflow

The workflows of all three frameworks are conceptually the same, containing the following steps:

* Load the transcript counts.
* Filter cells based on the number of transcripts and expressed genes.
* Normalize the transcript counts and scaled each cell to contain 10,000 transcripts.
* Identify variable genes.
* Scale the expression matrix and regressing out the number of transcripts.
* Perform cell cycle analysis.
* Perform principal component analysis.
* Identify clusters and build a cluster tree.
* Perform dimensional reduction.

Then, using the functions of cerebroApp, we add some more data:

* Calculate the percent of mitochondrial and ribosomal gene expression (`addPercentMtRibo()`).
* Get the most expressed genes in each sample and cluster (`getMostExpressedGenes()`).
* Get marker genes for each sample and cluster (`getMarkerGenes()`).
* Perform pathway enrichment analysis using the marker genes of each sample and cluster (`getEnrichedPathways()`).
* Perform gene set enrichment analysis for each sample and cluster (`performGeneSetEnrichmentAnalysis()`).

Next, we calculate trajectories of (1) all cells and (2) a subset of cells (those in G1 phase) using Monocle v2 and the variable features identified by Seurat.
We extract these trajectories from the respective Monocle objects and add them to our Seurat object through the `extractMonocleTrajectory()` function.

Lastly, from the Seurat object we export a Cerebro file (`.crb` extension) that can be loaded into Cerebro (`exportFromSeurat()`).

## How to reproduce

The example data sets were generated using the official Cerebro Docker image which was built in Docker ([Docker Hub](https://cloud.docker.com/u/romanhaa/repository/docker/romanhaa/cerebro)) and imported into [Singularity](https://singularity.lbl.gov/) (here I used Singularity 2.6.0).
The workflows for Seurat v2 and Seurat v3 are conceptually identical with some differences due to changes in the Seurat package.
Details and descriptions for the workflow can be found in the respective directory [Seurat v3](Seurat_v3).
