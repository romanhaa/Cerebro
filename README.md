# Cerebro

![Screenshot Cerebro: overview panel](screenshots/Screenshot_2019-03-01_at_15.22.15.png?raw=true "Cerebro: overview panel")

![Video Cerebro](screenshots/Cerebro_2019-03-12.2019-03-12_14_37_20.gif?raw=true "Cerebro")

Cerebro (cell report browser) is a standalone desktop application (currently available for macOS and Windows) which allows users to **interactively** visualize various parts of single cell transcriptomics data without requiring bioinformatic expertise.

The core is an [Shiny](https://shiny.rstudio.com/) application which is bottled into a standalone app using [Electron](https://electronjs.org/).
Therefore, it can also be run on web servers and Linux machines, requiring only R and a set of R packages.

Input data needs to be prepared using the [cerebroPrepare](https://github.com/romanhaa/cerebroPrepare) R package which was built specifically for this purpose.
It offers functionality to export a [Seurat](https://satijalab.org/seurat/) object to the correct format in a single step.
The file should be saved either with the `.rds` or `.crb` extension, indicating that internally it is an RDS object.
Furthermore, the cerebroPrepare package also provides functions to perform a set of (optional) analyses, e.g. pathway enrichment analysis based on marker gene lists of groups of cells.

The exported `.crb` file is then loaded into Cerebro and shows all available information.

Key features:

* Interactive 2D and 3D dimensional reductions.
* Sample and cluster overview panels.
* Tables of most expressed genes and marker genes for samples and clusters.
* Tables of enriched pathways in marker gene lists.
* Query gene(s) and gene sets from [MSigDB](https://http://software.broadinstitute.org/gsea/msigdb) and show their expression in dimensional reductions.
* All plots can be exported to PNG. In addition, 2D dimensional reductions can be exported to PDF.
* Tables can be downloaded in CSV or Excel format.

A basic example [Seurat](https://satijalab.org/seurat/) workflow and subsequent exporting can be found in the [`test_data`](test_data/) folder.
There you can also find the raw data and the output file that can be loaded into Cerebro.

## Interface panels in details

### Load data

Select `.rds` or `.crb` input file.
Shows number of cells, samples, clusters, as well as experiment name and organism.

### Overview

Shows 2D and 3D dimensional reductions.
Cells can be colored by meta data variables, automatically coloring the cells using a categorical or continuous scale.
Cells can be randomly down-sampled to improve performance.

### Samples

Shows sample-centric perspective of data.

* Composition of samples by cluster as table and plot.
* Distribution of number of transcripts and expressed genes by sample.
* Distribution of mitochondrial and ribosomal gene expression by sample (if it was computed with `cerebroPrepare`.
* Cell cycle by sample, either determined by the Seurat function or using Cyclone (if it was computed and assigned during exporting).

### Clusters

Shows cluster-centric perspective of data.
See info about `Samples` panel above for more details.

### Most expressed genes

If computed in `cerebroPrepare`, provides tables of most expressed genes by sample and cluster.

### Marker genes

If computed in `cerebroPrepare`, provides tables of marker genes by sample and cluster.

### Enriched pathways

If computed in `cerebroPrepare`, provides tables of enriched pathways in marker gene lists of samples and clusters.

### Gene expression

Allows to show the expression of specified genes (showing the average per cell if multiple genes) in the data set.
Calculation is triggered after pressing `SPACE` or `ENTER`.
Multiple genes must be submitted in separate lines or separated by either space, comma, semicolon.
Shows which genes are available or missing (or misspelled) in data set.
Expression levels are shown in dimensional reductions and as violin plots for every sample and cluster.
Average expression across all cells of the 50 most expressed genes (of the ones specified by the user) are shown as well to quickly spot which genes drive the color scale.

### Gene set expression

Basically the same as the gene expression panel except that it allows to select gene sets from [MSigDB](https://http://software.broadinstitute.org/gsea/msigdb) (requires internet connection).
Only available for human and mouse data.

### Gene ID conversion

Provides table that allows to convert gene IDs and names.
Includes GENCODE identifier, ENSEMBL identifier, HAVANA identifier, gene symbol and gene type.
Only available for mouse and human.
Based on GENCODE annotation version M16 (mouse) and version 27 (human).

### Analysis info

Overview of parameters that were used during the analysis, as long as they were provided.
Also shows list of mitochondrial and ribosomal genes present in the data set if computed with `cerebroPrepare`.

## Motivation

Single cell RNA-sequencing data is rich and complex.
Allowing experimental biologists to explore the results is beneficial for the iterative scientific process of performing analysis and deriving conclusions.
Cerebro provides an easy way to access the data without any bioinformatic expertise.

## Installation

Download release for your OS, unpack and run.
Currently, Cerebro is available only for macOS and Windows.

**Note:** Linux users have the option to launch the application through a Docker container or hosting the application using the `app.R` file and a local R installation.
A convenient IDE would be RStudio.

```R
install.packages("shiny")
library("shiny")
runApp("Cerebro/source/app.R")
```

## Test data

[Find documentation for a test data set here.](test_data/README.md)

## Technical notes

* Cerebro is a [Shiny](https://shiny.rstudio.com/) app that is bottled into a standalone application using [Electron](https://electronjs.org/).
* Plotting relies heavily on [`ggplot2`](https://ggplot2.tidyverse.org/) and [`plotly`](https://plot.ly/r/).
* Tables are built with [`formattable`](https://renkun-ken.github.io/formattable/).
* Access to MSigDB through [`msigdbr`](https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html)

## Building from source

### On macOS

To package Cerebro you need [Git](https://git-scm.com) and [Node.js](https://nodejs.org/en/download/) (which comes with [npm](http://npmjs.com)) installed on your computer.
Then, from the command line, run:

```bash
# clone this repository
git clone https://gitlab.com/romanhaa/Cerebro.git
# install Electron packager
npm install electron-packager --global
# go into the repository
cd Cerebro
# install dependencies
npm install
# run the app
npm start
# build the app
npm run package-mac
```

To build the Windows version under macOS it is necessary to install Wine.
I experienced problems with missing libraries of the stable version (4.0) so I recommend to use the developers version (4.4) using Homebrew:

```bash
brew tap caskroom/versions
brew update
brew install caskroom/versions/wine-devel
npm run package-win
```

### On Windows

If you're using Linux Bash for Windows, [see this guide](https://www.howtogeek.com/261575/how-to-run-graphical-linux-desktop-applications-from-windows-10s-bash-shell/) or use `node` from the command prompt.

## Troubleshooting

* If the app shows a blank/white window, press CMR+R (macOS) or CTRL+R (Windows) to refresh the page. Especially on slower machines it can happen that the interface loads before the Shiny application is launched.

## Credits

* Columbus Collaboratory layed out the basics of using Electron to create a standalone Shiny application: <https://github.com/ColumbusCollaboratory/electron-quick-start>
* Sample and cluster color palettes taken from https://flatuicolors.com/
* App icon made by [Kiranshastry](https://www.flaticon.com/authors/kiranshastry) from [https://www.flaticon.com](https://www.flaticon.com/) is licensed by [CC 3.0 BY](http://creativecommons.org/licenses/by/3.0/)

## Contribute

To report any bugs, submit patches, or request new features, please log an issue [through the issue tracker](https://github.com/romanhaa/Cerebro/issues/new). For direct inquiries, please send an email to <a href = "mailto: roman.hillje@ieo.it">roman.hillje@ieo.it</a>.

## License

Copyright (c) 2019 <COPYRIGHT HOLDER> http://www.ieo.it

[The MIT License (MIT)](LICENSE.md)

## To Do

* [ ] Create release and upload to GitHub.




