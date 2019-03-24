# Notes and changes

* Disable zoom in scatterD3.
  * In `d3-4.13.0.min.js`, change the `wheel.zoom` to `null`.
* Change color scale in scatterD3 and export functions to `YlGnBu`.
  * Required adding `d3-scale-chromatic.v.0.3.min.js` to scatterD3 and change object name inside.
  * Eventually consider making this a variable to choose by user, e.g. blue vs red scale.
* Make own box function so I don't have to repeat parameters all the time?
  * Same for action button.
  * Also modalDialog.
* Remove `as.vector()` for gene expression, not necessary.
* Check/update export scripts.
* Move to z-score instead of log-counts?
    * Not possible because occupy much more data.
* Use `formattable`.
* Make only RDS files selectable.
* Change text in some panels.
* Make input in gene expression panel case-independent.
* Make refresh shortcut local (instead of global).
* Make sample names be displayed instead of the long names (top expressed genes and marker genes).
* Border around plots.
* Show names in sample-by-cluster and cluster-by-sample plots and cell cycle in respective plots.
* Allow user to choose plotting order.
* Reorder some user inputs to save space.
* Show cluster tree if available.
* Hide some columns from tables in enriched pathways to have a cleaner view. They can still be shown manually.
* Use a box layout across the app.
* Use more attractive colors (from <flatuicolors.com>).
* Add info popup boxes with additional explanation for each/most panel(s).
* Reduce empty space (padding around boxes) across the app.
* Remove Z-score from marker gene tables since combined score should be more informative.
* Add color bar representation through formattable in marker gene tables.
* Remove legend in gene expression and gene set expression scatter plots because scale will still be visible in box plots.
* Add checks for presence of tables of top expressed genes and marker genes.
* Hide gene column in enriched pathway tab by default to save space.
* Update R version for Windows to 3.5.1.
* Add export buttons for dimensional reductions using ggplot instead of the scatterD3 function which crashes with too many cells.
* Manually modify code for scatterD3 legend to prevent sorting. A more consistent solution should be found however.
* If marker genes haven't already been checked in the pipeline, this will be done in the browser for murine and human samples.

# How to uninstall all non-base R packages

Taken from: https://www.r-bloggers.com/how-to-remove-all-user-installed-packages-in-r/

```r
# create a list of all installed packages
ip <- as.data.frame(installed.packages())

# we don't want to remove base or recommended packages either\
ip <- ip[!(ip[,"Priority"] %in% c("base", "recommended")),]

# determine the library where the packages are installed
path.lib <- unique(ip$LibPath)

# create a vector with all the names of the packages you want to remove
pkgs.to.remove <- ip[,1]

# remove the packages
sapply(pkgs.to.remove, remove.packages, lib = path.lib)
```

Then, re-install the libraries required by Cerebro.

```r
install.packages("BiocManager")
BiocManager::install("shiny")
setwd("../Cerebro/source/")
library("shiny")
runApp("app.R")
```
