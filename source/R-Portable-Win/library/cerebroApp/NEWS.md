# cerebroApp 1.2.2

## Fixes

- The title in the browser tab now correctly says "Cerebro" instead of containing some HTML code.
- Cluster trees should now be displayed correctly.
- `getEnrichedPathways()` no longer results in an error when marker genes are present but no database returns any enriched pathways, e.g. because there are too few marker genes. Thanks to @turkeyri for pointing it out and suggesting a solution!

# cerebroApp 1.2.1

## New features

- It is now possible to select cells in the dimensional reduction plots ('Overview', 'Gene expression', and 'Gene set expression' tabs) and retrieve additional info for them. For example, users can get tables of meta data or expression values and save them as a file for further analysis. Also, gene expression can be shown in the selected vs. non-selected cells.

## Minor changes

- Scales for expression levels by sample and cluster in "Gene expression" and "Gene set expression" tabs are now set to be from 0 to 1.2 times the highest value. This is to limit the violin plots which cannot be trimmed to the actual data range and will extend beyond, giving a false impression of negative values existing in the data.
- Hover info in expression by gene plot in "Gene expression" and "Gene set expression" tabs now show both the gene name and the mean expression value instead of just the gene name.

# cerebroApp 1.2.0

## New features

- New button for composition plots (e.g. samples by clusters or cell cycle) that allows to choose whether to scale by actual cell count (default) or percentage.
- New button for composition plots that allows to show/hide the respective table of numbers behind them.
- New tab "Color management": Users can now change the color assigned to each sample/cluster.
- "Gene expression" and "Gene set expression" panels: Users can now pick from a set of color scales and adjust the color range.
- The gene selection box in the "Gene expression" panel will now allow to view available genes and select them by clicking. It is not necessary anymore to hit Enter or Space to update the plot, this will be done automatically after providing new input.
- It is now possible to export assays other than `RNA` through the `assay` parameter in relevant functions.
- Launch old Cerebro interfaces through `version` parameter in `launchCerebro()`.
- We added a vignette which explains how to use cerebroApp and its functions.

## Minor changes

- Add citation info.
- Composition tables (e.g. samples by clusters or cell cycle) are now calculated in the Shiny app rather than being expected to be present in the `.crb file.
- Fix log message in `exportFromSeurat()` when extracting trajectories.
- The gene set selection box in the "Gene set expression" tab will not crash anymore when typing a sequence of letters that doesn't match any gene set names.
- Remove dependency on pre-assigned colors in the `.crb` file. If no colors have been assigned to samples and clusters when loading a data set, they will be assigned then.
- Update examples of functions and include mini-Seurat object and example gene set (GMT file) to run the examples.
- Modify pre-loaded data set in Cerebro interface to contain more data.
- When attempting to download genes in GO term "cell surface" in the `getMarkerGenes()` function, it tries at max. 3 times to contact the biomaRt server and continues without if all attempts failed. Sometimes the server does not respond which gave an error in previous versions of the function.
- Plenty of changes to meet Bioconductor guidelines (character count per line, replace `.` in dplyr pipes with `rlang::.data`, etc.).
- Reduce package size by compressing reference files, e.g. gene name/ID conversion tables.

# cerebroApp 1.1.0

- Release along with manuscript revision.

## New features

- New function `extractMonocleTrajectory()`: Users can extract data from trajectories calculated with Monocle v2.
- New tab "Trajectories": Allows visualization of trajectories calculated with Monocle v2.

# cerebroApp 1.0.0

- Public release along with manuscript submission to bioRxiv.
