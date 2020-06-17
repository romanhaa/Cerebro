# How to generate the example GMT file

The `GMT` file contains 20 genes present in the example Seurat PBMC data set to ensure the `performGeneSetEnrichmentAnalysis()` function will work properly.
Storing gene sets in the `GMT` format requires a name, a description (here we used a URL), followed by the genes which belong to the set.
All elements are separated by tabs.
One gene set per line.
