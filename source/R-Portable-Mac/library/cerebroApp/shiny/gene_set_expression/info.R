##----------------------------------------------------------------------------##
## Panel: Gene set expression
##----------------------------------------------------------------------------##

geneSetExpression_projection_info <- list(
  title = "Dimensional reduction",
  text = p(
    "Interactive projection of cells into 2-dimensional space based on their expression profile.",
    tags$ul(
      tags$li("Both tSNE and UMAP are frequently used algorithms for dimensional reduction in single cell transcriptomics. While they generally allow to make similar conclusions, some differences exist between the two (please refer to Google and/or literature, such as Becht E. et al., Dimensionality reduction for visualizing single-cell data using UMAP. Nature Biotechnology, 2018, 37, 38-44)."),
      tags$li("For human and murine data sets, all organism-specific gene sets from the MSigDB can be selected. If the experiment was performed in another organism, the murine gene sets will be available."),
      tags$li("Cell color reflects the average log-normalised expression of the genes in the selected gene set. Reported below the projection are the genes that are present and absent in this data set. Absent genes could either have been annotated with a different name or were not expressed in any of the cells. Matching of gene names is case-insensitive, that means Myc/MYC/myc are treated equally."),
      tags$li("Samples and clusters can be removed from the plot individually to highlight a contrast of interest."),
      tags$li("Cells can be plotted either randomly (which a more unbiased image) or in the order of expression (with highest expression plotted last), sometimes resulting in a more appealing figure."),
      tags$li("By default, the dot size is set to 15 without any transparency but both these attributes can be changed using the sliders on the left."),
      tags$li("The last 2 slider elements on the left can be used to resize the projection axes. This can be particularly useful when a projection contains a population of cell that is very far away from the rest and therefore creates a big empty space (which is not uncommon for UMAPs).")
    ),
    "The plot is interactive (drag and zoom) but depending on the computer of the user and the number of cells displayed it can become very slow."
  )
)

geneSetExpression_by_sample_info <- list(
  title = "Expression levels by sample",
  text = p("Average log-normalised expression of genes in selected gene set by sample.")
)

geneSetExpression_by_cluster_info <- list(
  title = "Expression levels by cluster",
  text = p("Average log-normalised expression of genes in selected gene set by cluster.")
)

