##----------------------------------------------------------------------------##
## Panel: Overview.
##----------------------------------------------------------------------------##

overview_projection_info <- list(
  title = "Dimensional reduction",
  text = p(
    "Interactive projection of cells into 2-dimensional space based on their expression profile.",
    tags$ul(
      tags$li("Both tSNE and UMAP are frequently used algorithms for dimensional reduction in single cell transcriptomics. While they generally allow to make similar conclusions, some differences exist between the two (please refer to Google and/or literature, such as Becht E. et al., Dimensionality reduction for visualizing single-cell data using UMAP. Nature Biotechnology, 2018, 37, 38-44)."),
      tags$li("Cells can be colored by the sample they came from, the cluster they were assigned, the number of transcripts or expressed genes, percentage of mitochondrial and ribosomal gene expression, an apoptotic score (calculated based on the expression of few marker genes; more info in the 'Sample info' tab on the left), or cell cycle status (determined using the Seurat and Cyclone method)."),
      tags$li("Confidence ellipses show the 95% confidence regions."),
      tags$li("Samples and clusters can be removed from the plot individually to highlight a contrast of interest."),
      tags$li("By default, the dot size is set to 15 without any transparency but both these attributes can be changed using the sliders on the left. The dot size can also be set to reflect the number of transcripts or expressed genes."),
      tags$li("The last 2 slider elements on the left can be used to resize the projection axes. This can be particularly useful when a projection contains a population of cell that is very far away from the rest and therefore creates a big empty space (which is not uncommon for UMAPs).")
    ),
    "The plot is interactive (drag and zoom) but depending on the computer of the user and the number of cells displayed it can become very slow."
  )
)




