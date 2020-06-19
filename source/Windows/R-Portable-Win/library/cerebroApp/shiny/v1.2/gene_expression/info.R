##----------------------------------------------------------------------------##
## Tab: Gene expression
##----------------------------------------------------------------------------##

expression_projection_info <- list(
  title = "Dimensional reduction",
  text = p(
    "Interactive projection of cells into 2-dimensional space based on their expression profile.",
    tags$ul(
      tags$li("Both tSNE and UMAP are frequently used algorithms for dimensional reduction in single cell transcriptomics. While they generally allow to make similar conclusions, some differences exist between the two (please refer to Google and/or literature, such as Becht E. et al., Dimensionality reduction for visualizing single-cell data using UMAP. Nature Biotechnology, 2018, 37, 38-44)."),
      tags$li("Cell color reflects the log-normalised expression of entered genes. If more than 1 gene is entered, the color reflects the average expression of all genes. Genes must be in separate lines or separated by a space, comma, or semicolon. Reported below the projection are the genes that are present and absent in this data set. Absent genes could either have been annotated with a different name or were not expressed in any of the cells. Matching of gene names is case-insensitive, that means Myc/MYC/myc are treated equally."),
      tags$li("Samples and clusters can be removed from the plot individually to highlight a contrast of interest."),
      tags$li("Cells can be plotted either randomly (which a more unbiased image) or in the order of expression (with highest expression plotted last), sometimes resulting in a more appealing figure."),
      tags$li("By default, the dot size is set to 15 without any transparency but both these attributes can be changed using the sliders on the left."),
      tags$li("The last 2 slider elements on the left can be used to resize the projection axes. This can be particularly useful when a projection contains a population of cell that is very far away from the rest and therefore creates a big empty space (which is not uncommon for UMAPs).")
    ),
    "The plot is interactive (drag and zoom) but depending on the computer of the user and the number of cells displayed it can become very slow."
  )
)

expression_details_selected_cells_info <- list(
  title = "Details of selected cells",
  text = p("Table containing (average) expression values of selected genes as well as selected meta data (sample, cluster, number of transcripts, number of expressed genes) for cells selected in the plot using the box or lasso selection tool. If you want the table to contain all cells in the data set, you must select all cells in the plot. The table can be saved to disk in CSV or Excel format for further analysis.")
)

expression_in_selected_cells_info <- list(
  title = "Expression levels in selected cells",
  text = p("This plot shows the log-normalised expression of selected genes for cells grouped by whether they were selected using the box or lasso selection tool. If more than 1 gene was provided, this reflects the average across all cells of each sample.")
)

expression_by_sample_info <- list(
  title = "Expression levels by sample",
  text = p("Log-normalised expression of genes inserted above by sample. If more than 1 gene was provided, this reflects the average across all cells of each sample.")
)

expression_by_cluster_info <- list(
  title = "Expression levels by cluster",
  text = p("Log-normalised expression of genes inserted above by cluster. If more than 1 gene was provided, this reflects the average across all cells of each cluster.")
)

expression_by_gene_info <- list(
  title = "Expression levels by gene",
  text = p("Log-normalised expression of 50 highest expressed genes inserted above. Shows mean across all cells.")
)
