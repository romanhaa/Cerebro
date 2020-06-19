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

overview_details_selected_cells_table_info <- list(
  title = "Details of selected cells",
  text = p("Table containing meta data (some columns may be hidden, check the 'Column visibility' button) for cells selected in the plot using the box or lasso selection tool. If you want the table to contain all cells in the data set, you must select all cells in the plot. The table can be saved to disk in CSV or Excel format for further analysis.")
)

overview_details_selected_cells_plot_info <- list(
  title = "Plot of selected cells",
  text = p("Depending on the variable selected to color cells in the dimensional reduction, this plot will show different things. If you select a categorical variable, e.g. 'sample' or 'cluster', you will get a bar plot showing which groups the cells selected with the box or lasso tool come from. Instead, if you select a continuous variable, e.g. the number of transcripts (nUMI), you will see a violin/box plot showing the distribution of that variable in the selected vs. non-selected cells.")
)
