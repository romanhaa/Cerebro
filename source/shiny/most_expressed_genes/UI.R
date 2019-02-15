##----------------------------------------------------------------------------##
## Panel: Most expressed genes.
##----------------------------------------------------------------------------##

tab_most_expressed_genes <- tabItem(
    tabName = "mostExpressedGenes",
    box(
      title = tagList(
        p(
          "Most expressed genes per sample",
          style = "padding-right: 5px; display: inline"
        ),
        actionButton(
          inputId = "most_expressed_genes_by_sample_info", label = "info",
          icon = NULL, class = "btn-xs",
          title = "Show additional information for this panel."
        )
      ),
      status = "primary", solidHeader = TRUE, width = 12,
      collapsible = TRUE,
      uiOutput("most_expressed_genes_by_sample_UI")
    ),
    box(
      title = tagList(
        p(
          "Most expressed genes per cluster",
          style = "padding-right: 5px; display: inline"
        ),
        actionButton(
          inputId = "most_expressed_genes_by_cluster_info", label = "info",
          icon = NULL, class = "btn-xs",
          title = "Show additional information for this panel."
        )
      ),
      status = "primary", solidHeader = TRUE, width = 12,
      collapsible = TRUE,
      uiOutput("most_expressed_genes_by_cluster_UI")
    )
  )