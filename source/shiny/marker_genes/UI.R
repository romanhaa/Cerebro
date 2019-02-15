##----------------------------------------------------------------------------##
## Panel: Marker genes.
##----------------------------------------------------------------------------##

tab_marker_genes <- tabItem(
    tabName = "markerGenes",
    box(
      title = tagList(
        p("Marker genes per sample",
          style = "padding-right: 5px; display: inline"
        ),
        actionButton(
          inputId = "marker_genes_by_sample_info", label = "info",
          icon = NULL, class = "btn-xs",
          title = "Show additional information for this panel."
        )
      ),
      status = "primary", solidHeader = TRUE, width = 12,
      collapsible = TRUE,
      uiOutput("marker_genes_by_sample_UI")
    ),
    box(
      title = tagList(
        p("Marker genes per cluster",
          style = "padding-right: 5px; display: inline"
        ),
        actionButton(
          inputId = "marker_genes_by_cluster_info", label = "info",
          icon = NULL, class = "btn-xs",
          title = "Show additional information for this panel."
        )
      ),
      status = "primary", solidHeader = TRUE, width = 12,
      collapsible = TRUE,
      uiOutput("marker_genes_by_cluster_UI")
    )
  )