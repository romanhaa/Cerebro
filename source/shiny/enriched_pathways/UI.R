tab_enriched_pathways <- tabItem(
    tabName = "enrichedPathways",
    box(
      title = tagList(
        p("Enriched pathways by sample",
          style = "padding-right: 5px; display: inline"
        ),
        actionButton(
          inputId = "enriched_pathways_by_sample_info", label = "info",
          icon = NULL, class = "btn-xs",
          title = "Show additional information for this panel."
        )
      ),
      status = "primary", solidHeader = TRUE, width = 12,
      collapsible = TRUE,
      uiOutput("enriched_pathways_by_sample_UI")
    ),
    box(
      title = tagList(
        p("Enriched pathways by cluster",
          style = "padding-right: 5px; display: inline"
        ),
        actionButton(
          inputId = "enriched_pathways_by_cluster_info", label = "info",
          icon = NULL, class = "btn-xs",
          title = "Show additional information for this panel."
        )
      ),
      status = "primary", solidHeader = TRUE, width = 12,
      collapsible = TRUE,
      uiOutput("enriched_pathways_by_cluster_UI")
    )
  )