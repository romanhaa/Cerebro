##--------------------------------------------------------------------------##
## Tab: Clusters.
##--------------------------------------------------------------------------##

tab_clusters <- tabItem(
    tabName = "clusters",
    box(
      title = tagList(
        p("Cluster tree", style = "padding-right: 5px; display: inline"),
        actionButton(
          inputId = "clusters_tree_info", label = "info", icon = NULL,
          class = "btn-xs",
          title = "Show additional information for this panel."
        )
      ),
      status = "primary", solidHeader = TRUE, width = 12,
      collapsible = TRUE,
      uiOutput("clusters_tree_UI")
    ),
    box(
      title = tagList(
        p(
          "Clusters by samples",
          style = "padding-right: 5px; display: inline"
        ),
        actionButton(
          inputId = "clusters_by_sample_info", label = "info", icon = NULL,
          class = "btn-xs",
          title = "Show additional information for this panel."
        )
      ),
      status = "primary", solidHeader = TRUE, width = 12,
      collapsible = TRUE,
      uiOutput("clusters_by_sample_UI")
    ),
    box(
      title = tagList(
        p(
          "Number of transcripts",
          style = "padding-right: 5px; display: inline"
        ),
        actionButton(
          inputId = "clusters_box_nUMI_info", label = "info", icon = NULL,
          class = "btn-xs",
          title = "Show additional information for this panel."
        )
      ),
      status = "primary", solidHeader = TRUE, width = 12,
      collapsible = TRUE,
      uiOutput("clusters_box_nUMI_UI")
    ),
    box(
      title = tagList(
        p(
          "Number of expressed genes",
          style = "padding-right: 5px; display: inline"
        ),
        actionButton(
          inputId = "clusters_box_nGene_info", label = "info", icon = NULL,
          class = "btn-xs",
          title = "Show additional information for this panel."
        )
      ),
      status = "primary", solidHeader = TRUE, width = 12,
      collapsible = TRUE,
      uiOutput("clusters_box_nGene_UI")
    ),
    box(
      title = tagList(
        p(
          "Mitochondrial gene expression",
          style = "padding-right: 5px; display: inline"
        ),
        actionButton(
          inputId = "clusters_box_percent_mt_info", label = "info",
          icon = NULL, class = "btn-xs",
          title = "Show additional information for this panel."
        )
      ),
      status = "primary", solidHeader = TRUE, width = 12,
      collapsible = TRUE,
      uiOutput("clusters_box_percent_mt_UI")
    ),
    box(
      title = tagList(
        p(
          "Ribosomal gene expression",
          style = "padding-right: 5px; display: inline"
        ),
        actionButton(
          inputId = "clusters_box_percent_ribo_info", label = "info",
          icon = NULL, class = "btn-xs",
          title = "Show additional information for this panel."
        )
      ),
      status = "primary", solidHeader = TRUE, width = 12,
      collapsible = TRUE,
      uiOutput("clusters_box_percent_ribo_UI")
    ),
    box(
      title = tagList(
        p(
          "Cell cycle analysis (Regev)",
          style = "padding-right: 5px; display: inline"
        ),
        actionButton(
          inputId = "clusters_by_cell_cycle_Regev_info", label = "info",
          icon = NULL, class = "btn-xs",
          title = "Show additional information for this panel."
        )
      ),
      status = "primary", solidHeader = TRUE, width = 12,
      collapsible = TRUE,
      uiOutput("clusters_by_cell_cycle_Regev_UI")
    ),
    box(
      title = tagList(
        p(
          "Cell cycle analysis (Cyclone)",
          style = "padding-right: 5px; display: inline"
        ),
        actionButton(
          inputId = "clusters_by_cell_cycle_Cyclone_info", label = "info",
          icon = NULL, class = "btn-xs",
          title = "Show additional information for this panel."
        )
      ),
      status = "primary", solidHeader = TRUE, width = 12,
      collapsible = TRUE,
      uiOutput("clusters_by_cell_cycle_Cyclone_UI")
    )
  )