##----------------------------------------------------------------------------##
## Tab: Enriched pathways
##----------------------------------------------------------------------------##

tab_enriched_pathways <- tabItem(
  tabName = "enrichedPathways",
  cerebroBox(
    title = tagList(
      boxTitle("Enriched pathways by sample"),
      cerebroInfoButton("enriched_pathways_by_sample_info")
    ),
    tagList(
      fluidRow(
        column(4,
          uiOutput("enriched_pathways_by_sample_select_sample_UI")
        ),
        column(8,
          uiOutput("enriched_pathways_by_sample_select_db_UI")
        )
      ),
      uiOutput("enriched_pathways_by_sample_UI")
    )
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Enriched pathways by cluster"),
      cerebroInfoButton("enriched_pathways_by_cluster_info")
    ),
    tagList(
      fluidRow(
        column(4,
          uiOutput("enriched_pathways_by_cluster_select_cluster_UI")
        ),
        column(8,
          uiOutput("enriched_pathways_by_cluster_select_db_UI")
        )
      ),
      uiOutput("enriched_pathways_by_cluster_UI")
    )
  )
)