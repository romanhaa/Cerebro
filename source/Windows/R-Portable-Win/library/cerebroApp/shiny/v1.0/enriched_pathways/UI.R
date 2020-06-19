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
    uiOutput("enriched_pathways_by_sample_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Enriched pathways by cluster"),
      cerebroInfoButton("enriched_pathways_by_cluster_info")
    ),
    uiOutput("enriched_pathways_by_cluster_UI")
  )
)
