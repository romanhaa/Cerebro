##----------------------------------------------------------------------------##
## Tab: Marker genes.
##----------------------------------------------------------------------------##

tab_marker_genes <- tabItem(
    tabName = "markerGenes",
    cerebroBox(
      title = tagList(
        boxTitle("Marker genes per sample"),
        cerebroInfoButton("marker_genes_by_sample_info")
      ),
      uiOutput("marker_genes_by_sample_UI")
    ),
    cerebroBox(
      title = tagList(
        boxTitle("Marker genes per cluster"),
        cerebroInfoButton("marker_genes_by_cluster_info")
      ),
      uiOutput("marker_genes_by_cluster_UI")
    )
  )