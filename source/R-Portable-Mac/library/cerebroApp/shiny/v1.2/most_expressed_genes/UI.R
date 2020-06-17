##----------------------------------------------------------------------------##
## Tab: Most expressed genes.
##----------------------------------------------------------------------------##

tab_most_expressed_genes <- tabItem(
  tabName = "mostExpressedGenes",
  cerebroBox(
    title = tagList(
      boxTitle("Most expressed genes per sample"),
      cerebroInfoButton("most_expressed_genes_by_sample_info")
    ),
    uiOutput("most_expressed_genes_by_sample_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Most expressed genes per cluster"),
      cerebroInfoButton("most_expressed_genes_by_cluster_info")
    ),
    uiOutput("most_expressed_genes_by_cluster_UI")
  )
)
