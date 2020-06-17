##--------------------------------------------------------------------------##
## Tab: Clusters.
##--------------------------------------------------------------------------##

tab_clusters <- tabItem(
  tabName = "clusters",
  cerebroBox(
    title = tagList(
      boxTitle("Cluster tree"),
      cerebroInfoButton("clusters_tree_info")
    ),
    uiOutput("clusters_tree_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Clusters by samples"),
      cerebroInfoButton("clusters_by_sample_info")
    ),
    tagList(
      uiOutput("clusters_by_sample_UI_buttons"),
      uiOutput("clusters_by_sample_UI_rest")
    )
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Number of transcripts"),
      cerebroInfoButton("clusters_nUMI_info")
    ),
    uiOutput("clusters_nUMI_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Number of expressed genes"),
      cerebroInfoButton("clusters_nGene_info")
    ),
    uiOutput("clusters_nGene_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Mitochondrial gene expression"),
      cerebroInfoButton("clusters_percent_mt_info")
    ),
    uiOutput("clusters_percent_mt_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Ribosomal gene expression"),
      cerebroInfoButton("clusters_percent_ribo_info")
    ),
    uiOutput("clusters_percent_ribo_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Cell cycle analysis by Seurat"),
      cerebroInfoButton("clusters_by_cell_cycle_seurat_info")
    ),
    tagList(
      uiOutput("clusters_by_cell_cycle_seurat_UI_buttons"),
      uiOutput("clusters_by_cell_cycle_seurat_UI_rest")
    )
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Cell cycle analysis by Cyclone"),
      cerebroInfoButton("clusters_by_cell_cycle_cyclone_info")
    ),
    tagList(
      uiOutput("clusters_by_cell_cycle_cyclone_UI_buttons"),
      uiOutput("clusters_by_cell_cycle_cyclone_UI_rest")
    )
  )
)
