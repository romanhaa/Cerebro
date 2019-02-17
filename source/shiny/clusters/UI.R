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
      uiOutput("clusters_by_sample_UI")
    ),
    cerebroBox(
      title = tagList(
        boxTitle("Number of transcripts"),
        cerebroInfoButton("clusters_box_nUMI_info")
      ),
      uiOutput("clusters_box_nUMI_UI")
    ),
    cerebroBox(
      title = tagList(
        boxTitle("Number of expressed genes"),
        cerebroInfoButton("clusters_box_nGene_info")
      ),
      uiOutput("clusters_box_nGene_UI")
    ),
    cerebroBox(
      title = tagList(
        boxTitle("Mitochondrial gene expression"),
        cerebroInfoButton("clusters_box_percent_mt_info")
      ),
      uiOutput("clusters_box_percent_mt_UI")
    ),
    cerebroBox(
      title = tagList(
        boxTitle("Ribosomal gene expression"),
        cerebroInfoButton("clusters_box_percent_ribo_info")
      ),
      uiOutput("clusters_box_percent_ribo_UI")
    ),
    cerebroBox(
      title = tagList(
        boxTitle("Cell cycle analysis (Regev)"),
        cerebroInfoButton("clusters_by_cell_cycle_Regev_info")
      ),
      uiOutput("clusters_by_cell_cycle_Regev_UI")
    ),
    cerebroBox(
      title = tagList(
        boxTitle("Cell cycle analysis (Cyclone)"),
        cerebroInfoButton("clusters_by_cell_cycle_Cyclone_info")
      ),
      uiOutput("clusters_by_cell_cycle_Cyclone_UI")
    )
  )