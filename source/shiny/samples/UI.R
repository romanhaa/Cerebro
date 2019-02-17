##----------------------------------------------------------------------------##
## Tab: Samples.
##----------------------------------------------------------------------------##

tab_samples <- tabItem(
  tabName = "samples",
  cerebroBox(
    title = tagList(
      boxTitle("Samples by cluster"),
      cerebroInfoButton("samples_by_cluster_info")
    ),
    uiOutput("samples_by_cluster_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Number of transcripts"),
      cerebroInfoButton("samples_box_nUMI_info")
    ),
    uiOutput("samples_box_nUMI_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Number of expressed genes"),
      cerebroInfoButton("samples_box_nGene_info")
    ),
    uiOutput("samples_box_nGene_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Mitochondrial gene expression"),
      cerebroInfoButton("samples_box_percent_mt_info")
    ),
    uiOutput("samples_box_percent_mt_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Ribosomal gene expression"),
      cerebroInfoButton("samples_box_percent_ribo_info")
    ),
    uiOutput("samples_box_percent_ribo_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Cell cycle analysis (Regev)"),
      cerebroInfoButton("samples_by_cell_cycle_Regev_info")
    ),
    uiOutput("samples_by_cell_cycle_Regev_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Cell cycle analysis (Cyclone)"),
      cerebroInfoButton("samples_by_cell_cycle_Cyclone_info")
    ),
    uiOutput("samples_by_cell_cycle_Cyclone_UI")
  )
)