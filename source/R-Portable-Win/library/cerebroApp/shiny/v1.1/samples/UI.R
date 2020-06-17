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
      cerebroInfoButton("samples_nUMI_info")
    ),
    uiOutput("samples_nUMI_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Number of expressed genes"),
      cerebroInfoButton("samples_nGene_info")
    ),
    uiOutput("samples_nGene_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Mitochondrial gene expression"),
      cerebroInfoButton("samples_percent_mt_info")
    ),
    uiOutput("samples_percent_mt_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Ribosomal gene expression"),
      cerebroInfoButton("samples_percent_ribo_info")
    ),
    uiOutput("samples_percent_ribo_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Cell cycle analysis by Seurat"),
      cerebroInfoButton("samples_by_cell_cycle_seurat_info")
    ),
    uiOutput("samples_by_cell_cycle_seurat_UI")
  ),
  cerebroBox(
    title = tagList(
      boxTitle("Cell cycle analysis by Cyclone"),
      cerebroInfoButton("samples_by_cell_cycle_cyclone_info")
    ),
    uiOutput("samples_by_cell_cycle_cyclone_UI")
  )
)
