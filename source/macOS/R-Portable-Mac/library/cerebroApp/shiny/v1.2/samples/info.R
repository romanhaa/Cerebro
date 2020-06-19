##----------------------------------------------------------------------------##
## Tab: Samples.
##----------------------------------------------------------------------------##

samples_overview_info <- list(
  title = "Overview of samples",
  text = p("Overview of samples with preliminary information.")
)

samples_by_cluster_info <- list(
  title = "Samples by cluster",
  text = p("Percentage bar plot representation of the table shown above. Allows to see which clusters contribute most strongly to each sample. Clusters can be removed from the plot by clicking on them in the legend.")
)

samples_nUMI_info <- list(
  title = "Number of transcripts",
  text = p("Violin plot of the number of transcripts (UMIs) found in each sample.")
)

samples_nGene_info <- list(
  title = "Number of expressed genes",
  text = p("Violin plot of the number of expressed genes found in each sample.")
)

samples_percent_mt_info <- list(
  title = "Mitochondrial gene expression",
  text = p("Violin plot of the percentage of mitochondrial gene expression found in each sample. This reflects the contribution of mitochondrial transcripts to the entire transcriptome in each cell. A list of all genes considered to be mitochondrial can be found in the 'Sample info' tab on the left.")
)

samples_percent_ribo_info <- list(
  title = "Ribosomal gene expression",
  text = p("Violin plot of the percentage of ribosomal gene expression found in each sample. This reflects the contribution of ribosomal transcripts to the entire transcriptome in each cell. A list of all genes considered to be ribosomal can be found in the 'Sample info' tab on the left.")
)

samples_by_cell_cycle_seurat_info <- list(
  title = "Cell cycle analysis (Seurat)",
  text = p("Cell cycle distribution by sample using the method embedded in the Seurat framework. For each cell, it calculates scores for both G2M and S phase based on lists of genes (see 'Analysis info' tab on the left) and assigns the cell cycle phase on the basis of these scores.")
)

samples_by_cell_cycle_cyclone_info <- list(
  title = "Cell cycle analysis (Cyclone)",
  text = p("Cell cycle distribution by sample using the machine learning-based Cyclone method published by Scialdone et al (2015). It assigns the cell cycle phase based on scores calculated using relative expression of lists of gene pairs. In contrast to the Seurat method, scores are calculated for G1 and G2M phase. Cells with a low score for both are assigned S phase. Inability to predict the cell cycle phase for a given cell with this method is most likely a result of very few expressed genes in the respective cell.")
)
