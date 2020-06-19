##----------------------------------------------------------------------------##
## Tab: Clusters.
##----------------------------------------------------------------------------##

clusters_tree_info <- list(
  title = "Cluster tree",
  text = p("Cluster tree reflecting the similarity of clusters based on their expression profiles. Instead of using the expression values, the correlation is calculated using the user-specified number of principal components (see 'Sample info' tab on the left).")
)

clusters_by_sample_info <- list(
  title = "Clusters by samples",
  text = p("Percentage bar plot representation of the table shown above. Allows to see which samples contribute most strongly to each cluster. Samples can be removed from the plot by clicking on them in the legend.")
)

clusters_nUMI_info <- list(
  title = "Number of transcripts",
  text = p("Violin plot of the number of transcripts (UMIs) found in each cluster.")
)

clusters_nGene_info <- list(
  title = "Number of expressed genes",
  text = p("Violin plot of the number of expressed genes found in each cluster.")
)

clusters_percent_mt_info <- list(
  title = "Mitochondrial gene expression",
  text = p("Violin plot of the percentage of mitochondrial gene expression found in each cluster. This reflects the contribution of mitochondrial transcripts to the entire transcriptome in each cell. A list of all genes considered to be mitochondrial can be found in the 'Sample info' tab on the left.")
)

clusters_percent_ribo_info <- list(
  title = "Ribosomal gene expression",
  text = p("Violin plot of the percentage of ribosomal gene expression found in each cluster. This reflects the contribution of ribosomal transcripts to the entire transcriptome in each cell. A list of all genes considered to be ribosomal can be found in the 'Sample info' tab on the left.")
)

clusters_by_cell_cycle_seurat_info <- list(
  title = "Cell cycle analysis (Seurat)",
  text = p("Cell cycle distribution by cluster using the method embedded in the Seurat framework. For each cell, it calculates scores for both G2M and S phase based on lists of genes (see 'Analysis info' tab on the left) and assigns the cell cycle phase on the basis of these scores.")
)

clusters_by_cell_cycle_cyclone_info <- list(
  title = "Cell cycle analysis (Cyclone)",
  text = p("Cell cycle distribution by cluster using the machine learning-based Cyclone method published by Scialdone et al (2015). It assigns the cell cycle phase based on scores calculated using relative expression of lists of gene pairs. In contrast to the Seurat method, scores are calculated for G1 and G2M phase. Cells with a low score for both are assigned S phase.")
)
