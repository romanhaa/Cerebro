##----------------------------------------------------------------------------##
## Panel: Most expressed genes.
##----------------------------------------------------------------------------##

most_expressed_genes_by_sample_info <- list(
  title = "Most expressed genes per sample",
  text = p("Table of top 100 most expressed genes in each sample. For example, if gene XY contributes with 5% to the total expression, that means 5% of all transcripts found in all cells of this sample come from that respective gene. These lists can help to identify/verify the dominant cell types.")
)

most_expressed_genes_by_cluster_info <- list(
  title = "Most expressed genes per cluster",
  text = p("Table of top 100 most expressed genes in each cluster. For example, if gene XY contributes with 5% to the total expression, that means 5% of all transcripts found in all cells of this cluster come from that respective gene. These lists can help to identify/verify the dominant cell types.")
)

