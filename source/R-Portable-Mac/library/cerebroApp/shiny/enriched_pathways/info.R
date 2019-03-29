##----------------------------------------------------------------------------##
## Panel: Enriched pathways
##----------------------------------------------------------------------------##

enriched_pathways_by_sample_info <- list(
  title = "Enriched pathways by sample",
  text = p("Using all marker genes identified for a respective sample, gene list enrichment analysis is performed using the Enrichr API, including gene ontology terms, KEGG and Wiki Pathways, BioCarta and many others. Terms are sorted based on the combined score. By default, the genes that overlap between the marker gene list and a term are not shown (for better visibility) but the column can be added using the 'Column visibility' button. For the details on the combined score is calculated, please refer to the Enrichr website and publication: http://amp.pharm.mssm.edu/Enrichr/.")
)

enriched_pathways_by_cluster_info <- list(
  title = "Enriched pathways by cluster",
  text = p("Using all marker genes identified for a respective cluster, gene list enrichment analysis is performed using the Enrichr API, including gene ontology terms, KEGG and Wiki Pathways, BioCarta and many others. Terms are sorted based on the combined score. By default, the genes that overlap between the marker gene list and a term are not shown (for better visibility) but the column can be added using the 'Column visibility' button. For the details on the combined score is calculated, please refer to the Enrichr website and publication: http://amp.pharm.mssm.edu/Enrichr/")
)
