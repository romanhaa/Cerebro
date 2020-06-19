##----------------------------------------------------------------------------##
## Panel: Enriched pathways
##----------------------------------------------------------------------------##

enriched_pathways_by_sample_info <- list(
  title = "Enriched pathways by sample",
  text = HTML(
    "<p>At the moment, Cerebro supports to perform pathway enrichment analysis through two methods which work in different ways: Enrichr, GSVA.<br>
    <br>
    <b>Enrichr</b><br>
    Using all marker genes identified for a respective sample, gene list enrichment analysis is performed using the Enrichr API, including gene ontology terms, KEGG and Wiki Pathways, BioCarta and many others. Terms are sorted based on the combined score. By default, the genes that overlap between the marker gene list and a term are not shown (for better visibility) but the column can be added using the 'Column visibility' button. For the details on the combined score is calculated, please refer to the <a target='_blank' href='http://amp.pharm.mssm.edu/Enrichr/'>Enrichr website</a> and publication.<br>
    <br>
    <b>GSVA</b><br>
    GSVA (Gene Set Variation Analysis) is a method to perform gene set enrichment analysis. Diaz-Mejia and colleagues found GSVA to perform well compared to other tools ('Evaluation of methods to assign cell type labels to cell clusters from single-cell RNA-sequencing data' F1000Research, 2019). Statistics (p-value and adj. p-value) are calculated as done by Diaz-Mejia et al. Columns with the gene lists for each term and the enrichment score are hidden by default but can be made visible through the 'Column visibility' button. More details about GSVA can be found on the <a target='_blank', href='https://bioconductor.org/packages/release/bioc/html/GSVA.html'>GSVA Bioconductor page</a>."
  )
)

enriched_pathways_by_cluster_info <- list(
  title = "Enriched pathways by cluster",
  text = HTML(
    "<p>At the moment, Cerebro supports to perform pathway enrichment analysis through two methods which work in different ways: Enrichr, GSVA.<br>
    <br>
    <b>Enrichr</b><br>
    Using all marker genes identified for a respective cluster, gene list enrichment analysis is performed using the Enrichr API, including gene ontology terms, KEGG and Wiki Pathways, BioCarta and many others. Terms are sorted based on the combined score. By default, the genes that overlap between the marker gene list and a term are not shown (for better visibility) but the column can be added using the 'Column visibility' button. For the details on the combined score is calculated, please refer to the <a target='_blank' href='http://amp.pharm.mssm.edu/Enrichr/'>Enrichr website</a> and publication/<br>
    <br>
    <b>GSVA</b><br>
    GSVA (Gene Set Variation Analysis) is a method to perform gene set enrichment analysis. Diaz-Mejia and colleagues found GSVA to perform well compared to other tools ('Evaluation of methods to assign cell type labels to cell clusters from single-cell RNA-sequencing data' F1000Research, 2019). Statistics (p-value and adj. p-value) are calculated as done by Diaz-Mejia et al. Columns with the gene lists for each term and the enrichment score are hidden by default but can be made visible through the 'Column visibility' button. More details about GSVA can be found on the <a target='_blank', href='https://bioconductor.org/packages/release/bioc/html/GSVA.html'>GSVA Bioconductor page</a>."
  )
)
