##----------------------------------------------------------------------------##
## Tab: Analysis info.
##----------------------------------------------------------------------------##

# general info
output[["sample_info_general"]] <- renderText({
  info <- paste0(
    "<strong><u>General</u></strong>",
    "<ul>",
      "<li><b>Date of analysis:</b> ",
      sample_data()$experiment$date_of_analysis,
      "<li><b>Experiment name:</b> ",
      sample_data()$experiment$experiment_name,
      "<li><b>Organism:</b> ",
      sample_data()$experiment$organism,
      "<li><b>Samples:</b> ",
      paste0(sample_data()$samples$overview$sample, collapse = ", "),
      "<li><b>Number of clusters found:</b> ",
      nrow(sample_data()$clusters$overview),
    "</ul>",
    "<strong><u>Parameters</u></strong>",
    "<ul>",
      "<li><b>Discard genes in fewer than X cells:</b> ",
      sample_data()$parameters$discard_genes_expressed_in_fewer_cells_than,
      "<li><b>Keep mitochondrial genes:</b> ",
      sample_data()$parameters$keep_mitochondrial_genes,
      "<li><b>Min/max # of UMI:</b> ",
      paste0(
        sample_data()$parameters$filtering$UMI_min, " / ",
        sample_data()$parameters$filtering$UMI_max
      ),
      "<li><b>Min/max # of expressed genes:</b> ",
      paste0(
        sample_data()$parameters$filtering$genes_min, " / ",
        sample_data()$parameters$filtering$genes_max
      ),
      "<li><b>Cluster resolution: </b>",
      sample_data()$parameters$cluster_resolution,
      "<li><b>Number of principal components: </b>",
      sample_data()$parameters$number_PCs,
      "<li><b>Variables to regress: </b>",
      sample_data()$parameters$variables_to_regress_out,
      "<li><b>tSNE perplexity: </b>",
      sample_data()$parameters$tSNE_perplexity,
    "</ul>",
    "<strong><u>Gene lists</u></strong>",
    "<ul>",
      "<li><b>Mitochondrial genes:</b> ",
      paste0(sample_data()$gene_lists$mitochondrial_genes, collapse = ", "),
      "<li><b>Ribosomal genes:</b> ",
      paste0(sample_data()$gene_lists$ribosomal_genes, collapse = ", "),
      "<li><b>S phase genes:</b> ",
      paste0(sample_data()$gene_lists$S_phase_genes, collapse = ", "),
      "<li><b>G2M phase genes:</b> ",
      paste0(sample_data()$gene_lists$G2M_phase_genes, collapse = ", "),
    "</ul>",
    "<strong><u>Marker genes</u></strong>",
    "<ul>",
      "<li><b>Only positive:</b> ",
      sample_data()$marker_genes$parameters$only_positive,
      "<li><b>Fraction of cells in group of interest that must express marker gene:</b> ",
      sample_data()$marker_genes$parameters$minimum_percentage,
      "<li><b>LogFC threshold:</b> ",
      sample_data()$marker_genes$parameters$logFC_threshold,
      "<li><b>p-value threshold:</b> ",
      sample_data()$marker_genes$parameters$p_value_threshold,
      "<li><b>Test:</b> ",
      sample_data()$marker_genes$parameters$test,
    "</ul>",
    "<strong><u>Pathway enrichment</u></strong>",
    "<ul>",
      "<li><b>Enrichr:</b>",
      "<ul>",
        "<li><b>Databases:</b> ",
        paste0(sample_data()$enriched_pathways$enrichr$parameters$databases, collapse = ", "),
        "<li><b>Adj. p-value cut-off:</b> ",
        sample_data()$enriched_pathways$enrichr$parameters$adj_p_cutoff,
        "<li><b>Max. terms:</b> ",
        sample_data()$enriched_pathways$enrichr$parameters$max_terms,
      "</ul>",
      "<li><b>GSVA:</b>",
      "<ul>",
        "<li><b>GMT file:</b> ",
        sample_data()$enriched_pathways$GSVA$parameters$GMT_file,
        "<li><b>p-value cut-off:</b> ",
        sample_data()$enriched_pathways$GSVA$parameters$thresh_p_val,
        "<li><b>q-value cut-off:</b> ",
        sample_data()$enriched_pathways$GSVA$parameters$thresh_q_val,
      "</ul>",
    "</ul>"
  )
  info_R_raw <- sample_data()$technical_info$R
  info_R <- c()
  for ( i in 1:length(info_R_raw) ) {
    info_R <- paste(info_R, "<br>", info_R_raw[i])
  }
  paste0(
    info,
    "<strong><u>Technical info (package versions)</u></strong>",
    "<ul>",
      "<li><strong>Seurat version:</strong> ",
      sample_data()$technical_info$seurat_version,
      "<li><strong>Session info:</strong> ",
    "</ul>",
    "<pre>",
    info_R,
    "</pre>"
  )
})

# R info
output[["sample_info_R"]] <- renderPrint({
  if ( !is.null(sample_data()$technical_info$R) ) {
    capture.output(sample_data()$technical_info$R)
  } else {
    print("Not available")
  }
})
