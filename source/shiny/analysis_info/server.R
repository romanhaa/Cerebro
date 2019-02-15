##----------------------------------------------------------------------------##
## Panel: Analysis info.
##----------------------------------------------------------------------------##

# general info
output$sample_info_general <- renderText({
  info <- paste0(
      "<ul>",
        "<li><b>Experiment name:</b> ",
        sample_data()$experiment$experiment_name,
        "<li><b>Organism:</b> ",
        sample_data()$experiment$organism,
        "<li><b>Samples:</b> ",
        paste0(sample_data()$samples$overview$sample, collapse = ", "),
        "<li><b># of clusters found:</b> ",
        nrow(sample_data()$clusters$overview),
        "<li><b>Min. cells:</b> ",
        sample_data()$parameters$discard_genes_expressed_in_fewer_cells_than,
        "<li><b>Keep mitochondrial genes:</b> ",
        sample_data()$parameters$keep_mitochondrial_genes,
        "<li><b>Variables to regress:</b> ",
        sample_data()$parameters$variables_to_regress_out,
        "<li><b>Number of principal components:</b> ",
        sample_data()$parameters$number_PCs,
        "<li><b>tSNE perplexity:</b> ",
        sample_data()$parameters$tSNE_perplexity,
        "<li><b>Cluster resolution:</b> ",
        sample_data()$parameters$cluster_resolution,
        "<li><b>enrichR databases:</b> ",
        paste0(sample_data()$parameters$enrichr_dbs, collapse = ", "),
        "<li><b>Mitochondrial genes:</b> ",
        paste0(sample_data()$gene_lists$mitochondrial_genes, collapse = ", "),
        "<li><b>Ribosomal genes:</b> ",
        paste0(sample_data()$gene_lists$ribosomal_genes, collapse = ", "),
        "<li><b>Genes used for apoptotic score:</b> ",
        paste0(sample_data()$gene_lists$apoptosis_genes, collapse = ", "),
        "<li><b>S phase genes:</b> ",
        paste0(sample_data()$gene_lists$S_phase_genes, collapse = ", "),
        "<li><b>G2M phase genes:</b> ",
        paste0(sample_data()$gene_lists$G2M_phase_genes, collapse = ", "),
        "<li><b>Min/max # of UMI:</b> ",
        paste0(
          sample_data()$parameters$filtering$UMI_min, " / ",
          sample_data()$parameters$filtering$UMI_max),
        "<li><b>Min/max # of expressed genes:</b> ",
        paste0(
          sample_data()$parameters$filtering$genes_min, " / ",
          sample_data()$parameters$filtering$genes_max),
        "<li><b>Min/max % of mitochonrial genes:</b> ",
        paste0(
          sample_data()$parameters$filtering$percent_mt_min, " / ",
          sample_data()$parameters$filtering$percent_mt_min),
        "<li><b>Min/max % of ribosomal genes:</b> ",
        paste0(
          sample_data()$parameters$filtering$percent_ribo_min, " / ",
          sample_data()$parameters$filtering$percent_ribo_min),
      "</ul>"
    )
  info_R_raw <- sample_data()$technical_info$R
  info_R <- c()
  for ( i in 1:length(info_R_raw) ) {
    info_R <- paste(info_R, "<br>", info_R_raw[i])
  }
  paste0(
    info,
    "<br><b>R environment and packages use in analysis:</b><br><pre>",
    info_R,
    "</pre>"
  )
})

# R info
output$sample_info_R <- renderPrint({
  if ( !is.null(sample_data()$technical_info$R) ) {
    capture.output(sample_data()$technical_info$R)
  } else {
    print("Not available")
  }
})
