##----------------------------------------------------------------------------##
## Tab: Gene id/symbol conversion.
##----------------------------------------------------------------------------##

output[["gene_info"]] <- DT::renderDataTable({
  if ( input[["geneIdConversion_organism"]] == "mouse" ) {
    conversion_table <- read.table("resources/mm10_gene_ID_name.txt",
      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else if ( input[["geneIdConversion_organism"]] == "human" ) {
    conversion_table <- read.table("resources/hg38_gene_ID_name.txt",
      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  }
  DT::datatable(
    conversion_table,
    filter = "none",
    selection = "multiple",
    escape = FALSE,
    autoHideNavigation = TRUE,
    rownames = FALSE,
    options = list(
      scrollX = FALSE,
      dom = "Bfrtip",
      lengthMenu = c(15, 30, 50, 100),
      pageLength = 50
    )
  )
})

# info box
observeEvent(input[["geneIdConversion_info"]], {
  showModal(
    modalDialog(
      geneIdConversion_info[["text"]],
      title = geneIdConversion_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})