##----------------------------------------------------------------------------##
## Panel: Gene id/symbol conversion.
##----------------------------------------------------------------------------##

tab_gene_id_conversion <- tabItem(
    tabName = "geneIdConversion",
    box(
      title = "Convert gene ID <-> gene symbol", status = "primary",
      solidHeader = TRUE, width = 12, collapsible = FALSE,
      tagList(
        selectInput(
          "geneIdConversion_organism", "Organism:",
          choices = c("mouse", "human")
        ),
        DT::dataTableOutput("gene_info")
      )
    )
  )