##----------------------------------------------------------------------------##
## Tab: Gene id/symbol conversion.
##----------------------------------------------------------------------------##

tab_gene_id_conversion <- tabItem(
  tabName = "geneIdConversion",
  cerebroBox(
    title = "Convert gene ID <-> gene symbol",
    tagList(
      selectInput(
        "geneIdConversion_organism",
        "Organism:",
        choices = c("mouse", "human")
      ),
      DT::dataTableOutput("gene_info")
    )
  )
)