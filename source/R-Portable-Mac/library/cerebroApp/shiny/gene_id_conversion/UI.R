##----------------------------------------------------------------------------##
## Tab: Gene id/symbol conversion.
##----------------------------------------------------------------------------##

tab_gene_id_conversion <- tabItem(
  tabName = "geneIdConversion",
  cerebroBox(
    title = tagList(
      boxTitle("Convert gene ID <-> gene symbol"),
      cerebroInfoButton("geneIdConversion_info")
    ),
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