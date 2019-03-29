##----------------------------------------------------------------------------##
## Tab: Gene expression
##----------------------------------------------------------------------------##

tab_gene_expression <- tabItem(
    tabName = "geneExpression",
    tags$script('
      $(document).on("keyup", function(event) {
        if ( event.keyCode == 13 ) {
          event.preventDefault();
          Shiny.onInputChange("keyPressed", Math.random());
        }
        if ( event.keyCode == 32 ) {
          Shiny.onInputChange("keyPressed", Math.random());
        }
      });
    '),
    tagList(
      fluidRow(
        column(width = 3, offset = 0, style = "padding: 0px;",
          cerebroBox(
            title = "Input parameters",
            tagList(
              uiOutput("expression_UI"),
              uiOutput("expression_scales")
            )
          )
        ),
        column(width = 9, offset = 0, style = "padding: 0px;",
          cerebroBox(
            title = tagList(
              boxTitle("Dimensional reduction"),
              actionButton(
                inputId = "expression_projection_info",
                label = "info",
                icon = NULL,
                class = "btn-xs",
                title = "Show additional information for this panel.",
                style = "margin-right: 5px"
              ),
              actionButton(
                inputId = "expression_projection_export",
                label = "export to PDF",
                icon = NULL,
                class = "btn-xs",
                title = "Export dimensional reduction to PDF file."
              )
            ),
            tagList(
              plotly::plotlyOutput(
                "expression_projection_plotly",
                width = "auto",
                height = "85vh"
              ),
              tags$br(),
              htmlOutput("expression_genes_displayed")
            )
          )
        )
      ),
      fluidRow(
        cerebroBox(
          title = tagList(
            boxTitle("Expression levels by sample"),
            cerebroInfoButton("expression_by_sample_info")
          ),
          plotly::plotlyOutput("expression_by_sample")
        )
      ),
      fluidRow(
        cerebroBox(
          title = tagList(
            boxTitle("Expression levels by cluster"),
            cerebroInfoButton("expression_by_cluster_info")
          ),
          plotly::plotlyOutput("expression_by_cluster")
        )
      ),
      fluidRow(
        cerebroBox(
          title = tagList(
            boxTitle("Expression levels by gene"),
            cerebroInfoButton("expression_by_gene_info")
          ),
          plotly::plotlyOutput("expression_by_gene")
        )
      )
    )
  )