##----------------------------------------------------------------------------##
## Panel: Gene expression
##----------------------------------------------------------------------------##

tab_gene_expression <- tabItem(
    tabName = "geneExpression",
    tags$script('
      $(document).on("keyup", function(e) {
        if (e.keyCode == 13) {
          Shiny.onInputChange("keyPressed", Math.random());
        }
      });
    '),
    tagList(
      fluidRow(
        column(width = 3, offset = 0, style = "padding: 0px;",
          box(
            title = "Input parameters", status = "primary",
            solidHeader = TRUE, width = 12, collapsible = TRUE,
            tagList(
              uiOutput("expression_UI"),
              uiOutput("expression_scales")
            )
          )
        ),
        column(width = 9, offset = 0, style = "padding: 0px;",
          box(
            title = tagList(
              p("Dimensional reduction",
                style = "padding-right: 5px; display: inline"
              ),
              actionButton(
                inputId = "expression_projection_info", label = "info",
                icon = NULL, class = "btn-xs",
                title = "Show additional information for this panel.",
                style = "margin-right: 5px"
              ),
              actionButton(
                inputId = "expression_projection_export",
                label = "export to PDF", icon = NULL, class = "btn-xs",
                title = "Export dimensional reduction to PDF file."
              )
            ),
            status = "primary", solidHeader = TRUE, width = 12,
            collapsible = TRUE,
            tagList(
              # scatterD3::scatterD3Output(
              #   "expression_projection", height = "720px"
              # ),
              plotly::plotlyOutput(
                "expression_projection_plotly", height = "720px"
              ),
              tags$br(),
              htmlOutput("expression_genes_displayed")
            )
          )
        )
      ),
      fluidRow(
        box(
          title = tagList(
            p("Expression levels by sample",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "expression_by_sample_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          plotly::plotlyOutput("expression_by_sample")
        )
      ),
      fluidRow(
        box(
          title = tagList(
            p("Expression levels by cluster",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "expression_by_cluster_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          plotly::plotlyOutput("expression_by_cluster")
        )
      )
    )
  )