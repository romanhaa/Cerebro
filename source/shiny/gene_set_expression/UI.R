##----------------------------------------------------------------------------##
## Panel: Gene set expression
##----------------------------------------------------------------------------##

tab_gene_set_expression <- tabItem(
    tabName = "geneSetExpression",
    tagList(
      fluidRow(
        column(width = 3, offset = 0, style = "padding:0px;",
          box(
            title = "Input parameters", status = "primary",
            solidHeader = TRUE, width = 12, collapsible = TRUE,
            tagList(
              uiOutput("geneSetExpression_UI"),
              uiOutput("geneSetExpression_scales")
            )
          )
        ),
        column(width = 9, offset = 0, style = "padding:0px;",
          box(
            title = tagList(
              p("Dimensional reduction",
                style = "padding-right: 5px; display: inline"
              ),
              actionButton(
                inputId = "geneSetExpression_projection_info",
                label = "info", icon = NULL, class = "btn-xs",
                title = "Show additional information for this panel.",
                style = "margin-right: 5px"
              ),
              actionButton(
                inputId = "geneSetExpression_projection_export",
                label = "export to PDF", icon = NULL, class = "btn-xs",
                title = "Export dimensional reduction to PDF file."
              )
            ),
            status = "primary", solidHeader = TRUE, width = 12,
            collapsible = TRUE,
            tagList(
              scatterD3::scatterD3Output(
                "geneSetExpression_projection", height = "720px"
              ),
              tags$br(),
              htmlOutput("geneSetExpression_genes_displayed")
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
              inputId = "geneSetExpression_by_sample_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          plotly::plotlyOutput("geneSetExpression_by_sample")
        )
      ),
      fluidRow(
        box(
          title = tagList(
            p("Expression levels by cluster",
            style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "geneSetExpression_by_cluster_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          plotly::plotlyOutput("geneSetExpression_by_cluster")
        )
      )
    )
  )