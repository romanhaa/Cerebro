tab_gene_set_expression <- tabItem(
    tabName = "geneSetExpression",
    tagList(
      fluidRow(
        column(width = 3, offset = 0, style = "padding:0px;",
          box(
            title = "Input parameters", status = "primary",
            solidHeader = TRUE, width = 12, collapsible = TRUE,
            tagList(
              uiOutput("geneSetexpression_UI"),
              uiOutput("geneSetexpression_scales")
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
                inputId = "geneSetexpression_projection_info",
                label = "info", icon = NULL, class = "btn-xs",
                title = "Show additional information for this panel.",
                style = "margin-right: 5px"
              ),
              actionButton(
                inputId = "geneSetexpression_projection_export",
                label = "export to PDF", icon = NULL, class = "btn-xs",
                title = "Export dimensional reduction to PDF file."
              )
            ),
            status = "primary", solidHeader = TRUE, width = 12,
            collapsible = TRUE,
            tagList(
              scatterD3::scatterD3Output(
                "geneSetexpression_projection", height = "720px"
              ),
              tags$br(),
              htmlOutput("geneSetexpression_genes_displayed")
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
              inputId = "geneSetexpression_by_sample_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          plotly::plotlyOutput("geneSetexpression_by_sample")
        )
      ),
      fluidRow(
        box(
          title = tagList(
            p("Expression levels by cluster",
            style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "geneSetexpression_by_cluster_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          plotly::plotlyOutput("geneSetexpression_by_cluster")
        )
      )
    )
  )