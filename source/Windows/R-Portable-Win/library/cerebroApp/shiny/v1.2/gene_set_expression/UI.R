##----------------------------------------------------------------------------##
## Tab: Gene set expression
##----------------------------------------------------------------------------##

tab_gene_set_expression <- tabItem(
  tabName = "geneSetExpression",
  tagList(
    fluidRow(
      column(width = 3, offset = 0, style = "padding:0px;",
        cerebroBox(
          title = "Input parameters",
          tagList(
            uiOutput("geneSetExpression_UI"),
            uiOutput("geneSetExpression_color_scale_range"),
            uiOutput("geneSetExpression_scales")
          )
        )
      ),
      column(width = 9, offset = 0, style = "padding:0px;",
        cerebroBox(
          title = tagList(
            boxTitle("Dimensional reduction"),
            actionButton(
              inputId = "geneSetExpression_projection_info",
              label = "info",
              icon = NULL,
              class = "btn-xs",
              title = "Show additional information for this panel.",
              style = "margin-right: 5px"
            ),
            actionButton(
              inputId = "geneSetExpression_projection_export",
              label = "export to PDF",
              icon = NULL,
              class = "btn-xs",
              title = "Export dimensional reduction to PDF file."
            )
          ),
          tagList(
            plotly::plotlyOutput(
              "geneSetExpression_projection",
              width = "auto",
              height = "85vh"
            ),
            tags$br(),
            htmlOutput("geneSetExpression_genes_displayed")
          )
        )
      )
    ),
    fluidRow(
      cerebroBox(
        title = tagList(
          boxTitle("Details of selected cells"),
          cerebroInfoButton("geneSetExpression_details_selected_cells_info")
        ),
        DT::dataTableOutput("geneSetExpression_details_selected_cells")
      )
    ),
    fluidRow(
      cerebroBox(
        title = tagList(
          boxTitle("Expression levels in selected cells"),
          cerebroInfoButton("geneSetExpression_in_selected_cells_info")
        ),
        plotly::plotlyOutput("geneSetExpression_in_selected_cells")
      )
    ),
    fluidRow(
      cerebroBox(
        title = tagList(
          boxTitle("Expression levels by sample"),
          cerebroInfoButton("geneSetExpression_by_sample_info")
        ),
        plotly::plotlyOutput("geneSetExpression_by_sample")
      )
    ),
    fluidRow(
      cerebroBox(
        title = tagList(
          boxTitle("Expression levels by cluster"),
          cerebroInfoButton("geneSetExpression_by_cluster_info")
        ),
        plotly::plotlyOutput("geneSetExpression_by_cluster")
      )
    ),
    fluidRow(
      cerebroBox(
        title = tagList(
          boxTitle("Expression levels by gene"),
          cerebroInfoButton("geneSetExpression_by_gene_info")
        ),
        plotly::plotlyOutput("geneSetExpression_by_gene")
      )
    )
  )
)
