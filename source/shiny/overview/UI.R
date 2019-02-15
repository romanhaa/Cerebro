##----------------------------------------------------------------------------##
## Panel: Overview.
##----------------------------------------------------------------------------##

tab_overview <- tabItem(
    tabName = "overview",
    tagList(
      fluidRow(
        column(width = 3, offset = 0, style = "padding: 0px;",
          box(
            title = "Input parameters", status = "primary",
            solidHeader = TRUE, width = 12, collapsible = TRUE,
            tagList(
              uiOutput("overview_UI"),
              uiOutput("overview_scales")
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
                inputId = "overview_projection_info", label = "info",
                icon = NULL, class = "btn-xs",
                title = "Show additional information for this panel.",
                style = "margin-right: 5px"
              ),
              actionButton(
                inputId = "overview_projection_export",
                label = "export to PDF", icon = NULL, class = "btn-xs",
                title = "Export dimensional reduction to PDF file."
              )
            ),
            status = "primary", solidHeader = TRUE, width = 12,
            collapsible = TRUE,
            # scatterD3::scatterD3Output(
            #   "overview_projection", height = "720px"
            # ),
            plotly::plotlyOutput("overview_projection_plotly", height = "720px")
          )
        )
      )
    )
  )