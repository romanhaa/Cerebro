##----------------------------------------------------------------------------##
## Tab: Overview.
##----------------------------------------------------------------------------##

tab_overview <- tabItem(
  tabName = "overview",
  tagList(
    fluidRow(
      column(width = 3, offset = 0, style = "padding: 0px;",
        cerebroBox(
          title = "Input parameters",
          tagList(
            uiOutput("overview_UI"),
            uiOutput("overview_scales")
          )
        )
      ),
      column(width = 9, offset = 0, style = "padding: 0px;",
        cerebroBox(
          title = tagList(
            boxTitle("Dimensional reduction"),
            actionButton(
              inputId = "overview_projection_info",
              label = "info",
              icon = NULL,
              class = "btn-xs",
              title = "Show additional information for this panel.",
              style = "margin-right: 5px"
            ),
            actionButton(
              inputId = "overview_projection_export",
              label = "export to PDF",
              icon = NULL,
              class = "btn-xs",
              title = "Export dimensional reduction to PDF file."
            )
          ),
          plotly::plotlyOutput(
            "overview_projection",
            width = "auto",
            height = "85vh"
          )
        )
      )
    )
  )
)
