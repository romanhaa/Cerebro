##----------------------------------------------------------------------------##
## Tab: Analysis info.
##----------------------------------------------------------------------------##

tab_analysis_info <- tabItem(
  tabName = "info",
  fluidPage(
    fluidRow(
      column(12,
        titlePanel("Information about samples and analysis"),
        htmlOutput("sample_info_general")
      )
    )
  )
)