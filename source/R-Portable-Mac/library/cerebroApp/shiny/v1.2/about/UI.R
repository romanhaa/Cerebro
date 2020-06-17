##----------------------------------------------------------------------------##
## Tab: About.
##----------------------------------------------------------------------------##

tab_about <- tabItem(
  tabName = "about",
  fluidPage(
    fluidRow(
      column(12,
        titlePanel("About Cerebro")
      ),
      column(8,
        htmlOutput("about"),
        uiOutput("preferences"),
        actionButton("browser", "browser"),
        tags$script("$('#browser').hide();")
      ),
      column(4,
        imageOutput("logo_Cerebro")
      )
    )
  )
)
