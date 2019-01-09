## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE, results='asis'-----------------------------------------
tbl <- data.frame(
  Interactions = c("draggable", "resizable", "selectable", "sortable"),
  State = c(
    "The position of the draggable element",
    "The dimension of the resizable element",
    "The selected items inside the selectable element",
    "The order of items inside the sortable element"
  )
)
knitr::kable(tbl)

## ---- eval=FALSE---------------------------------------------------------
#  ui <- fluidPage(
#    actionButton("save", "Save position"),
#    actionButton("restore", "Restore position"),
#    # create a draggable textInput
#    jqui_draggable(textInput("foo", "Textinput"))
#  )
#  
#  server <- function(input, output) {
#    # on save button clicked, save the current postion of the textInput
#    observeEvent(input$save, {
#      jqui_draggable("#foo", operation = "save")
#    })
#    # on restore button clicked, move the textInput back to the last saved position
#    observeEvent(input$restore, {
#      jqui_draggable("#foo", operation = "load")
#    })
#  }
#  
#  shinyApp(ui, server)

## ---- eval=FALSE---------------------------------------------------------
#  ui <- fluidPage(
#    actionButton("save", "Save order"),
#    actionButton("restore", "Restore order"),
#    orderInput("foo1", label = NULL, items = 1:3, connect = "foo2"),
#    orderInput("foo2", label = NULL, items = NULL, placeholder = "empty")
#  )
#  
#  server <- function(input, output) {
#    observeEvent(input$save, {
#      jqui_sortable("#foo1,#foo2", operation = "save")
#    })
#    observeEvent(input$restore, {
#      jqui_sortable("#foo1,#foo2", operation = "load")
#    })
#  }
#  
#  shinyApp(ui, server)

## ---- eval=FALSE---------------------------------------------------------
#  
#  ui <- function(request) {
#    fluidPage(
#      bookmarkButton(),
#      jqui_resizable(plotOutput('gg', width = '200px', height = '200px'))
#    )
#  }
#  
#  server <- function(input, output) {
#    output$gg <- renderPlot({
#      ggplot(mtcars, aes(x = cyl, y = mpg)) + geom_point()
#    })
#    # enable interaction state bookmarking
#    jqui_bookmarking()
#  }
#  
#  enableBookmarking(store = "url")
#  
#  shinyApp(ui, server)

