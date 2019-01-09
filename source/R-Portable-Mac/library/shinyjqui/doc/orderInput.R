## ---- include=FALSE------------------------------------------------------
library(shiny)
library(shinyjqui)

## ---- eval=FALSE---------------------------------------------------------
#  server <- function(input, output) {
#    output$order <- renderPrint({input$foo_order})
#  }
#  
#  ui <- fluidPage(
#    orderInput(inputId = 'foo', label = 'A simple example', items = c('A', 'B', 'C')),
#    verbatimTextOutput('order')
#  )
#  
#  shinyApp(ui, server)

## ---- eval=FALSE---------------------------------------------------------
#  # items in A can be dragged to B
#  orderInput('A', 'A', items = 1:3, connect = 'B')
#  # items in B can be dragged to A
#  orderInput('B', 'B', items = 4:6, connect = 'A')

## ---- eval=FALSE---------------------------------------------------------
#  # In source mode, items dragged to B are copied
#  orderInput('A', 'A', items = 1:3, connect = 'B', as_source = TRUE)
#  orderInput('B', 'B', items = 4:6)

## ---- eval=FALSE---------------------------------------------------------
#  orderInput('A', 'A', items = 1:3, connect = 'B')
#  orderInput('B', 'B', items = NULL, placeholder = 'Drag item here...')

## ---- eval=FALSE---------------------------------------------------------
#  orderInput('default', 'default', items = 1:3, item_class = 'default')
#  orderInput('primary', 'primary', items = 1:3, item_class = 'primary')
#  orderInput('success', 'success', items = 1:3, item_class = 'success')
#  orderInput('info', 'info', items = 1:3, item_class = 'info')
#  orderInput('warning', 'warning', items = 1:3, item_class = 'warning')
#  orderInput('danger', 'danger', items = 1:3, item_class = 'danger')

