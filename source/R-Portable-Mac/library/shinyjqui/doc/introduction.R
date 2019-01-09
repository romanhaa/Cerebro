## ---- include=FALSE------------------------------------------------------
library(shiny)
library(shinyjqui)

## ---- eval = FALSE-------------------------------------------------------
#  # create a draggable textInput in shiny ui
#  ui <- fluidPage(
#    jqui_draggable(textInput('foo', 'Input'))
#  )

## ---- eval = FALSE-------------------------------------------------------
#  # create a textInput in shiny ui
#  ui <- fluidPage(
#    textInput('foo', 'Input')
#  )
#  
#  # make the textInput with id "foo" draggable in shiny server
#  server <- function(input, output) {
#    jqui_draggable(ui = '#foo')
#  }

## ---- eval = FALSE-------------------------------------------------------
#  # in shiny ui, make each element in the tagList draggable
#  ui <- fluidPage(
#    jqui_draggable(
#      tagList(
#        selectInput('sel', 'Select', choices = month.abb),
#        checkboxGroupInput('chbox', 'Checkbox', choices = month.abb),
#        plotOutput('plot', width = '400px', height = '400px')
#      )
#    )
#  )

## ---- eval = FALSE-------------------------------------------------------
#  server <- function(input, output) {
#    jqui_draggable("#sel,#chbox,#plot")
#  }

## ---- echo=FALSE---------------------------------------------------------
draggable_shiny <- data.frame(
  Interaction_type = 'draggable',
  suffix = c('position', 'is_dragging'),
  `The_returned_shiny_input_value` = c(
    'A list of the element\'s left and top distances (in pixels) to its parent element',
    'TRUE or FALSE that indicate whether the element is dragging'
  )
)

droppable_shiny <- data.frame(
  Interaction_type = 'droppable',
  suffix = c('dragging', 'over', 'drop', 'dropped', 'out'),
  `The_returned_shiny_input_value` = c(
    'The id of an acceptable element that is now dragging',
    'The id of the last acceptable element that is dragged over',
    'The id of the last acceptable element that is dropped',
    'The ids of all acceptable elements that is currently dropped',
    'The id of the last acceptable element that is dragged out'
  )
)

resizable_shiny <- data.frame(
  Interaction_type = 'resizable',
  suffix = c('size', 'is_resizing'),
  `The_returned_shiny_input_value` = c(
    'A list of the element\'s current size',
    'TRUE or FALSE that indicate whether the element is resizing'
  )
)

selectable_shiny <- data.frame(
  Interaction_type = 'selectable',
  suffix = c('selected', 'is_selecting'),
  `The_returned_shiny_input_value` = c(
    'A dataframe containing the id and innerText of curently selected items',
    'TRUE or FALSE that indicate whether the element is selecting (e.g. during lasso selection)'
  )
)

sortable_shiny <- data.frame(
  Interaction_type = 'sortable',
  suffix = c('order'),
  `The_returned_shiny_input_value` = c(
    'A dataframe containing the id and innerText of items in the current order'
  )
)

knitr::kable(rbind(draggable_shiny, droppable_shiny, resizable_shiny, 
                   selectable_shiny, sortable_shiny))


## ---- eval = FALSE-------------------------------------------------------
#  shiny_opt = list(
#  
#    # define shiny input value input$id_suffix1
#    suffix1 = list(
#      # on event_type1 run callback1 and send the returned value to input$id_suffix1
#      event_type1 = JS(callback1),
#      # on event_type2 or event_type3 run callback2 and send the returned value to input$id_suffix1
#      `event_type2 event_type3` = JS(callback2),
#      ...
#    ),
#  
#    # define another shiny input value input$id_suffix2
#    suffix2 = list(
#      ...
#    ),
#  
#    # define more shiny input values
#  
#  )
#  
#  # pass the shiny option to draggable
#  jqui_draggable('#foo', options = list(
#    shiny = shiny_opt,
#    #other draggable-specific options
#  ))

## ---- eval = FALSE-------------------------------------------------------
#  # server
#  jqui_draggable('#foo', options = list(
#    shiny = list(
#      # By default, draggable element has a shiny input value showing the
#      # element's position (relative to the parent element). Here, another shiny
#      # input value (input$foo_offset) is added. It gives the element's offset
#      # (position relative to the document).
#      offset = list(
#        # return the updated offset value when the draggable is created or dragging
#        `dragcreate drag` = JS('function(event, ui) {return $(event.target).offset();}'),
#      )
#    )
#  ))

## ---- eval = FALSE-------------------------------------------------------
#  # drag only horizontally
#  jqui_draggable('#foo', options = list(axis = 'x'))
#  # make movement snapping to a 80 x 80 grid
#  jqui_draggable('#foo', options = list(grid = c(80, 80)))

## ---- eval = FALSE-------------------------------------------------------
#  jqui_droppable('#foo', options = list(
#    accept = '#bar', # jQuery selector to define which draggable element to monitor. Accept anything if not set.
#    classes = list(
#      `ui-droppable-active` = 'ui-state-focus', # change class when draggable element is dragging
#      `ui-droppable-hover` = 'ui-state-highlight' # change class when draggable element is dragging over
#    ),
#    drop = JS(
#      'function(event, ui){$(this).addClass("ui-state-active");}'
#    ) # a javascrip callback to change class when draggable element is dropped in
#  ))

## ---- eval = FALSE-------------------------------------------------------
#  # keep aspect ratio when resizing
#  jqui_resizable('#foo', options = list(aspectRatio = TRUE))
#  
#  # Limit the resizable element to a maximum or minimum height or width
#  jqui_resizable('#foo', options = list(minHeight = 100, maxHeight = 300,
#                                        minWidth = 200, maxWidth = 400))
#  
#  # make the two plotOutputs resize synchronously
#  jqui_resizable(plotOutput('plot1', width = '400px', height = '400px'),
#                    options = list(alsoResize = '#plot2')),
#  plotOutput('plot2', width = '400px', height = '400px')

## ---- eval = FALSE-------------------------------------------------------
#  # highlight the selected plotOutput
#  jqui_selectable(
#    div(
#      plotOutput('plot1', width = '400px', height = '400px'),
#      plotOutput('plot2', width = '400px', height = '400px')
#    ),
#    options = list(classes = list(`ui-selected` = 'ui-state-highlight'))
#  )

## ---- eval = FALSE-------------------------------------------------------
#  # change opacity while sorting
#  jqui_sortable('#foo', options = list(opacity = 0.5))
#  
#  # only items with class "items" inside the element become sortable
#  jqui_sortable('#foo', options = list(items = '> .items'))
#  
#  # connect two sortable elements, so that items in one element can be dragged to another
#  jqui_sortable('#foo1', options = list(connectWith = '#foo2'))
#  jqui_sortable('#foo2', options = list(connectWith = '#foo1'))
#  

## ---- echo=FALSE---------------------------------------------------------
get_jqui_effects()

## ---- echo=FALSE---------------------------------------------------------
func_intro <- data.frame(Functions = c('jqui_effect', 'jqui_show', 'jqui_hide', 'jqui_toggle'), 
                         Description = c('Let element(s) to show an animation immediately.',
                                         'Display hidden element(s) with an animation',
                                         'Hide element(s) with an animation',
                                         'Display or hide element(s) with an animation'),
                         Where_to_use = rep('server', times = 4),
                         stringsAsFactors = FALSE)
knitr::kable(func_intro, row.names = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  # ui
#  plotOutput('foo', width = '400px', height = '400px')
#  
#  # server
#  jqui_effect('#foo', effect = 'bounce') # bounces the plot
#  jqui_effect('#foo', effect = 'scale', options = list(percent = 50)) # scale to 50%
#  jqui_hide('#foo', effect = 'size', options = list(width = 200, height = 60)) # resize then hide
#  jqui_show('#foo', effect = 'clip') # show the plot by clipping

## ---- echo=FALSE---------------------------------------------------------
func_intro <- data.frame(Functions = c('jqui_add_class', 'jqui_remove_class', 'jqui_switch_class'), 
                         Description = c('Add class(es) to element(s) while animating all style changes.',
                                         'Remove class(es) from element(s) while animating all style changes.',
                                         'Add and remove class(es) to element(s) while animating all style changes.'),
                         Where_to_use = rep('server', times = 3),
                         stringsAsFactors = FALSE)
knitr::kable(func_intro, row.names = FALSE)

