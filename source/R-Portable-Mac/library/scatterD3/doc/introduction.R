## ---- include=FALSE------------------------------------------------------
library(scatterD3)

## ----basic, eval=FALSE---------------------------------------------------
#  library(scatterD3)
#  scatterD3(x = mtcars$wt, y = mtcars$mpg)

## ----basic_nse-----------------------------------------------------------
scatterD3(data = mtcars , x = wt, y = mpg)

## ----basic_cust----------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, 
          point_size = 35, point_opacity = 0.5, fixed = TRUE,
          colors = "#A94175")

## ----hover_cust----------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, 
          point_size = 100, point_opacity = 0.5,
          hover_size = 4, hover_opacity = 1)

## ----categorical---------------------------------------------------------
mtcars$cyl_fac <- paste(mtcars$cyl, "cylinders")
scatterD3(data = mtcars, x = cyl_fac, y = mpg)

## ----categorical_left_margin---------------------------------------------
scatterD3(data = mtcars, x = wt, y = cyl_fac, left_margin = 80)

## ----labels--------------------------------------------------------------
mtcars$names <- rownames(mtcars)
scatterD3(data = mtcars, x = wt, y = mpg, lab = names, labels_size = 9)

## ----mapping-------------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, col_var = cyl, symbol_var = gear)

## ----map_size------------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, col_var = cyl, size_var = hp, 
          size_range = c(10,1000), point_opacity = 0.7)

## ----map_custom_colors---------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, col_var = cyl,
          colors = c("4" = "#ECD078", "8" = "#C02942", "6" = "#53777A"))

## ----map_continuous_color------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, col_var = disp)

## ----opacity_var---------------------------------------------------------
scatterD3(data=mtcars, x=mpg, y=wt, opacity_var = drat)

## ----lines---------------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, 
          lines = data.frame(slope = -5.344, intercept = 37.285))

## ----lines_style---------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, 
          lines = data.frame(slope = 0, 
                             intercept = 30,
                             stroke = "red",
                             stroke_width = 5,
                             stroke_dasharray = "10,5"))

## ----lines_default-------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, fixed = TRUE, 
          lines = data.frame(slope = c(0, Inf), 
                             intercept = c(0, 0),
                             stroke = "#000",
                             stroke_width = 1,
                             stroke_dasharray = 5))

## ----log_scales----------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, col_var = cyl,
          x_log = TRUE, y_log = TRUE)

## ----axis_limits---------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, xlim=c(0,10), ylim=c(10,35))

## ----cust_labels---------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, col_var = cyl, symbol_var = gear,
          xlab = "Weight", ylab = "Mpg", col_lab = "Cylinders", symbol_lab = "Gears")

## ----cust_labels_size----------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, col_var = cyl,
          xlab = "Weight", ylab = "Mpg", 
          axes_font_size = "120%",
          legend_font_size = "14px")

## ----cust_left_margin----------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, col_var = cyl,
          left_margin = 80)

## ----caption_character---------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, col_var = cyl,
          caption = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nullam aliquam egestas pretium. Donec auctor semper vestibulum. Phasellus in tempor lacus. Maecenas vehicula, ipsum id malesuada placerat, diam lorem aliquet lectus, non lacinia quam leo quis eros.")

## ----caption_list--------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, col_var = cyl,
          caption = list(title = "Caption title",
                         subtitle = "Caption subtitle",
                         text = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nullam aliquam egestas pretium. Donec auctor semper vestibulum. Phasellus in tempor lacus. Maecenas vehicula, ipsum id malesuada placerat, diam lorem aliquet lectus, non lacinia quam leo quis eros."))

## ----cust_tooltips-------------------------------------------------------
tooltips <- paste("This is an incredible <strong>", rownames(mtcars),"</strong><br />with ", 
                  mtcars$cyl, "cylinders !")
scatterD3(data = mtcars, x = wt, y = mpg, tooltip_text = tooltips)

## ----urls----------------------------------------------------------------
mtcars$urls <- paste0("https://www.duckduckgo.com/?q=", rownames(mtcars))
scatterD3(data = mtcars, x = wt, y = mpg, lab = names, url_var = urls)

## ----click_callback------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg,
   click_callback = "function(id, index) {
   alert('scatterplot ID: ' + id + ' - Point index: ' + index) 
   }")

## ---- click_callback_shiny, eval=FALSE-----------------------------------
#  scatterD3(data = mtcars, x = wt, y = mpg,
#    click_callback = "function(id, index) {
#    if(id && typeof(Shiny) != 'undefined') {
#        Shiny.onInputChange('selected_point', index);
#    }
#  }")

## ----click_callback_shiny_ui, eval = FALSE-------------------------------
#  textOutput("click_selected")

## ----click_callback_shiny_server, eval = FALSE---------------------------
#  output$click_selected <- renderText(paste0("Clicked point : ", input$selected_point))

## ----zoom_callback-------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg,
   zoom_callback = "function(xmin, xmax, ymin, ymax) {
    var zoom = '<strong>Zoom</strong><br />xmin = ' + xmin + '<br />xmax = ' + xmax + '<br />ymin = ' + ymin + '<br />ymax = ' + ymax;
    document.getElementById('zoomExample').innerHTML = zoom;
   }")

## ----ellipses------------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, ellipses = TRUE)

## ----ellipses_col--------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, col_var = cyl, ellipses = TRUE)

## ----nomenu--------------------------------------------------------------
scatterD3(data = mtcars, x = wt, y = mpg, menu = FALSE)

## ----lasso---------------------------------------------------------------
mtcars$names <- rownames(mtcars)
scatterD3(data = mtcars, x = wt, y = mpg, lab = names, lasso = TRUE)

## ----lasso_callback------------------------------------------------------
mtcars$names <- rownames(mtcars)
scatterD3(data = mtcars,
          x = wt, y = mpg, lab = names, 
          lasso = TRUE,
          lasso_callback = "function(sel) {alert(sel.data().map(function(d) {return d.lab}).join('\\n'));}")

## ----labels_export-------------------------------------------------------
mtcars$names <- rownames(mtcars)
scatterD3(data = mtcars, x = wt, y = mpg, lab = names)

## ----labels_export_scatterD3, eval = FALSE-------------------------------
#  labels <- read.csv("scatterD3_labels.csv")
#  scatterD3(data = mtcars, x = wt, y = mpg, lab = names, labels_positions = labels)

## ----labels_export_ggplot2, eval = FALSE---------------------------------
#  labels <- read.csv("scatterD3_labels.csv")
#  library(ggplot2)
#  ggplot() +
#    geom_point(data = mtcars, aes(x=wt, y=mpg)) +
#    geom_text(data = labels,
#              aes(x = lab_x,
#                  y = lab_y,
#                  label = lab))

## ----cust_arrows---------------------------------------------------------
scatterD3(x = c(1, 0.9, 0.7, 0.2, -0.4, -0.5), xlab = "x",
          y = c(1, 0.1, -0.5, 0.5, -0.6, 0.7), ylab = "y",
          lab = LETTERS[1:6], type_var = c("point", rep("arrow", 5)),
          unit_circle = TRUE, fixed = TRUE, 
          xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2))

