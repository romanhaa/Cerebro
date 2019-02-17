##----------------------------------------------------------------------------##
## Tab: Overview.
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## UI elements.
##----------------------------------------------------------------------------##
output[["overview_UI"]] <- renderUI({
  tagList(
    selectInput(
      "overview_projection_to_display",
      label = "Projection",
      choices = names(sample_data()$projections)
    ),
    shinyWidgets::pickerInput(
      "overview_samples_to_display",
      label = "Samples to display",
      choices = sample_data()$sample_names,
      selected = sample_data()$sample_names,
      options = list("actions-box" = TRUE),
      multiple = TRUE
    ),
    shinyWidgets::pickerInput(
      "overview_clusters_to_display",
      label = "Clusters to display",
      choices = sample_data()$cluster_names,
      selected = sample_data()$cluster_names,
      options = list("actions-box" = TRUE),
      multiple = TRUE
    ),
    sliderInput(
      "overview_percentage_cells_to_show",
      label = "Show % of cells",
      min = scatter_plot_percentage_cells_to_show[["min"]],
      max = scatter_plot_percentage_cells_to_show[["max"]],
      step = scatter_plot_percentage_cells_to_show[["step"]],
      value = scatter_plot_percentage_cells_to_show[["default"]]
    ),
    selectInput(
      "overview_cell_color",
      label = "Color cells by",
      choices = names(sample_data()$cells)[! names(sample_data()$cells) %in% c("cell_barcode")]
    ),
    selectInput(
      "overview_cell_size_variable",
      label = "Change point size by",
      choices = c("None", "nUMI", "nGene"),
      selected = "None"
    ),
    sliderInput(
      "overview_cell_size_value",
      label = "Point size",
      min = scatter_plot_dot_size[["min"]],
      max = scatter_plot_dot_size[["max"]],
      step = scatter_plot_dot_size[["step"]],
      value = scatter_plot_dot_size[["default"]]
    ),
    checkboxInput(
      "overview_show_ellipses",
      label = "Confidence ellipses",
      value = FALSE
    ),
    sliderInput(
      "overview_cell_opacity",
      label = "Point opacity",
      min = scatter_plot_opacity[["min"]],
      max = scatter_plot_opacity[["max"]],
      step = scatter_plot_opacity[["step"]],
      value = scatter_plot_opacity[["default"]]
    )
  )
})

##----------------------------------------------------------------------------##
## UI elements for X and Y limits in projection.
##----------------------------------------------------------------------------##
output[["overview_scales"]] <- renderUI({
  projection_to_display <- if ( is.null(input[["overview_projection_to_display"]]) || is.na(input[["overview_projection_to_display"]]) ) {
    names(sample_data()$projections)[1] 
  } else {
    input[["overview_projection_to_display"]]
  }
  range_x_min <- round(min(sample_data()$projections[[ projection_to_display ]][,1])*1.1)
  range_x_max <- round(max(sample_data()$projections[[ projection_to_display ]][,1])*1.1)
  range_y_min <- round(min(sample_data()$projections[[ projection_to_display ]][,2])*1.1)
  range_y_max <- round(max(sample_data()$projections[[ projection_to_display ]][,2])*1.1)
  tagList(
    sliderInput(
      "overview_scale_x_manual_range",
      label = "X axis",
      min = range_x_min,
      max = range_x_max,
      value = c(range_x_min, range_x_max)
    ),
    sliderInput(
      "overview_scale_y_manual_range",
      label = "Y axis",
      min = range_y_min,
      max = range_y_max,
      value = c(range_y_min, range_y_max)
    )
  )
})

##----------------------------------------------------------------------------##
## Projection.
##----------------------------------------------------------------------------##
# output$overview_projection <- scatterD3::renderScatterD3({

#   # don't do anything before these inputs are selected
#   req(input[["overview_projection_to_display"]])
#   req(input$overview_samples_to_display)
#   req(input$overview_clusters_to_display)
#   req(input$overview_scale_x_manual_range)

#   # define which projection should be plotted
#   if ( is.null(input[["overview_projection_to_display"]]) || is.na(input[["overview_projection_to_display"]]) ) {
#     projection_to_display <- names(sample_data()$projections)[1]
#   } else {
#     projection_to_display <- input[["overview_projection_to_display"]]
#   }
  
#   # define which samples should be plotted
#   if ( is.null(input$overview_samples_to_display) || is.na(input$overview_samples_to_display) ) {
#     samples_to_display <- sample_data()$samples$overview$sample
#   } else {
#     samples_to_display <- input$overview_samples_to_display
#   }

#   # define which clusters should be plotted
#   if ( is.null(input$overview_clusters_to_display) || is.na(input$overview_clusters_to_display) ) {
#     clusters_to_display <- sample_data()$clusters$overview$cluster
#   } else {
#     clusters_to_display <- input$overview_clusters_to_display
#   }

#   # define which cells should be plotted
#   cells_to_display <- which(
#       grepl(
#         sample_data()$cells$sample,
#         pattern = paste0("^", samples_to_display, "$", collapse="|")
#       ) & 
#       grepl(
#         sample_data()$cells$cluster,
#         pattern = paste0("^", clusters_to_display, "$", collapse="|")
#       )
#     )

#   # randomly remove cells
#   if ( input$overview_percentage_cells_to_show < 100 ) {
#     number_of_cells_to_plot <- ceiling(
#       input$overview_percentage_cells_to_show / 100 * length(cells_to_display)
#     )
#     cells_to_display <- cells_to_display[ sample(1:length(cells_to_display), number_of_cells_to_plot) ]
#   }

#   # extract cells to plot
#   to_plot <- cbind(
#       sample_data()$projections[[ projection_to_display ]][ cells_to_display , ],
#       sample_data()$cells[ cells_to_display , ]
#     )
#   to_plot <- to_plot[ sample(1:nrow(to_plot)) , ]

#   # define variable used to color cells by
#   col_var <- to_plot[ , input[["overview_cell_color"]] ]

#   # define colors
#   if ( is.null(input[["overview_cell_color"]]) || is.na(input[["overview_cell_color"]]) ) {
#     colors <- NULL
#   } else if ( input[["overview_cell_color"]] == "sample" ) {
#     colors <- sample_data()$samples$colors
#   } else if ( input[["overview_cell_color"]] == "cluster" ) {
#     colors <- sample_data()$clusters$colors
#   } else if ( input[["overview_cell_color"]] %in% c("cell_cycle_Regev","cell_cycle_Cyclone") ) {
#     colors <- cell_cycle_colorset
#   } else if ( is.factor(to_plot[,input[["overview_cell_color"]]]) ) {
#     colors <- setNames(colors[1:length(levels(to_plot[,input[["overview_cell_color"]]]))], levels(to_plot[,input[["overview_cell_color"]]]))
#   } else if ( is.character(to_plot[,input[["overview_cell_color"]]]) ) {
#     colors <- colors
#   } else {
#     colors <- NULL
#   }

#   # define variable used for cell size
#   size_var <- if ( input$overview_cell_size_variable == "None" ) NULL else to_plot[ , input$overview_cell_size_variable ]

#   # plot
#   scatterD3::scatterD3(
#     x = to_plot[ , 1 ],
#     y = to_plot[ , 2 ],
#     xlab = colnames(to_plot)[ 1 ],
#     ylab = colnames(to_plot)[ 2 ],
#     xlim = c(
#         input$overview_scale_x_manual_range[1],
#         input$overview_scale_x_manual_range[2]
#       ),
#     ylim = c(
#         input$overview_scale_y_manual_range[1],
#         input$overview_scale_y_manual_range[2]
#       ),
#     point_size = input$overview_cell_size_value,
#     col_var = col_var,
#     col_lab = input[["overview_cell_color"]],
#     colors = colors,
#     ellipses = input$overview_show_ellipses,
#     size_var = size_var,
#     point_opacity = input$overview_cell_opacity,
#     hover_size = 4,
#     hover_opacity = 1,
#     transitions = FALSE,
#     menu = FALSE,
#     tooltip_text = paste0(
#       "<b>Cell</b>: ", to_plot[ , "cell_barcode" ], "<br/>",
#       "<b>Sample</b>: ", to_plot[ , "sample" ], "<br/>",
#       "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br/>",
#       "<b>Transcripts</b>: ", formatC(to_plot[ , "nUMI" ], format = "f", big.mark = ",", digits = 0), "<br/>",
#       "<b>Expressed genes</b>: ", formatC(to_plot[ , "nGene" ], format = "f", big.mark = ",", digits = 0), "<br/>"
#     )
#   )
# })

output[["overview_projection"]] <- plotly::renderPlotly({
  # don't do anything before these inputs are selected
  req(
    input[["overview_projection_to_display"]],
    input[["overview_samples_to_display"]],
    input[["overview_clusters_to_display"]],
    input[["overview_percentage_cells_to_show"]],
    input[["overview_cell_color"]],
    input[["overview_cell_size_variable"]],
    input[["overview_cell_size_value"]],
    input[["overview_cell_opacity"]],
    input[["overview_scale_x_manual_range"]],
    input[["overview_scale_y_manual_range"]]
  )

  projection_to_display <- input[["overview_projection_to_display"]]
  samples_to_display <- input[["overview_samples_to_display"]]
  clusters_to_display <- input[["overview_clusters_to_display"]]
  cells_to_display <- which(
      grepl(
        sample_data()$cells$sample,
        pattern = paste0("^", samples_to_display, "$", collapse="|")
      ) & 
      grepl(
        sample_data()$cells$cluster,
        pattern = paste0("^", clusters_to_display, "$", collapse="|")
      )
    )

  # randomly remove cells
  if ( input[["overview_percentage_cells_to_show"]] < 100 ) {
    number_of_cells_to_plot <- ceiling(
      input[["overview_percentage_cells_to_show"]] / 100 * length(cells_to_display)
    )
    cells_to_display <- cells_to_display[ sample(1:length(cells_to_display), number_of_cells_to_plot) ]
  }

  # extract cells to plot
  to_plot <- cbind(
      sample_data()$projections[[ projection_to_display ]][ cells_to_display , ],
      sample_data()$cells[ cells_to_display , ]
    )
  to_plot <- to_plot[ sample(1:nrow(to_plot)) , ]

  # define variable used to color cells by
  col_var <- to_plot[ , input[["overview_cell_color"]] ]

  # define colors
  colors <- if ( input[["overview_cell_color"]] == "sample" ) {
    sample_data()$samples$colors
  } else if ( input[["overview_cell_color"]] == "cluster" ) {
    sample_data()$clusters$colors
  } else if ( input[["overview_cell_color"]] %in% c("cell_cycle_Regev","cell_cycle_Cyclone") ) {
    cell_cycle_colorset
  } else if ( is.factor(to_plot[ , input[["overview_cell_color"]] ]) ) {
    setNames(colors[1:length(levels(to_plot[ , input[["overview_cell_color"]] ]))], levels(to_plot[ , input[["overview_cell_color"]] ]))
  } else if ( is.character(to_plot[ , input[["overview_cell_color"]] ]) ) {
    colors
  } else {
    NULL
  }

  # define variable used for cell size
  size_var <- if ( input[["overview_cell_size_variable"]] == "None" ) NULL else to_plot[ , input[["overview_cell_size_variable"]] ]

  if ( is.numeric(to_plot[ , input[["overview_cell_color"]] ]) ) {
    plotly::plot_ly(
      to_plot,
      x = ~to_plot[,1],
      y = ~to_plot[,2],
      type = "scattergl",
      mode = "markers",
      marker = list(
        colorbar = list(
          title = colnames(to_plot)[which(colnames(to_plot) == input[["overview_cell_color"]])]
        ),
        color = ~to_plot[ , input[["overview_cell_color"]] ],
        opacity = input[["overview_cell_opacity"]],
        colorscale = "YlGnBu",
        reversescale = TRUE,
        line = list(
          color = "rgb(196,196,196)",
          width = 1
        ),
        size = input[["overview_cell_size_value"]]
      ),
      hoverinfo = "text",
      text = ~paste(
        "<b>Cell</b>: ", to_plot[ , "cell_barcode" ], "<br>",
        "<b>Sample</b>: ", to_plot[ , "sample" ], "<br>",
        "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br>",
        "<b>Transcripts</b>: ", formatC(to_plot[ , "nUMI" ], format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>Expressed genes</b>: ", formatC(to_plot[ , "nGene" ], format = "f", big.mark = ",", digits = 0)
      ),
      height = 720
    ) %>%
    plotly::layout(
      xaxis = list(
        title = colnames(to_plot)[1],
        mirror = TRUE,
        showline = TRUE,
        zeroline = FALSE
      ),
      yaxis = list(
        title = colnames(to_plot)[2],
        mirror = TRUE,
        showline = TRUE,
        zeroline = FALSE
      ),
      hoverlabel = list(font = list(size = 11))
    )
  } else {
    plotly::plot_ly(
      to_plot,
      x = ~to_plot[,1],
      y = ~to_plot[,2],
      color = ~to_plot[ , input[["overview_cell_color"]] ],
      colors = colors,
      type = "scattergl",
      mode = "markers",
      marker = list(
        opacity = input[["overview_cell_opacity"]],
        line = list(
          color = "rgb(196,196,196)",
          width = 1
        ),
        size = input[["overview_cell_size_value"]]
      ),
      hoverinfo = "text",
      text = ~paste(
        "<b>Cell</b>: ", to_plot[ , "cell_barcode" ], "<br>",
        "<b>Sample</b>: ", to_plot[ , "sample" ], "<br>",
        "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br>",
        "<b>Transcripts</b>: ", formatC(to_plot[ , "nUMI" ], format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>Expressed genes</b>: ", formatC(to_plot[ , "nGene" ], format = "f", big.mark = ",", digits = 0)
      ),
      height = 720
    ) %>%
    plotly::layout(
      xaxis = list(
        title = colnames(to_plot)[1],
        mirror = TRUE,
        showline = TRUE,
        zeroline = FALSE
      ),
      yaxis = list(
        title = colnames(to_plot)[2],
        mirror = TRUE,
        showline = TRUE,
        zeroline = FALSE
      ),
      hoverlabel = list(font = list(size = 11))
    )
  }
})

##----------------------------------------------------------------------------##
## Info button.
##----------------------------------------------------------------------------##
observeEvent(input[["overview_projection_info"]], {
  showModal(
    modalDialog(
      overview_projection_info[["text"]],
      title = overview_projection_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Export projection.
##----------------------------------------------------------------------------##
observeEvent(input[["overview_projection_export"]], {
  req(
    input[["overview_projection_to_display"]],
    input[["overview_samples_to_display"]],
    input[["overview_clusters_to_display"]],
    input[["overview_percentage_cells_to_show"]],
    input[["overview_cell_color"]],
    input[["overview_cell_size_variable"]],
    input[["overview_cell_size_value"]],
    input[["overview_cell_opacity"]],
    input[["overview_scale_x_manual_range"]],
    input[["overview_scale_y_manual_range"]]
  )
  library("ggplot2")
  projection_to_display <- input[["overview_projection_to_display"]]
  samples_to_display <- input[["overview_samples_to_display"]]
  clusters_to_display <- input[["overview_clusters_to_display"]]
  cells_to_display <- which(
    grepl(
      sample_data()$cells$sample,
      pattern = paste0("^", samples_to_display, "$", collapse = "|")
    ) &
    grepl(
      sample_data()$cells$cluster,
      pattern = paste0("^", clusters_to_display, "$", collapse = "|")
    )
  )
  to_plot <- cbind(
    sample_data()$projections[[ projection_to_display ]][ cells_to_display , ],
    sample_data()$cells[ cells_to_display , ]
  )
  to_plot <- to_plot[ sample(1:nrow(to_plot)) , ]
  xlim <- c(
    input[["overview_scale_x_manual_range"]][1],
    input[["overview_scale_x_manual_range"]][2]
  )
  ylim <- c(
    input[["overview_scale_y_manual_range"]][1],
    input[["overview_scale_y_manual_range"]][2]
  )
  if ( is.factor(to_plot[ , input[["overview_cell_color"]] ]) || is.character(to_plot[ , input[["overview_cell_color"]] ]) ) {
    if ( input[["overview_cell_color"]] == "sample" ) {
      cols <- sample_data()$samples$colors
    } else if ( input[["overview_cell_color"]] == "cluster" ) {
      cols <- sample_data()$clusters$colors
    } else if ( input[["overview_cell_color"]] %in% c("cell_cycle_Regev","cell_cycle_Cyclone") ) {
      cols <- cell_cycle_colorset
    } else if ( is.factor(to_plot[ , input[["overview_cell_color"]] ]) ) {
      cols <- setNames(colors[1:length(levels(to_plot[ , input[["overview_cell_color"]] ]))], levels(to_plot[ , input[["overview_cell_color"]] ]))
    } else {
      cols <- colors
    }
    p <- ggplot(
        to_plot,
        aes_q(
          x = as.name(colnames(to_plot)[1]),
          y = as.name(colnames(to_plot)[2]),
          fill = as.name(input[["overview_cell_color"]])
        )
      ) +
      geom_point(
        shape = 21,
        size = input[["overview_cell_size_value"]]/3,
        stroke = 0.2,
        color = "#c4c4c4",
        alpha = input[["overview_cell_opacity"]]
      ) +
      scale_fill_manual(values = cols) +
      lims(x = xlim, y = ylim) +
      theme_bw()
  } else {
    p <- ggplot(
        to_plot,
        aes_q(
          x = as.name(colnames(to_plot)[1]),
          y = as.name(colnames(to_plot)[2]),
          fill = as.name(input[["overview_cell_color"]])
        )
      ) +
      geom_point(
        shape = 21,
        size = 3,
        stroke = 0.2,
        color = "#c4c4c4",
        alpha = input[["overview_cell_opacity"]]
      ) +
      scale_fill_distiller(
        palette = "YlGnBu",
        direction = 1,
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
      ) +
      lims(x = xlim, y = ylim) +
      theme_bw()
  }

  out_filename <- paste0(
      plot_export_path, "Cerebro_",
      gsub(
        sample_data()$experiment$experiment_name,
        pattern = " ", replacement = "_"
      ),
      "_overview_", input[["overview_projection_to_display"]], "_by_",
      gsub(
        input[["overview_cell_color"]],
        pattern = "\\.", replacement = "_"
      ),
      ".pdf"
    )

  pdf(NULL)
  ggsave(out_filename, p, height = 8, width = 11)

  if ( file.exists(out_filename) ) {
    shinyWidgets::sendSweetAlert(
      session = session,
      title = "Success!",
      text = paste0("Plot saved successfully as: ", out_filename),
      type = "success"
    )
  } else {
    shinyWidgets::sendSweetAlert(
      session = session,
      title = "Error!",
      text = "Sorry, it seems something went wrong...",
      type = "error"
    )
  }
})