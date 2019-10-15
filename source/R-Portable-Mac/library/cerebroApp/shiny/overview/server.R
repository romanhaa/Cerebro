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
      "overview_dot_color",
      label = "Color cells by",
      choices = names(sample_data()$cells)[! names(sample_data()$cells) %in% c("cell_barcode")]
    ),
    sliderInput(
      "overview_dot_size",
      label = "Dot size",
      min = scatter_plot_dot_size[["min"]],
      max = scatter_plot_dot_size[["max"]],
      step = scatter_plot_dot_size[["step"]],
      value = scatter_plot_dot_size[["default"]]
    ),
    sliderInput(
      "overview_dot_opacity",
      label = "Dot opacity",
      min = scatter_plot_dot_opacity[["min"]],
      max = scatter_plot_dot_opacity[["max"]],
      step = scatter_plot_dot_opacity[["step"]],
      value = scatter_plot_dot_opacity[["default"]]
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
  range_x_min <- sample_data()$projections[[ projection_to_display ]][,1] %>% min() %>% "*"(ifelse(.<0, 1.1, 0.9)) %>% round()
  range_x_max <- sample_data()$projections[[ projection_to_display ]][,1] %>% max() %>% "*"(ifelse(.<0, 0.9, 1.1)) %>% round()
  range_y_min <- sample_data()$projections[[ projection_to_display ]][,2] %>% min() %>% "*"(ifelse(.<0, 1.1, 0.9)) %>% round()
  range_y_max <- sample_data()$projections[[ projection_to_display ]][,2] %>% max() %>% "*"(ifelse(.<0, 0.9, 1.1)) %>% round()
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

output[["overview_projection"]] <- plotly::renderPlotly({
  # don't do anything before these inputs are selected
  req(
    input[["overview_projection_to_display"]],
    input[["overview_samples_to_display"]],
    input[["overview_clusters_to_display"]],
    input[["overview_percentage_cells_to_show"]],
    input[["overview_dot_color"]],
    input[["overview_dot_size"]],
    input[["overview_dot_opacity"]],
    input[["overview_scale_x_manual_range"]],
    input[["overview_scale_y_manual_range"]]
  )

  projection_to_display <- input[["overview_projection_to_display"]]
  samples_to_display <- input[["overview_samples_to_display"]]
  clusters_to_display <- input[["overview_clusters_to_display"]]
  cells_to_display <- which(
      (sample_data()$cells$sample %in% samples_to_display) &
      (sample_data()$cells$cluster %in% clusters_to_display)
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

  # define colors
  colors <- if ( input[["overview_dot_color"]] == "sample" ) {
    sample_data()$samples$colors
  } else if ( input[["overview_dot_color"]] == "cluster" ) {
    sample_data()$clusters$colors
  } else if ( input[["overview_dot_color"]] %in% c("cell_cycle_seurat","cell_cycle_cyclone") ) {
    cell_cycle_colorset
  } else if ( is.factor(to_plot[[ input[["overview_dot_color"]] ]]) ) {
    setNames(colors[1:length(levels(to_plot[ , input[["overview_dot_color"]] ]))], levels(to_plot[ , input[["overview_dot_color"]] ]))
  } else if ( is.character(to_plot[[ input[["overview_dot_color"]] ]]) ) {
    colors
  } else {
    NULL
  }

  # define variable used for cell size
  # size_var <- if ( input[["overview_cell_size_variable"]] == "None" ) NULL else to_plot[ , input[["overview_cell_size_variable"]] ]

  if ( ncol(sample_data()$projections[[ projection_to_display ]]) == 3 ) {
    if ( is.numeric(to_plot[[ input[["overview_dot_color"]] ]]) ) {
      plotly::plot_ly(
        to_plot,
        x = ~to_plot[,1],
        y = ~to_plot[,2],
        z = ~to_plot[,3],
        type = "scatter3d",
        mode = "markers",
        marker = list(
          colorbar = list(
            title = input[["overview_dot_color"]]
          ),
          color = ~to_plot[ , input[["overview_dot_color"]] ],
          opacity = input[["overview_dot_opacity"]],
          colorscale = "YlGnBu",
          reversescale = TRUE,
          line = list(
            color = "rgb(196,196,196)",
            width = 1
          ),
          size = input[["overview_dot_size"]]
        ),
        hoverinfo = "text",
        text = ~paste(
          "<b>Cell</b>: ", to_plot[ , "cell_barcode" ], "<br>",
          "<b>Sample</b>: ", to_plot[ , "sample" ], "<br>",
          "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br>",
          "<b>Transcripts</b>: ", formatC(to_plot[ , "nUMI" ], format = "f", big.mark = ",", digits = 0), "<br>",
          "<b>Expressed genes</b>: ", formatC(to_plot[ , "nGene" ], format = "f", big.mark = ",", digits = 0)
        )
      ) %>%
      plotly::layout(
        scene = list(
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
          zaxis = list(
            title = colnames(to_plot)[3],
            mirror = TRUE,
            showline = TRUE,
            zeroline = FALSE
          )
        ),
        hoverlabel = list(
          font = list(
            size = 11
          )
        )
      )
    } else {
      plotly::plot_ly(
        to_plot,
        x = ~to_plot[,1],
        y = ~to_plot[,2],
        z = ~to_plot[,3],
        color = ~to_plot[ , input[["overview_dot_color"]] ],
        colors = colors,
        type = "scatter3d",
        mode = "markers",
        marker = list(
          opacity = input[["overview_dot_opacity"]],
          line = list(
            color = "rgb(196,196,196)",
            width = 1
          ),
          size = input[["overview_dot_size"]]
        ),
        hoverinfo = "text",
        text = ~paste(
          "<b>Cell</b>: ", to_plot[ , "cell_barcode" ], "<br>",
          "<b>Sample</b>: ", to_plot[ , "sample" ], "<br>",
          "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br>",
          "<b>Transcripts</b>: ", formatC(to_plot[ , "nUMI" ], format = "f", big.mark = ",", digits = 0), "<br>",
          "<b>Expressed genes</b>: ", formatC(to_plot[ , "nGene" ], format = "f", big.mark = ",", digits = 0)
        )
      ) %>%
      plotly::layout(
        scene = list(
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
          zaxis = list(
            title = colnames(to_plot)[3],
            mirror = TRUE,
            showline = TRUE,
            zeroline = FALSE
          )
        ),
        hoverlabel = list(font = list(size = 11))
      )
    }
  } else {
    if ( is.numeric(to_plot[[ input[["overview_dot_color"]] ]]) ) {
      plot <- plotly::plot_ly(
        to_plot,
        x = ~to_plot[,1],
        y = ~to_plot[,2],
        type = "scatter",
        mode = "markers",
        marker = list(
          colorbar = list(
            title = input[["overview_dot_color"]]
          ),
          color = ~to_plot[ , input[["overview_dot_color"]] ],
          opacity = input[["overview_dot_opacity"]],
          colorscale = "YlGnBu",
          reversescale = TRUE,
          line = list(
            color = "rgb(196,196,196)",
            width = 1
          ),
          size = input[["overview_dot_size"]]
        ),
        hoverinfo = "text",
        text = ~paste(
          "<b>Cell</b>: ", to_plot[ , "cell_barcode" ], "<br>",
          "<b>Sample</b>: ", to_plot[ , "sample" ], "<br>",
          "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br>",
          "<b>Transcripts</b>: ", formatC(to_plot[ , "nUMI" ], format = "f", big.mark = ",", digits = 0), "<br>",
          "<b>Expressed genes</b>: ", formatC(to_plot[ , "nGene" ], format = "f", big.mark = ",", digits = 0)
        )
      ) %>%
      plotly::layout(
        xaxis = list(
          title = colnames(to_plot)[1],
          mirror = TRUE,
          showline = TRUE,
          zeroline = FALSE,
          range = input[["overview_scale_x_manual_range"]]
        ),
        yaxis = list(
          title = colnames(to_plot)[2],
          mirror = TRUE,
          showline = TRUE,
          zeroline = FALSE,
          range = input[["overview_scale_y_manual_range"]]
        ),
        hoverlabel = list(font = list(size = 11))
      )
      if ( preferences[["use_webgl"]] == TRUE ) {
        plot %>% plotly::toWebGL()
      } else {
        plot
      }
    } else {
      plot <- plotly::plot_ly(
        to_plot,
        x = ~to_plot[,1],
        y = ~to_plot[,2],
        color = ~to_plot[[ input[["overview_dot_color"]] ]],
        colors = colors,
        type = "scatter",
        mode = "markers",
        marker = list(
          opacity = input[["overview_dot_opacity"]],
          line = list(
            color = "rgb(196,196,196)",
            width = 1
          ),
          size = input[["overview_dot_size"]]
        ),
        hoverinfo = "text",
        text = ~paste(
          "<b>Cell</b>: ", to_plot[[ "cell_barcode" ]], "<br>",
          "<b>Sample</b>: ", to_plot[[ "sample" ]], "<br>",
          "<b>Cluster</b>: ", to_plot[[ "cluster" ]], "<br>",
          "<b>Transcripts</b>: ", formatC(to_plot[[ "nUMI" ]], format = "f", big.mark = ",", digits = 0), "<br>",
          "<b>Expressed genes</b>: ", formatC(to_plot[[ "nGene" ]], format = "f", big.mark = ",", digits = 0)
        )
      ) %>%
      plotly::layout(
        xaxis = list(
          title = colnames(to_plot)[1],
          mirror = TRUE,
          showline = TRUE,
          zeroline = FALSE,
          range = input[["overview_scale_x_manual_range"]]
        ),
        yaxis = list(
          title = colnames(to_plot)[2],
          mirror = TRUE,
          showline = TRUE,
          zeroline = FALSE,
          range = input[["overview_scale_y_manual_range"]]
        ),
        hoverlabel = list(font = list(size = 11))
      )
      if ( preferences[["use_webgl"]] == TRUE ) {
        plot %>% plotly::toWebGL()
      } else {
        plot
      }
    }
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
    input[["overview_dot_color"]],
    input[["overview_dot_size"]],
    input[["overview_dot_opacity"]],
    input[["overview_scale_x_manual_range"]],
    input[["overview_scale_y_manual_range"]]
  )
  library("ggplot2")
  if ( exists("plot_export_path") ) {
    projection_to_display <- input[["overview_projection_to_display"]]
    samples_to_display <- input[["overview_samples_to_display"]]
    clusters_to_display <- input[["overview_clusters_to_display"]]
    cells_to_display <- which(
        (sample_data()$cells$sample %in% samples_to_display) &
        (sample_data()$cells$cluster %in% clusters_to_display)
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

    if ( ncol(sample_data()$projections[[ projection_to_display ]]) == 3 ) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Sorry!",
        text = "It's currently not possible to create PDF plots from 3D dimensional reductions. Please use the PNG export button in the panel or a 2D dimensional reduction instead.",
        type = "error"
      )
    } else {
      if ( is.factor(to_plot[ , input[["overview_dot_color"]] ]) || is.character(to_plot[ , input[["overview_dot_color"]] ]) ) {
        if ( input[["overview_dot_color"]] == "sample" ) {
          if ( !is.null(sample_data()$samples$colors) ) {
            cols <- sample_data()$samples$colors
          } else {
            cols <- colors
          }
        } else if ( input[["overview_dot_color"]] == "cluster" ) {
          if ( !is.null(sample_data()$clusters$colors) ) {
            cols <- sample_data()$clusters$colors
          } else {
            cols <- colors
          }
        } else if ( input[["overview_dot_color"]] %in% c("cell_cycle_seurat","cell_cycle_cyclone") ) {
          cols <- cell_cycle_colorset
        } else if ( is.factor(to_plot[ , input[["overview_dot_color"]] ]) ) {
          cols <- setNames(colors[1:length(levels(to_plot[ , input[["overview_dot_color"]] ]))], levels(to_plot[ , input[["overview_dot_color"]] ]))
        } else {
          cols <- colors
        }
        p <- ggplot(
            to_plot,
            aes_q(
              x = as.name(colnames(to_plot)[1]),
              y = as.name(colnames(to_plot)[2]),
              fill = as.name(input[["overview_dot_color"]])
            )
          ) +
          geom_point(
            shape = 21,
            size = input[["overview_dot_size"]]/3,
            stroke = 0.2,
            color = "#c4c4c4",
            alpha = input[["overview_dot_opacity"]]
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
              fill = as.name(input[["overview_dot_color"]])
            )
          ) +
          geom_point(
            shape = 21,
            size = input[["overview_dot_size"]]/3,
            stroke = 0.2,
            color = "#c4c4c4",
            alpha = input[["overview_dot_opacity"]]
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
            input[["overview_dot_color"]],
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
    }
  } else {
    shinyWidgets::sendSweetAlert(
      session = session,
      title = "Error!",
      text = "Sorry, we couldn't find a place to store the figure. Please submit an issue on GitHub @ https://github.com/romanhaa/cerebroApp",
      type = "error"
    )
  }
})