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
  if ( input[["overview_dot_color"]] == "sample" ) {
    colors_this_plot <- reactive_colors()$samples
  } else if ( input[["overview_dot_color"]] == "cluster" ) {
    colors_this_plot <- reactive_colors()$clusters
  } else if ( input[["overview_dot_color"]] %in% c("cell_cycle_seurat","cell_cycle_cyclone","Phase") ) {
    colors_this_plot <- cell_cycle_colorset
  } else if ( is.factor(to_plot[[ input[["overview_dot_color"]] ]]) ) {
    colors_this_plot <- setNames(
      default_colorset[1:length(levels(to_plot[ , input[["overview_dot_color"]] ]))],
      levels(to_plot[ , input[["overview_dot_color"]] ])
    )
  } else if ( is.character(to_plot[[ input[["overview_dot_color"]] ]]) ) {
    colors_this_plot <- setNames(
      default_colorset[1:length(unique(to_plot[ , input[["overview_dot_color"]] ]))],
      unique(to_plot[ , input[["overview_dot_color"]] ])
    )
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
        ),
        source = "overview_projection"
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
        colors = colors_this_plot,
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
        ),
        source = "overview_projection"
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
        ),
        source = "overview_projection"
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
        colors = colors_this_plot,
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
        ),
        source = "overview_projection"
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
          colors_this_plot <- reactive_colors()$samples
        } else if ( input[["overview_dot_color"]] == "cluster" ) {
          colors_this_plot <- reactive_colors()$clusters
        } else if ( input[["overview_dot_color"]] %in% c("cell_cycle_seurat","cell_cycle_cyclone") ) {
          colors_this_plot <- cell_cycle_colorset
        } else if ( is.factor(to_plot[ , input[["overview_dot_color"]] ]) ) {
          colors_this_plot <- setNames(
            default_colorset[1:length(levels(to_plot[ , input[["overview_dot_color"]] ]))],
            levels(to_plot[ , input[["overview_dot_color"]] ])
          )
        } else {
          colors_this_plot <- default_colorset
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
          scale_fill_manual(values = colors_this_plot) +
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

##----------------------------------------------------------------------------##
## Table for details of selected cells.
##----------------------------------------------------------------------------##

output[["overview_details_selected_cells_table"]] <- DT::renderDataTable(server = FALSE, {
  req(
    input[["overview_projection_to_display"]]
  )

  ## extract cells to plot
  to_plot <- cbind(
      sample_data()$projections[[ input[["overview_projection_to_display"]] ]],
      sample_data()$cells
    )

  if ( is.null(plotly::event_data("plotly_selected", source = "overview_projection")) ) {
    table <- tibble(
      cell_barcode = character(),
      sample = character(),
      cluster = character(),
      nUMI = numeric(),
      nGene = numeric(),
      cell_cycle = character()
    )
  } else if ( length(plotly::event_data("plotly_selected", source = "overview_projection")) == 0 ) {
    table <- tibble(
      cell_barcode = character(),
      sample = character(),
      cluster = character(),
      nUMI = numeric(),
      nGene = numeric(),
      cell_cycle = character()
    )
  } else {
    selected_cells <- plotly::event_data("plotly_selected", source = "overview_projection") %>%
      dplyr::mutate(identifier = paste0(x, '-', y))
    table <- to_plot %>%
      dplyr::rename(X1 = 1, X2 = 2) %>%
      dplyr::mutate(identifier = paste0(X1, '-', X2)) %>%
      dplyr::filter(identifier %in% selected_cells$identifier) %>%
      dplyr::select(-c(X1,X2,identifier)) %>%
      dplyr::select(cell_barcode, sample, cluster, nUMI, nGene, everything())
    if ( nrow(table) == 0 ) {
      table <- tibble(
        cell_barcode = character(),
        sample = character(),
        cluster = character(),
        nUMI = numeric(),
        nGene = numeric(),
        cell_cycle = character()
      )
    } else {
      table <- table %>%
        dplyr::mutate(
          nUMI = formattable::comma(nUMI, big.mark = ',', digits = 0),
          nGene = formattable::comma(nGene, big.mark = ',', digits = 0)
        )
    }
  }
  table %>%
  dplyr::rename(
    'Cell barcode' = cell_barcode,
    'Sample' = sample,
    'Cluster' = cluster
  ) %>%
  formattable::formattable(list(
    'nUMI' = formattable::color_tile("white", "orange"),
    'nGene' = formattable::color_tile("white", "orange")
  )) %>%
  formattable::as.datatable(
    filter = "top",
    selection = "none",
    escape = FALSE,
    autoHideNavigation = TRUE,
    rownames = FALSE,
    extensions = c("Buttons"),
    class = "cell-border stripe",
    options = list(
      columnDefs = list(list(visible = FALSE, targets = c(6:ncol(table)-1))),
      dom = "Bfrtip",
      lengthMenu = c(15, 30, 50, 100),
      pageLength = 15,
      buttons = list(
        "colvis",
        list(
          extend = "collection",
          text = "Download",
          buttons = list(
            list(
              extend = "csv",
              filename = "overview_details_of_selected_cells",
              title = "Details of selected cells"
            ),
            list(
              extend = "excel",
              filename = "overview_details_of_selected_cells",
              title = "Details of selected cells"
            )
          )
        )
      )
    )
  ) %>%
  DT::formatStyle(
    columns = c('Sample','Cluster'),
    textAlign = 'center'
  ) %>%
  DT::formatStyle(
    columns = c('nUMI', 'nGene'),
    textAlign = 'right'
  )
})

# info box
observeEvent(input[["overview_details_selected_cells_table_info"]], {
  showModal(
    modalDialog(
      overview_details_selected_cells_table_info$text,
      title = overview_details_selected_cells_table_info$title,
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Plot for selected cells.
## - in sync with selected color variable
##   - if categorical: number of cells in each group
##   - if numerical: box/violin plot
##----------------------------------------------------------------------------##

output[["overview_details_selected_cells_plot"]] <- plotly::renderPlotly({
  req(
    input[["overview_projection_to_display"]],
    input[["overview_dot_color"]]
  )

  ## extract cells to plot
  to_plot <- cbind(
      sample_data()$projections[[ input[["overview_projection_to_display"]] ]],
      sample_data()$cells
    )

  if (
    is.null(plotly::event_data("plotly_selected", source = "overview_projection")) |
    length(plotly::event_data("plotly_selected", source = "overview_projection")) == 0
  ) {
    to_plot <- to_plot %>% dplyr::mutate(group = 'not selected')
  } else {
    selected_cells <- plotly::event_data("plotly_selected", source = "overview_projection") %>%
      dplyr::mutate(identifier = paste0(x, '-', y))
    to_plot <- to_plot %>%
      dplyr::rename(X1 = 1, X2 = 2) %>%
      dplyr::mutate(
        identifier = paste0(X1, '-', X2),
        group = ifelse(identifier %in% selected_cells$identifier, 'selected', 'not selected'),
        group = factor(group, levels = c('selected', 'not selected'))
      )
  }

  color_variable <- input[["overview_dot_color"]]

  ## if the selected coloring variable is categorical, represent the selected
  ## cells in a bar chart
  if (
    is.factor(to_plot[[ color_variable ]]) ||
    is.character(to_plot[[ color_variable ]])
  ) {
    ## calculate number of cells in each group
    t <- to_plot %>%
      dplyr::filter(group == 'selected') %>%
      dplyr::select(!!! rlang::syms(color_variable)) %>%
      dplyr::group_by_at(1) %>%
      dplyr::tally() %>%
      dplyr::ungroup()

    ## convert factor to character to avoid empty bars when selecting cells of
    ## certain groups
    t[[1]] <- as.character(t[[1]])

    ## create color assignment for groups
    if ( color_variable == 'sample' ) {
      colors_this_plot <- reactive_colors()$samples
    } else if ( color_variable == 'cluster' ) {
      colors_this_plot <- reactive_colors()$clusters
    } else {
      colors_this_plot <- setNames(
        default_colorset[1:length(t[[ 1 ]])],
        t[[ 1 ]]
      )
    }

    ## make bar chart
    plotly::plot_ly(
      t,
      x = ~t[[1]],
      y = ~t[[2]],
      type = "bar",
      color = ~t[[1]],
      colors = colors_this_plot,
      source = "subset",
      showlegend = FALSE,
      hoverinfo = "y"
    ) %>%
    plotly::layout(
      title = "",
      xaxis = list(
        title = "",
        mirror = TRUE,
        showline = TRUE
      ),
      yaxis = list(
        title = "Number of cells",
        hoverformat = ".0f",
        mirror = TRUE,
        showline = TRUE
      ),
      dragmode = "select",
      hovermode = "compare"
    )

  ## if the selected coloring variable is not categorical but continuous
  } else {
    ## remove unnecessary columns
    t <- to_plot %>%
      dplyr::select(group, !!! rlang::syms(color_variable))

    ## create violin/box plot
    plotly::plot_ly(
      t,
      x = ~t[[1]],
      y = ~t[[2]],
      type = "violin",
      box = list(
        visible = TRUE
      ),
      meanline = list(
        visible = TRUE
      ),
      color = ~t[[1]],
      colors = setNames(c('#e74c3c','#7f8c8d'),c('selected', 'not selected')),
      source = "subset",
      showlegend = FALSE,
      hoverinfo = "y",
      marker = list(
        size = 5
      )
    ) %>%
    plotly::layout(
      title = "",
      xaxis = list(
        title = "",
        mirror = TRUE,
        showline = TRUE
      ),
      yaxis = list(
        title = colnames(t)[2],
        hoverformat = ".0f",
        mirror = TRUE,
        showline = TRUE
      ),
      dragmode = "select",
      hovermode = "compare"
    )
  }
})

# info box
observeEvent(input[["overview_details_selected_cells_plot_info"]], {
  showModal(
    modalDialog(
      overview_details_selected_cells_plot_info$text,
      title = overview_details_selected_cells_plot_info$title,
      easyClose = TRUE,
      footer = NULL
    )
  )
})
