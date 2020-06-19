##----------------------------------------------------------------------------##
## Tab: Trajectory.
##----------------------------------------------------------------------------##

# what needs to be done
# - let user show gene expression as well? probably more complicated

# UI element: display results or alternative text
output[["trajectory_UI"]] <- renderUI({
  if ( length(sample_data()$trajectory$monocle2) > 0 ) {
    tagList(
      fluidRow(
        column(width = 3, offset = 0, style = "padding: 0px;",
          cerebroBox(
            title = "Input parameters",
            tagList(
              uiOutput("trajectory_input")
            )
          )
        ),
        column(width = 9, offset = 0, style = "padding: 0px;",
          cerebroBox(
            title = tagList(
              boxTitle("Trajectory"),
              actionButton(
                inputId = "trajectory_projection_info",
                label = "info",
                icon = NULL,
                class = "btn-xs",
                title = "Show additional information for this panel.",
                style = "margin-right: 5px"
              ),
              actionButton(
                inputId = "trajectory_projection_export",
                label = "export to PDF",
                icon = NULL,
                class = "btn-xs",
                title = "Export trajectory to PDF file."
              )
            ),
            plotly::plotlyOutput(
              "trajectory_projection",
              width = "auto",
              height = "85vh"
            )
          )
        )
      ),
      fluidRow(
        cerebroBox(
          title = tagList(
            boxTitle("Distribution along pseudotime"),
            cerebroInfoButton("trajectory_density_info")
          ),
          plotly::plotlyOutput("trajectory_density_plot")
        )
      ),
      fluidRow(
        cerebroBox(
          title = tagList(
            boxTitle("Number of cells by state"),
            cerebroInfoButton("trajectory_cells_by_state_info")
          ),
          DT::dataTableOutput("trajectory_number_of_cells_by_state_table")
        )
      ),
      fluidRow(
        cerebroBox(
          title = tagList(
            boxTitle("States by sample"),
            cerebroInfoButton("states_by_sample_info")
          ),
          tagList(
            uiOutput("states_by_sample_UI_buttons"),
            uiOutput("states_by_sample_UI_rest")
          )
        )
      ),
      fluidRow(
        cerebroBox(
          title = tagList(
            boxTitle("States by cluster"),
            cerebroInfoButton("states_by_cluster_info")
          ),
          tagList(
            uiOutput("states_by_cluster_UI_buttons"),
            uiOutput("states_by_cluster_UI_rest")
          )
        )
      ),
      fluidRow(
        cerebroBox(
          title = tagList(
            boxTitle("States by cell cycle (Seurat)"),
            cerebroInfoButton("states_by_cell_cycle_seurat_info")
          ),
          shiny::uiOutput("states_by_cell_cycle_seurat_UI")
        )
      ),
      fluidRow(
        cerebroBox(
          title = tagList(
            boxTitle("Number of transcripts by state"),
            cerebroInfoButton("states_nUMI_info")
          ),
          plotly::plotlyOutput("states_nUMI_plot")
        )
      ),
      fluidRow(
        cerebroBox(
          title = tagList(
            boxTitle("Number of expressed genes by state"),
            cerebroInfoButton("states_nGene_info")
          ),
          plotly::plotlyOutput("states_nGene_plot")
        )
      )
    )
  } else {
    cerebroBox(title = "Trajectory", textOutput("trajectory_missing"))
  }
})

# alternative text
output[["trajectory_missing"]] <- renderText({
  "No trajectories available to display."
})

##----------------------------------------------------------------------------##
## UI elements.
##----------------------------------------------------------------------------##
output[["trajectory_input"]] <- renderUI({
  tagList(
    selectInput(
      "trajectory_to_display",
      label = "Trajectory",
      choices = names(sample_data()$trajectory$monocle2)
    ),
    shinyWidgets::pickerInput(
      "trajectory_samples_to_display",
      label = "Samples to display",
      choices = sample_data()$sample_names,
      selected = sample_data()$sample_names,
      options = list("actions-box" = TRUE),
      multiple = TRUE
    ),
    shinyWidgets::pickerInput(
      "trajectory_clusters_to_display",
      label = "Clusters to display",
      choices = sample_data()$cluster_names,
      selected = sample_data()$cluster_names,
      options = list("actions-box" = TRUE),
      multiple = TRUE
    ),
    sliderInput(
      "trajectory_percentage_cells_to_show",
      label = "Show % of cells",
      min = scatter_plot_percentage_cells_to_show[["min"]],
      max = scatter_plot_percentage_cells_to_show[["max"]],
      step = scatter_plot_percentage_cells_to_show[["step"]],
      value = scatter_plot_percentage_cells_to_show[["default"]]
    ),
    selectInput(
      "trajectory_dot_color",
      label = "Color cells by",
      choices = c("state","pseudotime",names(sample_data()$cells)[! names(sample_data()$cells) %in% c("cell_barcode")])
    ),
    sliderInput(
      "trajectory_dot_size",
      label = "Dot size",
      min = scatter_plot_dot_size[["min"]],
      max = scatter_plot_dot_size[["max"]],
      step = scatter_plot_dot_size[["step"]],
      value = scatter_plot_dot_size[["default"]]
    ),
    sliderInput(
      "trajectory_dot_opacity",
      label = "Dot opacity",
      min = scatter_plot_dot_opacity[["min"]],
      max = scatter_plot_dot_opacity[["max"]],
      step = scatter_plot_dot_opacity[["step"]],
      value = scatter_plot_dot_opacity[["default"]]
    )
  )
})

##----------------------------------------------------------------------------##
## Projection.
##----------------------------------------------------------------------------##
output[["trajectory_projection"]] <- plotly::renderPlotly({
  # don't do anything before these inputs are selected
  req(
    input[["trajectory_to_display"]],
    input[["trajectory_samples_to_display"]],
    input[["trajectory_clusters_to_display"]],
    input[["trajectory_percentage_cells_to_show"]],
    input[["trajectory_dot_color"]],
    input[["trajectory_dot_size"]],
    input[["trajectory_dot_opacity"]]
  )

  trajectory_to_display <- input[["trajectory_to_display"]]
  samples_to_display <- input[["trajectory_samples_to_display"]]
  clusters_to_display <- input[["trajectory_clusters_to_display"]]
  cells_to_display <- which(
    (sample_data()$cells$sample %in% samples_to_display) &
    (sample_data()$cells$cluster %in% clusters_to_display)
  )

  # randomly remove cells
  if ( input[["trajectory_percentage_cells_to_show"]] < 100 ) {
    number_of_cells_to_plot <- ceiling(
      input[["trajectory_percentage_cells_to_show"]] / 100 * length(cells_to_display)
    )
    cells_to_display <- cells_to_display[ sample(1:length(cells_to_display), number_of_cells_to_plot) ]
  }

  # extract cells to plot
  to_plot <- cbind(
      sample_data()$trajectory$monocle2[[ trajectory_to_display ]][["meta"]][ cells_to_display , ],
      sample_data()$cells[ cells_to_display , ]
    ) %>%
    dplyr::filter(!is.na(pseudotime))
  to_plot <- to_plot[ sample(1:nrow(to_plot)) , ]

  color_variable <- input[["trajectory_dot_color"]]

  # convert edges of trajectory into list format to plot with plotly
  trajectory_edges <- sample_data()$trajectory$monocle2[[trajectory_to_display]][["edges"]]
  trajectory_lines <- list()
  for (i in 1:nrow(trajectory_edges) ) {
    line = list(
      type = "line",
      line = list(color = "black"),
      xref = "x",
      yref = "y",
      x0 = trajectory_edges$source_dim_1[i],
      y0 = trajectory_edges$source_dim_2[i],
      x1 = trajectory_edges$target_dim_1[i],
      y1 = trajectory_edges$target_dim_2[i]
    )
    trajectory_lines <- c(trajectory_lines, list(line))
  }

  if ( is.factor(to_plot[[ color_variable ]]) || is.character(to_plot[[ color_variable ]]) ) {
    if ( color_variable == "sample" ) {
      colors_this_plot <- reactive_colors()$sampless
    } else if ( color_variable == "cluster" ) {
      colors_this_plot <- reactive_colors()$clusters
    } else if ( color_variable %in% c("cell_cycle_seurat","cell_cycle_cyclone") ) {
      colors_this_plot <- cell_cycle_colorset
    } else if ( is.factor(to_plot[[ color_variable ]]) ) {
      colors_this_plot <- setNames(
        default_colorset[1:length(levels(to_plot[[ color_variable ]]))],
        levels(to_plot[[ color_variable ]])
      )
    } else {
      colors_this_plot <- default_colorset
    }
    plot <- plotly::plot_ly(
      to_plot,
      x = ~DR_1,
      y = ~DR_2,
      color = ~to_plot[[ color_variable ]],
      colors = colors_this_plot,
      type = "scatter",
      mode = "markers",
      marker = list(
        opacity = input[["trajectory_dot_opacity"]],
        line = list(
          color = "rgb(196,196,196)",
          width = 1
        ),
        size = input[["trajectory_dot_size"]]
      ),
      hoverinfo = "text",
      text = ~paste(
        "<b>Cell</b>: ", to_plot[ , "cell_barcode" ], "<br>",
        "<b>Sample</b>: ", to_plot[ , "sample" ], "<br>",
        "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br>",
        "<b>Transcripts</b>: ", formatC(to_plot[ , "nUMI" ], format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>Expressed genes</b>: ", formatC(to_plot[ , "nGene" ], format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>State</b>: ", to_plot[ , "state" ], "<br>",
        "<b>Pseudotime</b>: ", round(to_plot[ , "pseudotime" ], 3)
      )
    ) %>%
    plotly::layout(
      shapes = trajectory_lines,
      xaxis = list(
        mirror = TRUE,
        showline = TRUE,
        zeroline = FALSE,
        range = range(to_plot$DR_1) * 1.1
      ),
      yaxis = list(
        mirror = TRUE,
        showline = TRUE,
        zeroline = FALSE,
        range = range(to_plot$DR_2) * 1.1
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
      data = to_plot,
      x = ~DR_1,
      y = ~DR_2,
      type = "scatter",
      mode = "markers",
      marker = list(
        colorbar = list(
          title = colnames(to_plot)[which(colnames(to_plot) == color_variable)]
        ),
        color = ~to_plot[[ color_variable ]],
        opacity = input[["trajectory_dot_opacity"]],
        colorscale = "YlGnBu",
        reversescale = TRUE,
        line = list(
          color = "rgb(196,196,196)",
          width = 1
        ),
        size = input[["trajectory_dot_size"]]
      ),
      hoverinfo = "text",
      text = ~paste(
        "<b>Cell</b>: ", to_plot[ , "cell_barcode" ], "<br>",
        "<b>Sample</b>: ", to_plot[ , "sample" ], "<br>",
        "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br>",
        "<b>Transcripts</b>: ", formatC(to_plot[ , "nUMI" ], format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>Expressed genes</b>: ", formatC(to_plot[ , "nGene" ], format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>State</b>: ", to_plot[ , "state" ], "<br>",
        "<b>Pseudotime</b>: ", round(to_plot[ , "pseudotime" ], 3)
      )
    ) %>%
    plotly::layout(
      shapes = trajectory_lines,
      xaxis = list(
        title = colnames(to_plot)[1],
        mirror = TRUE,
        showline = TRUE,
        zeroline = FALSE,
        range = range(to_plot$DR_1) * 1.1
      ),
      yaxis = list(
        title = colnames(to_plot)[2],
        mirror = TRUE,
        showline = TRUE,
        zeroline = FALSE,
        range = range(to_plot$DR_2) * 1.1
      ),
      hoverlabel = list(font = list(size = 11))
    )
    if ( preferences$use_webgl == TRUE ) {
      plotly::toWebGL(plot)
    } else {
      plot
    }
  }
})

##----------------------------------------------------------------------------##
## Info button.
##----------------------------------------------------------------------------##
observeEvent(input[["trajectory_projection_info"]], {
  showModal(
    modalDialog(
      trajectory_projection_info[["text"]],
      title = trajectory_projection_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Export projection.
##----------------------------------------------------------------------------##
observeEvent(input[["trajectory_projection_export"]], {
  req(
    input[["trajectory_to_display"]],
    input[["trajectory_samples_to_display"]],
    input[["trajectory_clusters_to_display"]],
    input[["trajectory_percentage_cells_to_show"]],
    input[["trajectory_dot_color"]],
    input[["trajectory_dot_size"]],
    input[["trajectory_dot_opacity"]]
  )
  library("ggplot2")
  if ( exists("plot_export_path") ) {
    trajectory_to_display <- input[["trajectory_to_display"]]
    samples_to_display <- input[["trajectory_samples_to_display"]]
    clusters_to_display <- input[["trajectory_clusters_to_display"]]
    cells_to_display <- which(
        (sample_data()$cells$sample %in% samples_to_display) &
        (sample_data()$cells$cluster %in% clusters_to_display)
      )
    to_plot <- cbind(
        sample_data()$trajectory$monocle2[[ trajectory_to_display ]][[ "meta" ]][ cells_to_display , ],
        sample_data()$cells[ cells_to_display , ]
      ) %>%
      dplyr::filter(!is.na(pseudotime))
    to_plot <- to_plot[ sample(1:nrow(to_plot)) , ]

    color_variable <- input[["trajectory_dot_color"]]

    if ( is.factor(to_plot[[ color_variable ]]) || is.character(to_plot[[ color_variable ]]) ) {
      if ( color_variable == "sample" ) {
        colors_this_plot <- reactive_colors()$samples
      } else if ( color_variable == "cluster" ) {
        colors_this_plot <- reactive_colors()$clusters
      } else if ( color_variable %in% c("cell_cycle_seurat","cell_cycle_cyclone") ) {
        colors_this_plot <- cell_cycle_colorset
      } else if ( is.factor(to_plot[ , color_variable ]) ) {
        colors_this_plot <- setNames(
          default_colorset[1:length(levels(to_plot[ , color_variable ]))],
          levels(to_plot[ , color_variable ])
        )
      } else {
        colors_this_plot <- default_colorset
      }
      p <- ggplot() +
        geom_point(
          data = to_plot,
          aes_string(x = colnames(to_plot)[1], y = colnames(to_plot)[2], fill = color_variable),
          shape = 21,
          size = input[["trajectory_dot_size"]]/3,
          stroke = 0.2,
          color = "#c4c4c4",
          alpha = input[["trajectory_dot_opacity"]]
        ) +
        geom_segment(
          data = sample_data()$trajectory$monocle2[[ trajectory_to_display ]]$edges,
          aes(source_dim_1, source_dim_2, xend = target_dim_1, yend = target_dim_2),
          size = 0.75, linetype = "solid", na.rm = TRUE
        ) +
        scale_fill_manual(values = colors_this_plot) +
        theme_bw()
    } else {
      p <- ggplot() +
        geom_point(
          data = to_plot,
          aes_string(x = colnames(to_plot)[1], y = colnames(to_plot)[2], fill = color_variable),
          shape = 21,
          size = input[["trajectory_dot_size"]]/3,
          stroke = 0.2,
          color = "#c4c4c4",
          alpha = input[["trajectory_dot_opacity"]]
        ) +
        geom_segment(
          data = sample_data()$trajectory$monocle2[[ trajectory_to_display ]]$edges,
          aes(source_dim_1, source_dim_2, xend = target_dim_1, yend = target_dim_2),
          size = 0.75, linetype = "solid", na.rm = TRUE
        ) +
        scale_fill_distiller(
          palette = "YlGnBu",
          direction = 1,
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
        ) +
        theme_bw()
    }

    out_filename <- paste0(
        plot_export_path, "Cerebro_",
        gsub(
          sample_data()$experiment$experiment_name,
          pattern = " ", replacement = "_"
        ),
        "_trajectory_", trajectory_to_display, "_by_",
        gsub(
          color_variable,
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
## Distribution along pseudotime.
##----------------------------------------------------------------------------##

output[["trajectory_density_plot"]] <- plotly::renderPlotly({
  # don't do anything before these inputs are selected
  req(
    input[["trajectory_to_display"]],
    input[["trajectory_samples_to_display"]],
    input[["trajectory_clusters_to_display"]],
    input[["trajectory_dot_color"]]
  )

  trajectory_to_display <- input[["trajectory_to_display"]]
  samples_to_display <- input[["trajectory_samples_to_display"]]
  clusters_to_display <- input[["trajectory_clusters_to_display"]]
  cells_to_display <- which(
      (sample_data()$cells$sample %in% samples_to_display) &
      (sample_data()$cells$cluster %in% clusters_to_display)
    )

  # extract cells to plot
  to_plot <- cbind(
      sample_data()$trajectory$monocle2[[ trajectory_to_display ]][["meta"]][ cells_to_display , ],
      sample_data()$cells[ cells_to_display , ]
    ) %>%
    dplyr::filter(!is.na(pseudotime))
  to_plot <- to_plot[ sample(1:nrow(to_plot)) , ]

  color_variable <- input[["trajectory_dot_color"]]

  if ( is.factor(to_plot[[ color_variable ]]) || is.character(to_plot[[ color_variable ]]) ) {
    if ( color_variable == "sample" ) {
      colors_this_plot <- reactive_colors()$samples
    } else if ( color_variable == "cluster" ) {
      colors_this_plot <- reactive_colors()$clusters
    } else if ( color_variable %in% c("cell_cycle_seurat","cell_cycle_cyclone") ) {
      colors_this_plot <- cell_cycle_colorset
    } else if ( is.factor(to_plot[[ color_variable ]]) ) {
      colors_this_plot <- setNames(
        default_colorset[1:length(levels(to_plot[[ color_variable ]]))],
        levels(to_plot[[ color_variable ]])
      )
    } else {
      colors_this_plot <- default_colorset
    }
    p <- ggplot(to_plot, aes_string(x = "pseudotime", fill = color_variable)) +
      geom_density(alpha = 0.4, color = "black") +
      theme_bw() +
      labs(x = "Pseudotime", y = "Density") +
      scale_fill_manual(values = colors_this_plot) +
      guides(fill = guide_legend(override.aes = list(alpha = 1)))
    plotly::ggplotly(p, tooltip = "text") %>%
    plotly::style(
      hoveron = "fill"
    )
  } else {
    colors_this_plot <- setNames(
      default_colorset[1:length(levels(sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]]$state))],
      levels(sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]]$state)
    )
    plot <- plotly::plot_ly(
      data = to_plot,
      x = ~pseudotime,
      y = ~to_plot[[ color_variable ]],
      type = "scatter",
      mode = "markers",
      color = ~state,
      colors = colors_this_plot,
      marker = list(
        opacity = input[["trajectory_dot_opacity"]],
        line = list(
          color = "rgb(196,196,196)",
          width = 1
        ),
        size = input[["trajectory_dot_size"]]
      ),
      hoverinfo = "text",
      text = ~paste(
        "<b>Cell</b>: ", to_plot[ , "cell_barcode" ], "<br>",
        "<b>Sample</b>: ", to_plot[ , "sample" ], "<br>",
        "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br>",
        "<b>Transcripts</b>: ", formatC(to_plot[ , "nUMI" ], format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>Expressed genes</b>: ", formatC(to_plot[ , "nGene" ], format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>State</b>: ", to_plot[ , "state" ], "<br>",
        "<b>Pseudotime</b>: ", round(to_plot[ , "pseudotime" ], 3)
      )
    ) %>%
    plotly::layout(
      xaxis = list(
        title = "Pseudotime",
        mirror = TRUE,
        showline = TRUE,
        zeroline = FALSE
      ),
      yaxis = list(
        title = color_variable,
        mirror = TRUE,
        showline = TRUE,
        zeroline = FALSE
      ),
      hoverlabel = list(font = list(size = 11))
    )
    if ( preferences$use_webgl == TRUE ) {
      plotly::toWebGL(plot)
    } else {
      plot
    }
  }
})

# info button
observeEvent(input[["trajectory_density_info"]], {
  showModal(
    modalDialog(
      trajectory_density_info[["text"]],
      title = trajectory_density_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Table of cells by state.
##----------------------------------------------------------------------------##

# table
output[["trajectory_number_of_cells_by_state_table"]] <- DT::renderDataTable(server = FALSE, {
  req(,
    input[["trajectory_to_display"]],
    input[["trajectory_dot_color"]]
  )
  if ( is.numeric(sample_data()$cells[[ input[["trajectory_dot_color"]] ]]) || input[["trajectory_dot_color"]] == "state" || input[["trajectory_dot_color"]] == "pseudotime" ) {
    table <- sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]] %>%
      dplyr::filter(!is.na(pseudotime)) %>%
      dplyr::group_by(state) %>%
      dplyr::summarize(total_cell_count = n())
  } else {
    table <- cbind(
        sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]],
        sample_data()$cells[[ input[["trajectory_dot_color"]] ]]
      ) %>%
      dplyr::filter(!is.na(pseudotime)) %>%
      dplyr::rename(selected_group = 5) %>%
      dplyr::group_by(state, selected_group) %>%
      dplyr::summarize(count = n()) %>%
      tidyr::spread(selected_group, count, fill = 0) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(total_cell_count = rowSums(.[,2:ncol(.)])) %>%
      dplyr::select(state,total_cell_count,everything())
  }
  DT::datatable(
    table,
    filter = "none",
    selection = "none",
    escape = FALSE,
    autoHideNavigation = TRUE,
    rownames = FALSE,
    extensions = c("Buttons"),
    class = "cell-border stripe",
    options = list(
      dom = "Brti",
      pageLength = 100,
      buttons = list(
        "colvis",
        list(
          extend = "collection",
          text = "Download",
          buttons = list(
            list(
              extend = "csv",
              filename = "cells_by_state",
              title = "Cells by state"
            ),
            list(
              extend = "excel",
              filename = "cells_by_state",
              title = "Cells by state"
            )
          )
        )
      )
    )
  )
})

# info button
observeEvent(input[["trajectory_number_of_cells_by_state_info"]], {
  showModal(
    modalDialog(
      trajectory_number_of_cells_by_state_info[["text"]],
      title = trajectory_number_of_cells_by_state_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## States by sample.
##----------------------------------------------------------------------------##

# UI element: buttons
output[["states_by_sample_UI_buttons"]] <- renderUI({
  tagList(
    shinyWidgets::materialSwitch(
      inputId = "states_by_sample_select_metric_for_bar_plot",
      label = "Show composition in percent [%]:",
      status = "primary",
      inline = TRUE
    ),
    shinyWidgets::materialSwitch(
      inputId = "states_by_sample_show_table",
      label = "Show table:",
      status = "primary",
      inline = TRUE
    )
  )
})

# UI element: rest
output[["states_by_sample_UI_rest"]] <- renderUI({
  tagList(
    plotly::plotlyOutput("states_by_sample_plot"),
    {
      if ( !is.null(input[["states_by_sample_show_table"]]) && input[["states_by_sample_show_table"]] == TRUE ) {
        DT::dataTableOutput("states_by_sample_table")
      }
    }
  )
})

# bar plot
output[["states_by_sample_plot"]] <- plotly::renderPlotly({
  req(input[["trajectory_to_display"]])
  # merge meta data with trajectory info
  cell_count_by_state_total <- cbind(
      sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]],
      sample_data()$cells
    ) %>%
    dplyr::filter(!is.na(pseudotime)) %>%
    dplyr::group_by(state) %>%
    dplyr::summarize(total = n()) %>%
    dplyr::ungroup()
  # make plot
  temp_data <- cbind(
      sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]],
      sample_data()$cells
    ) %>%
    dplyr::filter(!is.na(pseudotime)) %>%
    dplyr::group_by(state, sample) %>%
    dplyr::summarize(count = n()) %>%
    tidyr::spread(sample, count, fill = 0) %>%
    dplyr::ungroup() %>%
    reshape2::melt(id.vars = "state") %>%
    dplyr::left_join(., cell_count_by_state_total, by = "state") %>%
    dplyr::rename(sample = variable, cells = value)
  if ( input[["states_by_sample_select_metric_for_bar_plot"]] != TRUE ) {
    temp_data %>%
    plotly::plot_ly(
      x = ~state,
      y = ~cells,
      type = "bar",
      color = ~sample,
      colors = reactive_colors()$samples,
      hoverinfo = "text",
      text = ~paste0("<b>Sample ", .$sample, ": </b>", formatC(.$cells, big.mark = ','))
    ) %>%
    plotly::layout(
      xaxis = list(
        title = "",
        mirror = TRUE,
        showline = TRUE
      ),
      yaxis = list(
        title = "Number of cells",
        hoverformat = ".2f",
        mirror = TRUE,
        zeroline = FALSE,
        showline = TRUE
      ),
      barmode = "stack",
      hovermode = "compare"
    )
  } else {
    temp_data %>%
    dplyr::mutate(pct = cells / total * 100) %>%
    plotly::plot_ly(
      x = ~state,
      y = ~pct,
      type = "bar",
      color = ~sample,
      colors = reactive_colors()$samples,
      hoverinfo = "text",
      text = ~paste0("<b>Sample ", .$sample, ": </b>", format(round(.$pct, 1), nsmall = 1), "%")
    ) %>%
    plotly::layout(
      xaxis = list(
        title = "",
        mirror = TRUE,
        showline = TRUE
      ),
      yaxis = list(
        title = "Percentage (%)",
        range = c(0,100),
        hoverformat = ".2f",
        mirror = TRUE,
        zeroline = FALSE,
        showline = TRUE
      ),
      barmode = "stack",
      hovermode = "compare"
    )
  }
})

# table
output[["states_by_sample_table"]] <- DT::renderDataTable({
  req(input[["trajectory_to_display"]])
  # generate table
  temp_table <- cbind(
      sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]],
      sample_data()$cells
    ) %>%
    dplyr::filter(!is.na(pseudotime)) %>%
    dplyr::group_by(state,sample) %>%
    dplyr::summarize(count = n()) %>%
    tidyr::spread(sample, count, fill = 0) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    dplyr::select(c(state, 'total_cell_count', dplyr::everything()))
  if ( input[["states_by_sample_select_metric_for_bar_plot"]] == TRUE ) {
    # normalize counts to 100% percent
    for ( i in 3:ncol(temp_table) ) {
      temp_table[,i] <- round(temp_table[,i] / temp_table$total_cell_count * 100, digits = 1)
    }
  }
  # process table and convert to DT
  temp_table %>%
  rename(
    State = state,
    "# of cells" = total_cell_count
  ) %>%
  DT::datatable(
    filter = "none",
    selection = "none",
    escape = FALSE,
    autoHideNavigation = TRUE,
    rownames = FALSE,
    class = "cell-border stripe",
    options = list(
      scrollX = TRUE,
      sDom = '<"top">lrt<"bottom">ip',
      lengthMenu = c(15, 30, 50, 100),
      pageLength = 15
    )
  )
})

# info button
observeEvent(input[["states_by_sample_info"]], {
  showModal(
    modalDialog(
      states_by_sample_info[["text"]],
      title = states_by_sample_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## States by cluster.
##----------------------------------------------------------------------------##

# UI element: buttons
output[["states_by_cluster_UI_buttons"]] <- renderUI({
  tagList(
    shinyWidgets::materialSwitch(
      inputId = "states_by_cluster_select_metric_for_bar_plot",
      label = "Show composition in percent [%]:",
      status = "primary",
      inline = TRUE
    ),
    shinyWidgets::materialSwitch(
      inputId = "states_by_cluster_show_table",
      label = "Show table:",
      status = "primary",
      inline = TRUE
    )
  )
})

# UI element: rest
output[["states_by_cluster_UI_rest"]] <- renderUI({
  tagList(
    plotly::plotlyOutput("states_by_cluster_plot"),
    {
      if ( !is.null(input[["states_by_cluster_show_table"]]) && input[["states_by_cluster_show_table"]] == TRUE ) {
        DT::dataTableOutput("states_by_cluster_table")
      }
    }
  )
})

# bar plot
output[["states_by_cluster_plot"]] <- plotly::renderPlotly({
  req(input[["trajectory_to_display"]])
  # merge meta data with trajectory info
  cell_count_by_state_total <- cbind(
      sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]],
      sample_data()$cells
    ) %>%
    dplyr::filter(!is.na(pseudotime)) %>%
    dplyr::group_by(state) %>%
    dplyr::summarize(total = n()) %>%
    dplyr::ungroup()
  # make plot
  temp_data <- cbind(
      sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]],
      sample_data()$cells
    ) %>%
    dplyr::filter(!is.na(pseudotime)) %>%
    dplyr::group_by(state, cluster) %>%
    dplyr::summarize(count = n()) %>%
    tidyr::spread(cluster, count, fill = 0) %>%
    dplyr::ungroup() %>%
    reshape2::melt(id.vars = "state") %>%
    dplyr::left_join(., cell_count_by_state_total, by = "state") %>%
    dplyr::rename(cluster = variable, cells = value)
  if ( input[["states_by_cluster_select_metric_for_bar_plot"]] != TRUE ) {
    temp_data %>%
    plotly::plot_ly(
      x = ~state,
      y = ~cells,
      type = "bar",
      color = ~cluster,
      colors = reactive_colors()$clusters,
      hoverinfo = "text",
      text = ~paste0("<b>Cluster ", .$cluster, ": </b>", formatC(.$cells, big.mark = ','))
    ) %>%
    plotly::layout(
      xaxis = list(
        title = "",
        mirror = TRUE,
        showline = TRUE
      ),
      yaxis = list(
        title = "Number of cells",
        hoverformat = ".2f",
        mirror = TRUE,
        zeroline = FALSE,
        showline = TRUE
      ),
      barmode = "stack",
      hovermode = "compare"
    )
  } else {
    temp_data %>%
    dplyr::mutate(pct = cells / total * 100) %>%
    plotly::plot_ly(
      x = ~state,
      y = ~pct,
      type = "bar",
      color = ~cluster,
      colors = reactive_colors()$clusters,
      hoverinfo = "text",
      text = ~paste0("<b>Cluster ", .$cluster, ": </b>", format(round(.$pct, 1), nsmall = 1), "%")
    ) %>%
    plotly::layout(
      xaxis = list(
        title = "",
        mirror = TRUE,
        showline = TRUE
      ),
      yaxis = list(
        title = "Percentage (%)",
        range = c(0,100),
        hoverformat = ".2f",
        mirror = TRUE,
        zeroline = FALSE,
        showline = TRUE
      ),
      barmode = "stack",
      hovermode = "compare"
    )
  }
})

# table
output[["states_by_cluster_table"]] <- DT::renderDataTable({
  req(input[["trajectory_to_display"]])
  # generate table
  temp_table <- cbind(
      sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]],
      sample_data()$cells
    ) %>%
    dplyr::filter(!is.na(pseudotime)) %>%
    dplyr::group_by(state,cluster) %>%
    dplyr::summarize(count = n()) %>%
    tidyr::spread(cluster, count, fill = 0) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    dplyr::select(c(state, 'total_cell_count', dplyr::everything()))
  if ( input[["states_by_cluster_select_metric_for_bar_plot"]] == TRUE ) {
    # normalize counts to 100% percent
    for ( i in 3:ncol(temp_table) ) {
      temp_table[,i] <- round(temp_table[,i] / temp_table$total_cell_count * 100, digits = 1)
    }
  }
  # process table and convert to DT
  temp_table %>%
  rename(
    State = state,
    "# of cells" = total_cell_count
  ) %>%
  DT::datatable(
    filter = "none",
    selection = "none",
    escape = FALSE,
    autoHideNavigation = TRUE,
    rownames = FALSE,
    class = "cell-border stripe",
    options = list(
      scrollX = TRUE,
      sDom = '<"top">lrt<"bottom">ip',
      lengthMenu = c(15, 30, 50, 100),
      pageLength = 15
    )
  )
})

# info button
observeEvent(input[["states_by_cluster_info"]], {
  showModal(
    modalDialog(
      states_by_cluster_info[["text"]],
      title = states_by_cluster_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## States by cell cycle status (Seurat).
##----------------------------------------------------------------------------##

# UI element
output[["states_by_cell_cycle_seurat_UI"]] <- renderUI({
  if ( "cell_cycle_seurat" %in% colnames(sample_data()$cells) ) {
    tagList(
      uiOutput("states_by_cell_cycle_seurat_UI_buttons"),
      uiOutput("states_by_cell_cycle_seurat_UI_rest")
    )
  } else {
    textOutput("states_by_cell_cycle_seurat_text")
  }
})

# UI element: buttons
output[["states_by_cell_cycle_seurat_UI_buttons"]] <- renderUI({
  tagList(
    shinyWidgets::materialSwitch(
      inputId = "states_by_cell_cycle_seurat_select_metric_for_bar_plot",
      label = "Show composition in percent [%]:",
      status = "primary",
      inline = TRUE
    ),
    shinyWidgets::materialSwitch(
      inputId = "states_by_cell_cycle_seurat_show_table",
      label = "Show table:",
      status = "primary",
      inline = TRUE
    )
  )
})

# UI element: rest
output[["states_by_cell_cycle_seurat_UI_rest"]] <- renderUI({
  tagList(
    plotly::plotlyOutput("states_by_cell_cycle_seurat_plot"),
    {
      if ( !is.null(input[["states_by_cell_cycle_seurat_show_table"]]) && input[["states_by_cell_cycle_seurat_show_table"]] == TRUE ) {
        DT::dataTableOutput("states_by_cell_cycle_seurat_table")
      }
    }
  )
})

output[["states_by_cell_cycle_seurat_plot"]] <- plotly::renderPlotly({
  req(input[["trajectory_to_display"]])
  # merge meta data with trajectory info
  cell_count_by_state_total <- cbind(
      sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]],
      sample_data()$cells
    ) %>%
    dplyr::filter(!is.na(pseudotime)) %>%
    dplyr::group_by(state) %>%
    dplyr::summarize(total = n()) %>%
    dplyr::ungroup()
  # make plot
  temp_data <- cbind(
      sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]],
      sample_data()$cells
    ) %>%
    dplyr::filter(!is.na(pseudotime)) %>%
    dplyr::group_by(state, cell_cycle_seurat) %>%
    dplyr::summarize(count = n()) %>%
    tidyr::spread(cell_cycle_seurat, count, fill = 0) %>%
    dplyr::ungroup() %>%
    reshape2::melt(id.vars = "state") %>%
    dplyr::mutate(
      variable = factor(variable, levels = c("G1", "S", "G2M")),
    ) %>%
    dplyr::left_join(., cell_count_by_state_total, by = "state") %>%
    dplyr::rename(cell_cycle_seurat = variable, cells = value)
  if ( input[["states_by_cell_cycle_seurat_select_metric_for_bar_plot"]] != TRUE ) {
    temp_data %>%
    plotly::plot_ly(
      x = ~state,
      y = ~cells,
      type = "bar",
      color = ~cell_cycle_seurat,
      colors = cell_cycle_colorset,
      hoverinfo = "text",
      text = ~paste0("<b>", .$cell_cycle_seurat, ": </b>", formatC(.$cells, big.mark = ','))
    ) %>%
    plotly::layout(
      xaxis = list(
        title = "",
        mirror = TRUE,
        showline = TRUE
      ),
      yaxis = list(
        title = "Number of cells",
        hoverformat = ".2f",
        mirror = TRUE,
        zeroline = FALSE,
        showline = TRUE
      ),
      barmode = "stack",
      hovermode = "compare"
    )
  } else {
    temp_data %>%
    dplyr::mutate(pct = cells / total * 100) %>%
    plotly::plot_ly(
      x = ~state,
      y = ~pct,
      type = "bar",
      color = ~cell_cycle_seurat,
      colors = cell_cycle_colorset,
      hoverinfo = "text",
      text = ~paste0("<b>", .$cell_cycle_seurat, ": </b>", format(round(.$pct, 1), nsmall = 1), "%")
    ) %>%
    plotly::layout(
      xaxis = list(
        title = "",
        mirror = TRUE,
        showline = TRUE
      ),
      yaxis = list(
        title = "Percentage (%)",
        range = c(0,100),
        hoverformat = ".2f",
        mirror = TRUE,
        zeroline = FALSE,
        showline = TRUE
      ),
      barmode = "stack",
      hovermode = "compare"
    )
  }
})

# table
output[["states_by_cell_cycle_seurat_table"]] <- DT::renderDataTable({
  req(input[["trajectory_to_display"]])
  # generate table
  temp_table <- cbind(
      sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]],
      sample_data()$cells
    ) %>%
    dplyr::filter(!is.na(pseudotime)) %>%
    dplyr::group_by(state,cell_cycle_seurat) %>%
    dplyr::summarize(count = n()) %>%
    tidyr::spread(cell_cycle_seurat, count, fill = 0) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    dplyr::select(c(state, 'total_cell_count', 'G1', 'S', 'G2M', dplyr::everything()))
  if ( input[["states_by_cell_cycle_seurat_select_metric_for_bar_plot"]] == TRUE ) {
    # normalize counts to 100% percent
    for ( i in 3:ncol(temp_table) ) {
      temp_table[,i] <- round(temp_table[,i] / temp_table$total_cell_count * 100, digits = 1)
    }
  }
  # process table and convert to DT
  temp_table %>%
  rename(
    State = state,
    "# of cells" = total_cell_count
  ) %>%
  DT::datatable(
    filter = "none",
    selection = "none",
    escape = FALSE,
    autoHideNavigation = TRUE,
    rownames = FALSE,
    class = "cell-border stripe",
    options = list(
      scrollX = TRUE,
      sDom = '<"top">lrt<"bottom">ip',
      lengthMenu = c(15, 30, 50, 100),
      pageLength = 15
    )
  )
})

# alternative text
output[["states_by_cell_cycle_seurat_text"]] <- renderText({
  "Data not available."
})

# info button
observeEvent(input[["states_by_cell_cycle_seurat_info"]], {
  showModal(
    modalDialog(
      states_by_cell_cycle_seurat_info[["text"]],
      title = states_by_cell_cycle_seurat_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## nUMI by state.
##----------------------------------------------------------------------------##

# violin plot
output[["states_nUMI_plot"]] <- plotly::renderPlotly({
  req(input[["trajectory_to_display"]])
  colors_this_plot <- setNames(
    default_colorset[1:length(levels(sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]]$state))],
    levels(sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]]$state)
  )
  cbind(
    sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]],
    sample_data()$cells
  ) %>%
  dplyr::filter(!is.na(pseudotime)) %>%
  plotly::plot_ly(
    x = ~state,
    y = ~nUMI,
    type = "violin",
    box = list(
      visible = TRUE
    ),
    meanline = list(
      visible = TRUE
    ),
    color = ~state,
    colors = colors_this_plot,
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
      title = "Number of transcripts",
      hoverformat = ".0f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# info button
observeEvent(input[["states_nUMI_info"]], {
  showModal(
    modalDialog(
      states_nUMI_info[["text"]],
      title = states_nUMI_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## nGene by state.
##----------------------------------------------------------------------------##

# violin plot
output[["states_nGene_plot"]] <- plotly::renderPlotly({
  req(input[["trajectory_to_display"]])
  colors_this_plot <- setNames(
    default_colorset[1:length(levels(sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]]$state))],
    levels(sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]]$state)
  )
  cbind(
    sample_data()$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]],
    sample_data()$cells
  ) %>%
  dplyr::filter(!is.na(pseudotime)) %>%
  plotly::plot_ly(
    x = ~state,
    y = ~nGene,
    type = "violin",
    box = list(
      visible = TRUE
    ),
    meanline = list(
      visible = TRUE
    ),
    color = ~state,
    colors = colors_this_plot,
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
      title = "Number of expressed genes",
      hoverformat = ".0f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# info button
observeEvent(input[["states_nGene_info"]], {
  showModal(
    modalDialog(
      states_nGene_info[["text"]],
      title = states_nGene_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})
