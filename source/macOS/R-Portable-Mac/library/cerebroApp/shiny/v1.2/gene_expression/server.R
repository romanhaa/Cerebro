##----------------------------------------------------------------------------##
## Tab: Gene expression
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## UI elements for projection.
##----------------------------------------------------------------------------##
output[["expression_UI"]] <- renderUI({
  tagList(
    selectizeInput(
      'expression_genes_input',
      label = 'Gene(s)',
      choices = rownames(sample_data()$expression),
      options = list(create = TRUE), multiple = TRUE
    ),
    selectInput(
      "expression_projection_to_display",
      label = "Projection",
      choices = names(sample_data()$projections)
    ),
    shinyWidgets::pickerInput(
      "expression_samples_to_display",
      label = "Samples to display",
      choices = sample_data()$sample_names,
      selected = sample_data()$sample_names,
      options = list("actions-box" = TRUE),
      multiple = TRUE
    ),
    shinyWidgets::pickerInput(
      "expression_clusters_to_display",
      label = "Clusters to display",
      choices = sample_data()$cluster_names,
      selected = sample_data()$cluster_names,
      options = list("actions-box" = TRUE),
      multiple = TRUE
    ),
    sliderInput(
      "expression_percentage_cells_to_show",
      label = "Show % of cells",
      min = scatter_plot_percentage_cells_to_show[["min"]],
      max = scatter_plot_percentage_cells_to_show[["max"]],
      step = scatter_plot_percentage_cells_to_show[["step"]],
      value = scatter_plot_percentage_cells_to_show[["default"]]
    ),
    selectInput(
      "expression_projection_plotting_order",
      label = "Plotting order",
      choices = c("Random", "Highest expression on top"),
      selected = "Random"
    ),
    sliderInput(
      "expression_projection_dot_size",
      label = "Point size",
      min = scatter_plot_dot_size[["min"]],
      max = scatter_plot_dot_size[["max"]],
      step = scatter_plot_dot_size[["step"]],
      value = scatter_plot_dot_size[["default"]]
    ),
    sliderInput(
      "expression_projection_dot_opacity",
      label = "Point opacity",
      min = scatter_plot_dot_opacity[["min"]],
      max = scatter_plot_dot_opacity[["max"]],
      step = scatter_plot_dot_opacity[["step"]],
      value = scatter_plot_dot_opacity[["default"]]
    ),
    selectInput(
      "expression_projection_color_scale",
      label = "Color scale",
      choices = c("YlGnBu", "YlOrRd","Blues","Greens","Reds","RdBu","viridis"),
      selected = "YlGnBu"
    )
  )
})

##----------------------------------------------------------------------------##
## UI element for color scale range in projection.
##----------------------------------------------------------------------------##
output[["expression_color_scale_range"]] <- renderUI({
  range <- range(gene_expression_plot_data()$level)
  if ( range[1] == 0 & range[2] == 0 ) {
    range[2] = 1
  } else {
    range[1] <- range[1] %>% round(digits = 2)
    range[2] <- range[2] %>% round(digits = 2)
  }
  tagList(
    sliderInput(
      "expression_projection_color_scale_range",
      label = "Range of color scale",
      min = range[1],
      max = range[2],
      value = c(range[1], range[2])
    )
  )
})

##----------------------------------------------------------------------------##
## UI element for X and Y scales in projection.
##----------------------------------------------------------------------------##
output[["expression_scales"]] <- renderUI({
  req(input[["expression_projection_to_display"]])
  projection_to_display <- input[["expression_projection_to_display"]]
  range_x_min <- sample_data()$projections[[ projection_to_display ]][,1] %>% min() %>% "*"(ifelse(.<0, 1.1, 0.9)) %>% round()
  range_x_max <- sample_data()$projections[[ projection_to_display ]][,1] %>% max() %>% "*"(ifelse(.<0, 0.9, 1.1)) %>% round()
  range_y_min <- sample_data()$projections[[ projection_to_display ]][,2] %>% min() %>% "*"(ifelse(.<0, 1.1, 0.9)) %>% round()
  range_y_max <- sample_data()$projections[[ projection_to_display ]][,2] %>% max() %>% "*"(ifelse(.<0, 0.9, 1.1)) %>% round()
  tagList(
    sliderInput(
      "expression_projection_scale_x_manual_range",
      label = "Range of X axis",
      min = range_x_min,
      max = range_x_max,
      value = c(range_x_min, range_x_max)
    ),
    sliderInput(
      "expression_projection_scale_y_manual_range",
      label = "Range of Y axis",
      min = range_y_min,
      max = range_y_max,
      value = c(range_y_min, range_y_max)
    )
  )
})

##----------------------------------------------------------------------------##
## Reactive data: Genes from user.
##----------------------------------------------------------------------------##

# cannot use req() because it delays initialization and plot is updated only with button press so plot doesn't initialize at all
genesToPlot <- reactive({
  genesToPlot <- list()
  if ( is.null(input[["expression_genes_input"]]) ) {
    genesToPlot[["genes_to_display"]] <- character()
  } else {
    genesToPlot[["genes_to_display"]] <- input[["expression_genes_input"]] %>%
      strsplit(",| |;|\n") %>%
      unlist() %>%
      gsub(pattern = " ", replacement = "", fixed = TRUE) %>%
      unique() %>%
      .[. != ""]
  }
  genesToPlot[["genes_to_display_here"]] <- rownames(sample_data()$expression)[ match(tolower(genesToPlot[["genes_to_display"]]), tolower(rownames(sample_data()$expression))) ]
  genesToPlot[["genes_to_display_present"]] <- na.omit(genesToPlot[["genes_to_display_here"]])
  genesToPlot[["genes_to_display_missing"]] <- genesToPlot[["genes_to_display"]][ which(is.na(genesToPlot[["genes_to_display_here"]])) ]
  return(genesToPlot)
})

# select genes to be displayed
output[["expression_genes_displayed"]] <- renderText({
  paste0(
    "<b>Showing expression for ",
    length(genesToPlot()[["genes_to_display_present"]]), " gene(s):</b><br>",
    paste0(genesToPlot()[["genes_to_display_present"]], collapse = ", "),
    "<br><b>",
    length(genesToPlot()[["genes_to_display_missing"]]),
    " gene(s) are not in data set: </b><br>",
    paste0(genesToPlot()[["genes_to_display_missing"]], collapse = ", ")
  )
})

# data to plot
gene_expression_plot_data <- reactive({
  req(
    input[["expression_projection_to_display"]],
    input[["expression_samples_to_display"]],
    input[["expression_clusters_to_display"]],
    input[["expression_percentage_cells_to_show"]],
    input[["expression_projection_plotting_order"]]
  )
  projection_to_display <- input[["expression_projection_to_display"]]
  samples_to_display <- input[["expression_samples_to_display"]]
  clusters_to_display <- input[["expression_clusters_to_display"]]
  percentage_cells_show <- input[["expression_percentage_cells_to_show"]]
  plot_order <- input[["expression_projection_plotting_order"]]
  # check which cells to display
  cells_to_display <- which(
      (sample_data()$cells$sample %in% samples_to_display) &
      (sample_data()$cells$cluster %in% clusters_to_display)
    )
  # randomly remove cells
  if ( percentage_cells_show < 100 ) {
    number_of_cells_to_plot <- ceiling(
      percentage_cells_show / 100 * length(cells_to_display)
    )
    cells_to_display <- cells_to_display[ sample(1:length(cells_to_display), number_of_cells_to_plot) ]
  }
  plot <- cbind(
      sample_data()$projections[[ projection_to_display ]][ cells_to_display , ],
      sample_data()$cells[ cells_to_display , ]
    )
  if ( length(genesToPlot()$genes_to_display_present) == 0 ) {
    plot$level <- 0
  } else if ( length(genesToPlot()$genes_to_display_present) == 1 ) {
    plot$level <- genesToPlot()$genes_to_display_present %>%
      sample_data()$expression[ . , cells_to_display ]
  } else {
    plot$level <- genesToPlot()$genes_to_display_present %>%
      sample_data()$expression[ . , cells_to_display ] %>%
      Matrix::colMeans()
  }
  if ( plot_order == "Random" ) {
    plot <- sample(1:nrow(plot), nrow(plot)) %>%
      plot[ . , ]
  } else if ( plot_order == "Highest expression on top" ) {
    plot <- plot[ order(plot$level, decreasing = FALSE) , ]
  }
  return(plot)
})

##----------------------------------------------------------------------------##
## Projection.
##----------------------------------------------------------------------------##

output[["expression_projection"]] <- plotly::renderPlotly({
  req(
    input[["expression_projection_to_display"]],
    input[["expression_projection_dot_size"]],
    input[["expression_projection_dot_opacity"]],
    input[["expression_projection_color_scale"]],
    input[["expression_projection_color_scale_range"]],
    input[["expression_projection_scale_x_manual_range"]],
    input[["expression_projection_scale_y_manual_range"]]
  )
  if ( input[["expression_projection_color_scale"]] == 'viridis' ) {
    color_scale <- 'Viridis'
  } else {
    color_scale <- input[["expression_projection_color_scale"]]
  }
  if ( ncol(sample_data()$projections[[ input[["expression_projection_to_display"]] ]]) == 3 ) {
    plotly::plot_ly(
      gene_expression_plot_data(),
      x = gene_expression_plot_data()[,1],
      y = gene_expression_plot_data()[,2],
      z = gene_expression_plot_data()[,3],
      type = "scatter3d",
      mode = "markers",
      marker = list(
        colorbar = list(
          title = "Expression"
        ),
        color = ~level,
        opacity = input[["expression_projection_dot_opacity"]],
        colorscale = color_scale,
        cauto = FALSE,
        cmin = input[["expression_projection_color_scale_range"]][1],
        cmax = input[["expression_projection_color_scale_range"]][2],
        reversescale = TRUE,
        line = list(
          color = "rgb(196,196,196)",
          width = 1
        ),
        size = input[["expression_projection_dot_size"]]
      ),
      hoverinfo = "text",
      text = ~paste(
        "<b>Cell</b>: ", gene_expression_plot_data()$cell_barcode, "<br>",
        "<b>Sample</b>: ", gene_expression_plot_data()$sample, "<br>",
        "<b>Cluster</b>: ", gene_expression_plot_data()$cluster, "<br>",
        "<b>Transcripts</b>: ", formatC(gene_expression_plot_data()$nUMI, format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>Expressed genes</b>: ", formatC(gene_expression_plot_data()$nGene, format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>Expression level</b>: ", formatC(gene_expression_plot_data()$level, format = "f", big.mark = ",", digits = 3)
      ),
      source = "expression_projection"
    ) %>%
    plotly::layout(
      scene = list(
        xaxis = list(
          title = colnames(gene_expression_plot_data())[1],
          mirror = TRUE,
          showline = TRUE,
          zeroline = FALSE
        ),
        yaxis = list(
          title = colnames(gene_expression_plot_data())[2],
          mirror = TRUE,
          showline = TRUE,
          zeroline = FALSE
        ),
        zaxis = list(
          title = colnames(gene_expression_plot_data())[3],
          mirror = TRUE,
          showline = TRUE,
          zeroline = FALSE
        )
      ),
      hoverlabel = list(
        font = list(
          size = 11,
          color = "black"
        ),
        bgcolor = "lightgrey"
      )
    )
  } else {
    plot <- plotly::plot_ly(
      gene_expression_plot_data(),
      x = gene_expression_plot_data()[,1],
      y = gene_expression_plot_data()[,2],
      type = "scatter",
      mode = "markers",
      marker = list(
        colorbar = list(
          title = "Expression"
        ),
        color = ~level,
        opacity = input[["expression_projection_dot_opacity"]],
        colorscale = color_scale,
        cauto = FALSE,
        cmin = input[["expression_projection_color_scale_range"]][1],
        cmax = input[["expression_projection_color_scale_range"]][2],
        reversescale = TRUE,
        line = list(
          color = "rgb(196,196,196)",
          width = 1
        ),
        size = input[["expression_projection_dot_size"]]
      ),
      hoverinfo = "text",
      text = ~paste(
        "<b>Cell</b>: ", gene_expression_plot_data()$cell_barcode, "<br>",
        "<b>Sample</b>: ", gene_expression_plot_data()$sample, "<br>",
        "<b>Cluster</b>: ", gene_expression_plot_data()$cluster, "<br>",
        "<b>Transcripts</b>: ", formatC(gene_expression_plot_data()$nUMI, format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>Expressed genes</b>: ", formatC(gene_expression_plot_data()$nGene, format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>Expression level</b>: ", formatC(gene_expression_plot_data()$level, format = "f", big.mark = ",", digits = 3)
      ),
      source = "expression_projection"
    ) %>%
    plotly::layout(
      xaxis = list(
        title = colnames(gene_expression_plot_data())[1],
        mirror = TRUE,
        showline = TRUE,
        zeroline = FALSE,
        range = c(
          input[["expression_projection_scale_x_manual_range"]][1],
          input[["expression_projection_scale_x_manual_range"]][2]
        )
      ),
      yaxis = list(
        title = colnames(gene_expression_plot_data())[2],
        mirror = TRUE,
        showline = TRUE,
        zeroline = FALSE,
        range = c(
          input[["expression_projection_scale_y_manual_range"]][1],
          input[["expression_projection_scale_y_manual_range"]][2]
        )
      ),
      dragmode = "pan",
      hoverlabel = list(
        font = list(
          size = 11,
          color = "black"
        ),
        bgcolor = "lightgrey"
      )
    )
    if ( preferences$use_webgl == TRUE ) {
      plot %>% plotly::toWebGL()
    } else {
      plot
    }
  }
})

##----------------------------------------------------------------------------##
## Info box.
##----------------------------------------------------------------------------##
observeEvent(input[["expression_projection_info"]], {
  showModal(
    modalDialog(
      expression_projection_info$text,
      title = expression_projection_info$title,
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Export function.
##----------------------------------------------------------------------------##
observeEvent(input[["expression_projection_export"]], {
  req(
    input[["expression_projection_to_display"]],
    input[["expression_projection_plotting_order"]],
    input[["expression_projection_dot_size"]],
    input[["expression_projection_dot_opacity"]],
    input[["expression_projection_color_scale"]],
    input[["expression_projection_color_scale_range"]],
    input[["expression_projection_scale_x_manual_range"]],
    input[["expression_projection_scale_y_manual_range"]]
  )
  library("ggplot2")
  if ( exists("plot_export_path") ) {
    xlim <- c(
      input[["expression_projection_scale_x_manual_range"]][1],
      input[["expression_projection_scale_x_manual_range"]][2]
    )
    ylim <- c(
      input[["expression_projection_scale_y_manual_range"]][1],
      input[["expression_projection_scale_y_manual_range"]][2]
    )
    if ( length(genesToPlot()$genes_to_display_present) == 0 ) {
      out_filename <- paste0(
        plot_export_path, "Cerebro_",
        sample_data()$experiment$experiment_name, "_gene_expression_none"
      )
    } else if ( length(genesToPlot()$genes_to_display_present) == 1 ) {
      out_filename <- paste0(
        plot_export_path, "Cerebro_",
        sample_data()$experiment$experiment_name, "_gene_expression_",
        genesToPlot()$genes_to_display_present, "_",
        input[["expression_projection_to_display"]]
      )
    } else {
      out_filename <- paste0(
        plot_export_path, "Cerebro_",
        sample_data()$experiment$experiment_name, "_gene_expression_",
        genesToPlot()$genes_to_display_present[1],
        "_and_others_", input[["expression_projection_to_display"]]
      )
    }

    if ( input[["expression_projection_plotting_order"]] == "Random" ) {
      out_filename <- paste0(out_filename, "_random_order.pdf")
    } else if ( input[["expression_projection_plotting_order"]] == "Highest expression on top" ) {
      out_filename <- paste0(out_filename, "_highest_expression_on_top.pdf")
    }

    if ( ncol(sample_data()$projections[[ input[["expression_projection_to_display"]] ]]) == 3 ) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Sorry!",
        text = "It's currently not possible to create PDF plots from 3D dimensional reductions. Please use the PNG export button in the panel or a 2D dimensional reduction instead.",
        type = "error"
      )
    } else {
      p <- ggplot(
          gene_expression_plot_data(),
          aes_q(
            x = as.name(colnames(gene_expression_plot_data())[1]),
            y = as.name(colnames(gene_expression_plot_data())[2]),
            fill = as.name("level")
          )
        ) +
        geom_point(
          shape = 21,
          size = input[["expression_projection_dot_size"]]/3,
          stroke = 0.2,
          color = "#c4c4c4",
          alpha = input[["expression_projection_dot_opacity"]]
        ) +
        lims(x = xlim, y = ylim) +
        theme_bw()

        if ( input[["expression_projection_color_scale"]] == 'viridis' ) {
          p <- p + viridis::scale_fill_viridis(
            option = "viridis",
            limits = input[["expression_projection_color_scale_range"]],
            oob = scales::squish,
            direction = -1,
            name = "Log-normalised\nexpression",
            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
          )
        } else {
          p <- p + scale_fill_distiller(
            palette = input[["expression_projection_color_scale"]],
            limits = input[["expression_projection_color_scale_range"]],
            oob = scales::squish,
            direction = 1,
            name = "Log-normalised\nexpression",
            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
          )
        }

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

output[["expression_details_selected_cells"]] <- DT::renderDataTable(server = FALSE, {
  ## if no selection has been made, return empty table
  if ( is.null(plotly::event_data("plotly_selected", source = "expression_projection")) ) {
    #table <- gene_expression_plot_data()
    #print('no selection made')
    table <- tibble(
      cell_barcode = character(),
      level = numeric(),
      sample = character(),
      cluster = character(),
      nUMI = numeric(),
      nGene = numeric()
    )
  ## if no cells are selected, return empty table
  } else if ( length(plotly::event_data("plotly_selected", source = "expression_projection")) == 0 ) {
    #table <- gene_expression_plot_data()
    #print('selection made but no cells inside')
    table <- tibble(
      cell_barcode = character(),
      level = numeric(),
      sample = character(),
      cluster = character(),
      nUMI = numeric(),
      nGene = numeric()
    )
  # #if at least 1 cell has been selected
  } else {
    #print('selection made and some cells inside')
    ## get info of selected cells and create identifier from X-Y coordinates
    selected_cells <- plotly::event_data("plotly_selected", source = "expression_projection") %>%
      dplyr::mutate(identifier = paste0(x, '-', y))
    ## filter out non-selected cells with X-Y identifier and select some meta
    ## data
    table <- gene_expression_plot_data() %>%
      dplyr::rename(X1 = 1, X2 = 2) %>%
      dplyr::mutate(identifier = paste0(X1, '-', X2)) %>%
      dplyr::filter(identifier %in% selected_cells$identifier) %>%
      dplyr::select(cell_barcode, level, sample, cluster, nUMI, nGene)
    ## if no cells match the selection (e.g. when changing dimensional
    ## reduction), return empty table
    if ( nrow(table) == 0 ) {
      #print('selection made and some cells inside but filtering returned 0 cells')
      table <- tibble(
        cell_barcode = character(),
        level = numeric(),
        sample = character(),
        cluster = character(),
        nUMI = numeric(),
        nGene = numeric()
      )
    ## if some cells have been selected, format numbers
    } else {
      table <- table %>%
        dplyr::mutate(
          level = round(level, digits = 3),
          nUMI = formattable::comma(nUMI, big.mark = ',', digits = 0),
          nGene = formattable::comma(nGene, big.mark = ',', digits = 0)
        )
    }
  }
  table %>%
  dplyr::rename(
    'Cell barcode' = cell_barcode,
    'Expression of selected genes' = level,
    'Sample' = sample,
    'Cluster' = cluster,
    '# of transcripts' = nUMI,
    '# of expressed genes' = nGene
  ) %>%
  formattable::formattable(list(
    'Expression of selected genes' = formattable::color_tile("white", "orange"),
    '# of transcripts' = formattable::color_tile("white", "orange"),
    '# of expressed genes' = formattable::color_tile("white", "orange")
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
              filename = "gene_expression_details_of_selected_cells",
              title = "Gene expression details of selected cells"
            ),
            list(
              extend = "excel",
              filename = "gene_expression_details_of_selected_cells",
              title = "Gene expression details of selected cells"
            )
          )
        )
      )
    )
  ) %>%
  DT::formatStyle(
    columns = c('Expression of selected genes', '# of transcripts', '# of expressed genes'),
    textAlign = 'right'
  ) %>%
  DT::formatStyle(
    columns = 'Sample',
    textAlign = 'center'#,
    # backgroundColor = DT::styleEqual(
    #   names(reactive_colors()$samples),
    #   reactive_colors()$samples
    # )
  ) %>%
    DT::formatStyle(
    columns = 'Cluster',
    textAlign = 'center'#,
    # color = DT::styleEqual(
    #   names(reactive_colors()$clusters),
    #   reactive_colors()$clusters
    # )
  )
})

# info box
observeEvent(input[["expression_details_selected_cells_info"]], {
  showModal(
    modalDialog(
      expression_details_selected_cells_info$text,
      title = expression_details_selected_cells_info$title,
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Expression in selected cells.
##----------------------------------------------------------------------------##

# violin + box plot
output[["expression_in_selected_cells"]] <- plotly::renderPlotly({
  if (
    is.null(plotly::event_data("plotly_selected", source = "expression_projection")) |
    length(plotly::event_data("plotly_selected", source = "expression_projection")) == 0
  ) {
    data <- gene_expression_plot_data() %>% dplyr::mutate(group = 'not selected')
  } else {
    selected_cells <- plotly::event_data("plotly_selected", source = "expression_projection") %>%
      dplyr::mutate(identifier = paste0(x, '-', y))
    data <- gene_expression_plot_data() %>%
      dplyr::rename(X1 = 1, X2 = 2) %>%
      dplyr::mutate(
        identifier = paste0(X1, '-', X2),
        group = ifelse(identifier %in% selected_cells$identifier, 'selected', 'not selected'),
        group = factor(group, levels = c('selected', 'not selected'))
      ) %>%
      dplyr::select(group, level)
  }
  plotly::plot_ly(
    data,
    x = ~group,
    y = ~level,
    type = "violin",
    box = list(
      visible = TRUE
    ),
    meanline = list(
      visible = TRUE
    ),
    color = ~group,
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
      title = "Expression level",
      range = c(0, max(data$level) * 1.2),
      hoverformat = ".2f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# info box
observeEvent(input[["expression_in_selected_cells_info"]], {
  showModal(
    modalDialog(
      expression_in_selected_cells_info$text,
      title = expression_in_selected_cells_info$title,
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Expression by sample.
##----------------------------------------------------------------------------##

# box plot
output[["expression_by_sample"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    gene_expression_plot_data(),
    x = ~sample,
    y = ~level,
    type = "violin",
    box = list(
      visible = TRUE
    ),
    meanline = list(
      visible = TRUE
    ),
    color = ~sample,
    colors = reactive_colors()$samples,
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
      title = "Expression level",
      range = c(0, max(gene_expression_plot_data()$level) * 1.2),
      hoverformat = ".2f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# info box
observeEvent(input[["expression_by_sample_info"]], {
  showModal(
    modalDialog(
      expression_by_sample_info$text,
      title = expression_by_sample_info$title,
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Expression by cluster.
##----------------------------------------------------------------------------##

# box plot
output[["expression_by_cluster"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    gene_expression_plot_data(),
    x = ~cluster,
    y = ~level,
    type = "violin",
    box = list(
      visible = TRUE
    ),
    meanline = list(
      visible = TRUE
    ),
    color = ~cluster,
    colors = reactive_colors()$clusters,
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
      title = "Expression level",
      range = c(0, max(gene_expression_plot_data()$level) * 1.2),
      hoverformat = ".2f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode =  "select",
    hovermode = "compare"
  )
})

# info box
observeEvent(input[["expression_by_cluster_info"]], {
  showModal(
    modalDialog(
      expression_by_cluster_info$text,
      title = expression_by_cluster_info$title,
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Expression by gene.
##----------------------------------------------------------------------------##

# bar plot
output[["expression_by_gene"]] <- plotly::renderPlotly({
  req(
    input[["expression_projection_color_scale"]]
  )
  if ( length(genesToPlot()$genes_to_display_present) == 0 ) {
    expression_levels <- data.frame(
      "gene" = character(),
      "expression" = integer()
    )
  } else if ( length(genesToPlot()$genes_to_display_present) == 1 ) {
    expression_levels <- data.frame(
      "gene" = genesToPlot()$genes_to_display_present,
      "expression" = mean(sample_data()$expression[ genesToPlot()$genes_to_display_present , ])
    )
  } else {
    expression_levels <- data.frame(
      "gene" = rownames(sample_data()$expression[ genesToPlot()$genes_to_display_present , ]),
      "expression" = Matrix::rowMeans(sample_data()$expression[ genesToPlot()$genes_to_display_present , ])
    ) %>%
    arrange(-expression) %>%
    top_n(50, expression)
  }
  # color scale
  if ( input[["expression_projection_color_scale"]] == 'viridis' ) {
    color_scale <- 'Viridis'
  } else {
    color_scale <- input[["expression_projection_color_scale"]]
  }
  plotly::plot_ly(
    expression_levels,
    x = ~gene,
    y = ~expression,
    text = ~paste0(
      expression_levels$gene, ': ',
      format(expression_levels$expression, digits = 3)
    ),
    type = "bar",
    marker = list(
      color = ~expression,
      colorscale = color_scale,
      reversescale = TRUE,
      line = list(
        color = "rgb(196,196,196)",
        width = 1
      )
    ),
    hoverinfo = "text",
    showlegend = FALSE
  ) %>%
  plotly::layout(
    title = "",
    xaxis = list(
      title = "",
      type = "category",
      categoryorder = "array",
      categoryarray = expression_levels$gene,
      mirror = TRUE,
      showline = TRUE
    ),
    yaxis = list(
      title = "Expression level",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# info box
observeEvent(input[["expression_by_gene_info"]], {
  showModal(
    modalDialog(
      expression_by_gene_info[["text"]],
      title = expression_by_gene_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})
