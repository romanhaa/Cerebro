##----------------------------------------------------------------------------##
## Tab: Gene expression
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## UI elements for projection.
##----------------------------------------------------------------------------##
output[["expression_UI"]] <- renderUI({
  tagList(
    textAreaInput(
      "expression_genes_input",
      label = "Gene(s):",
      value = "",
      placeholder = "Insert genes here.",
      resize = "vertical"
    ),
    selectInput(
      "expression_projection_to_display",
      label = "Projection:",
      choices = names(sample_data()$projections)
    ),
    shinyWidgets::pickerInput(
      "expression_samples_to_display",
      label = "Samples to display:",
      choices = sample_data()$sample_names,
      selected = sample_data()$sample_names,
      options = list("actions-box" = TRUE),
      multiple = TRUE
    ),
    shinyWidgets::pickerInput(
      "expression_clusters_to_display",
      label = "Clusters to display:",
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
      label = "Plotting order:",
      choices = c("Random", "Highest expression on top"),
      selected = "Random"
    ),
    sliderInput(
      "expression_projection_dot_size",
      label = "Point size:",
      min = scatter_plot_dot_size[["min"]],
      max = scatter_plot_dot_size[["max"]],
      step = scatter_plot_dot_size[["step"]],
      value = scatter_plot_dot_size[["default"]]
    ),
    sliderInput(
      "expression_projection_dot_opacity",
      label = "Point opacity:",
      min = scatter_plot_dot_opacity[["min"]],
      max = scatter_plot_dot_opacity[["max"]],
      step = scatter_plot_dot_opacity[["step"]],
      value = scatter_plot_dot_opacity[["default"]]
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
      label = "X axis",
      min = range_x_min,
      max = range_x_max,
      value = c(range_x_min, range_x_max)
    ),
    sliderInput(
      "expression_projection_scale_y_manual_range",
      label = "Y axis",
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
genesToPlot <- eventReactive(input[["keyPressed"]], ignoreNULL = FALSE, {
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

output[["expression_projection_plotly"]] <- plotly::renderPlotly({
  req(
    input[["expression_projection_dot_opacity"]],
    input[["expression_projection_dot_size"]],
    input[["expression_projection_scale_x_manual_range"]],
    input[["expression_projection_scale_y_manual_range"]]
  )
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
        colorscale = "YlGnBu",
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
      )
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
        colorscale = "YlGnBu",
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
      )
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
        scale_fill_distiller(
          palette = "YlGnBu",
          direction = 1,
          name = "Log-normalised\nexpression",
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
        ) +
        lims(x = xlim, y = ylim) +
        theme_bw()

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
    colors = sample_data()$samples$colors,
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
    colors = sample_data()$clusters$colors,
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
  plotly::plot_ly(
    expression_levels,
    x = ~gene,
    y = ~expression,
    text = ~gene,
    type = "bar",
    marker = list(
      color = ~expression,
      colorscale = "YlGnBu",
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


