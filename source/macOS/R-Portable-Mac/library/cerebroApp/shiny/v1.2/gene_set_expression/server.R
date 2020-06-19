##----------------------------------------------------------------------------##
## Tab: Gene set expression
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## Reactive data: Gene sets available.
##----------------------------------------------------------------------------##
geneSets <- reactive({
  if ( sample_data()$experiment$organism == "mm" ) {
    msigdbr::msigdbr(species = "Mus musculus")
  } else if ( sample_data()$experiment$organism == "hg" ) {
    msigdbr::msigdbr(species = "Homo sapiens")
  } else {
    msigdbr::msigdbr(species = "Mus musculus")
  }
})

##----------------------------------------------------------------------------##
## UI element: Which genes are available in data set.
##----------------------------------------------------------------------------##
output[["geneSetExpression_genes_displayed"]] <- renderText({
  genes_to_display_text <- paste0(
      "<br><b>Total unique genes in gene set:</b><br>",
      length(geneSetData()$genes_to_display),
      "<br><b>Showing expression for ",
      length(geneSetData()$genes_to_display_present),
      " genes:</b><br>",
      paste0(geneSetData()$genes_to_display_present, collapse = ", "),
      "<br><b>",
      length(geneSetData()$genes_to_display_missing),
      " gene(s) are not in data set: </b><br>",
      paste0(geneSetData()$genes_to_display_missing, collapse = ", ")
    )
  if ( sample_data()$experiment$organism != "mm" & sample_data()$experiment$organism != "hg" ) {
    paste0(
      '<br><b><font color="red">Note:</b> Gene sets are available for human and mouse only. Organism for loaded samples is either not set or none of the two. Mouse gene sets are loaded and can be used.</font><br>',
      genes_to_display_text
    )
  } else {
    genes_to_display_text
  }
})

##----------------------------------------------------------------------------##
## UI element for projection.
##----------------------------------------------------------------------------##
output[["geneSetExpression_UI"]] <- renderUI({
  tagList(
    selectizeInput(
      'geneSetExpression_select_geneSet',
      label = 'Gene set',
      choices = c("-", unique(geneSets()$gs_name)),
      options = list(create = TRUE), multiple = FALSE, selected = "-"
    ),
    selectInput(
      "geneSetExpression_projection_to_display",
      label = "Projection",
      choices = names(sample_data()$projections)
    ),
    shinyWidgets::pickerInput(
      "geneSetExpression_samples_to_display",
      label = "Samples to display",
      choices = sample_data()$sample_names,
      selected = sample_data()$sample_names,
      options = list("actions-box" = TRUE),
      multiple = TRUE
    ),
    shinyWidgets::pickerInput(
      "geneSetExpression_clusters_to_display",
      label = "Clusters to display",
      choices = sample_data()$cluster_names,
      selected = sample_data()$cluster_names,
      options = list("actions-box" = TRUE),
      multiple = TRUE
    ),
    sliderInput(
      "geneSetExpression_percentage_cells_to_show",
      label = "Show % of cells",
      min = scatter_plot_percentage_cells_to_show[["min"]],
      max = scatter_plot_percentage_cells_to_show[["max"]],
      step = scatter_plot_percentage_cells_to_show[["step"]],
      value = scatter_plot_percentage_cells_to_show[["default"]]
    ),
    selectInput(
      "geneSetExpression_projection_plotting_order",
      label = "Plotting order",
      choices = c("Random", "Highest expression on top")
    ),
    sliderInput(
      "geneSetExpression_projection_dot_size",
      label = "Dot size",
      min = scatter_plot_dot_size[["min"]],
      max = scatter_plot_dot_size[["max"]],
      step = scatter_plot_dot_size[["step"]],
      value = scatter_plot_dot_size[["default"]]
    ),
    sliderInput(
      "geneSetExpression_projection_dot_opacity",
      label = "Dot opacity",
      min = scatter_plot_dot_opacity[["min"]],
      max = scatter_plot_dot_opacity[["max"]],
      step = scatter_plot_dot_opacity[["step"]],
      value = scatter_plot_dot_opacity[["default"]]
    ),
    selectInput(
      "geneSetExpression_projection_color_scale",
      label = "Color scale",
      choices = c("YlGnBu", "YlOrRd","Blues","Greens","Reds","RdBu","viridis"),
      selected = "YlGnBu"
    )
  )
})

##----------------------------------------------------------------------------##
## UI element for color scale range in projection.
##----------------------------------------------------------------------------##
output[["geneSetExpression_color_scale_range"]] <- renderUI({
  range <- range(geneSetExpression_plot_data()$level)
  if ( range[1] == 0 & range[2] == 0 ) {
    range[2] = 1
  } else {
    range[1] <- range[1] %>% round(digits = 2)
    range[2] <- range[2] %>% round(digits = 2)
  }
  tagList(
    sliderInput(
      "geneSetExpression_projection_color_scale_range",
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
output[["geneSetExpression_scales"]] <- renderUI({
  req(input[["geneSetExpression_projection_to_display"]])
  projection_to_display <- input[["geneSetExpression_projection_to_display"]]
  range_x_min <- sample_data()$projections[[ projection_to_display ]][,1] %>% min() %>% "*"(ifelse(.<0, 1.1, 0.9)) %>% round()
  range_x_max <- sample_data()$projections[[ projection_to_display ]][,1] %>% max() %>% "*"(ifelse(.<0, 0.9, 1.1)) %>% round()
  range_y_min <- sample_data()$projections[[ projection_to_display ]][,2] %>% min() %>% "*"(ifelse(.<0, 1.1, 0.9)) %>% round()
  range_y_max <- sample_data()$projections[[ projection_to_display ]][,2] %>% max() %>% "*"(ifelse(.<0, 0.9, 1.1)) %>% round()
  tagList(
    sliderInput(
      "geneSetExpression_projection_scale_x_manual_range",
      label = "Range of X axis",
      min = range_x_min,
      max = range_x_max,
      value = c(range_x_min, range_x_max)
    ),
    sliderInput(
      "geneSetExpression_projection_scale_y_manual_range",
      label = "Range of Y axis",
      min = range_y_min,
      max = range_y_max,
      value = c(range_y_min, range_y_max)
    )
  )
})

##----------------------------------------------------------------------------##
## Reactive data: Genes of gene set in data set.
##----------------------------------------------------------------------------##
geneSetData <- reactive({
  req(input[["geneSetExpression_select_geneSet"]])
  geneSetData <- list()
  if ( input[["geneSetExpression_select_geneSet"]] == "-" ) {
    geneSetData$genes_to_display_present <- NULL
  } else {
    geneSetData[["genes_to_display"]] <- geneSets()[ which(geneSets()$gs_name == input[["geneSetExpression_select_geneSet"]]) , "gene_symbol" ]$gene_symbol
    geneSetData[["genes_to_display"]] <- unique(geneSetData[["genes_to_display"]])
    geneSetData[["genes_to_display_here"]] <- rownames(sample_data()$expression)[ match(tolower(geneSetData[["genes_to_display"]]), tolower(rownames(sample_data()$expression))) ]
    geneSetData[["genes_to_display_present"]] <- geneSetData[["genes_to_display_here"]][ which(!is.na(geneSetData[["genes_to_display_here"]])) ]
    geneSetData[["genes_to_display_missing"]] <- geneSetData[["genes_to_display"]][ which(is.na(geneSetData[["genes_to_display_here"]])) ]
  }
  geneSetData
})

##----------------------------------------------------------------------------##
## Reactive data: Data to plot.
##----------------------------------------------------------------------------##
geneSetExpression_plot_data <- reactive({
  req(
    input[["geneSetExpression_projection_to_display"]],
    input[["geneSetExpression_samples_to_display"]],
    input[["geneSetExpression_clusters_to_display"]],
    input[["geneSetExpression_percentage_cells_to_show"]],
    input[["geneSetExpression_projection_plotting_order"]]
  )
  projection_to_display <- input[["geneSetExpression_projection_to_display"]]
  samples_to_display <- input[["geneSetExpression_samples_to_display"]]
  clusters_to_display <- input[["geneSetExpression_clusters_to_display"]]
  percentage_cells_show <- input[["geneSetExpression_percentage_cells_to_show"]]
  plot_order <- input[["geneSetExpression_projection_plotting_order"]]
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
  if ( length(geneSetData()$genes_to_display_present) == 0 ) {
    plot$level <- 0
  } else if ( length(geneSetData()$genes_to_display_present) == 1 ) {
    plot$level <- geneSetData()$genes_to_display_present %>%
      sample_data()$expression[ . , cells_to_display ]
  } else {
    plot$level <- geneSetData()$genes_to_display_present %>%
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
output[["geneSetExpression_projection"]] <- plotly::renderPlotly({
  req(
    input[["geneSetExpression_projection_to_display"]],
    input[["geneSetExpression_projection_dot_size"]],
    input[["geneSetExpression_projection_dot_opacity"]],
    input[["geneSetExpression_projection_color_scale"]],
    input[["geneSetExpression_projection_color_scale_range"]],
    input[["geneSetExpression_projection_scale_x_manual_range"]],
    input[["geneSetExpression_projection_scale_y_manual_range"]]
  )
  if ( input[["geneSetExpression_projection_color_scale"]] == 'viridis' ) {
    color_scale <- 'Viridis'
  } else {
    color_scale <- input[["geneSetExpression_projection_color_scale"]]
  }
  if ( ncol(sample_data()$projections[[ input[["geneSetExpression_projection_to_display"]] ]]) == 3 ) {
    plotly::plot_ly(
      geneSetExpression_plot_data(),
      x = geneSetExpression_plot_data()[,1],
      y = geneSetExpression_plot_data()[,2],
      z = geneSetExpression_plot_data()[,3],
      type = "scatter3d",
      mode = "markers",
      marker = list(
        colorbar = list(
          title = "Expression"
        ),
        color = ~level,
        opacity = input[["geneSetExpression_projection_dot_opacity"]],
        colorscale = color_scale,
        cauto = FALSE,
        cmin = input[["geneSetExpression_projection_color_scale_range"]][1],
        cmax = input[["geneSetExpression_projection_color_scale_range"]][2],
        reversescale = TRUE,
        line = list(
          color = "rgb(196,196,196)",
          width = 1
        ),
        size = input[["geneSetExpression_projection_dot_size"]]
      ),
      hoverinfo = "text",
      text = ~paste(
        "<b>Cell</b>: ", geneSetExpression_plot_data()$cell_barcode, "<br>",
        "<b>Sample</b>: ", geneSetExpression_plot_data()$sample, "<br>",
        "<b>Cluster</b>: ", geneSetExpression_plot_data()$cluster, "<br>",
        "<b>Transcripts</b>: ", formatC(geneSetExpression_plot_data()$nUMI, format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>Expressed genes</b>: ", formatC(geneSetExpression_plot_data()$nGene, format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>Expression level</b>: ", formatC(geneSetExpression_plot_data()$level, format = "f", big.mark = ",", digits = 3)
      ),
      source = "geneSetExpression_projection"
    ) %>%
    plotly::layout(
      scene = list(
        xaxis = list(
          title = colnames(geneSetExpression_plot_data())[1],
          mirror = TRUE,
          showline = TRUE,
          zeroline = FALSE
        ),
        yaxis = list(
          title = colnames(geneSetExpression_plot_data())[2],
          mirror = TRUE,
          showline = TRUE,
          zeroline = FALSE
        ),
        zaxis = list(
          title = colnames(geneSetExpression_plot_data())[3],
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
      geneSetExpression_plot_data(),
      x = geneSetExpression_plot_data()[,1],
      y = geneSetExpression_plot_data()[,2],
      type = "scatter",
      mode = "markers",
      marker = list(
        colorbar = list(
          title = "Expression"
        ),
        color = ~level,
        opacity = input[["geneSetExpression_projection_dot_opacity"]],
        colorscale = color_scale,
        cauto = FALSE,
        cmin = input[["geneSetExpression_projection_color_scale_range"]][1],
        cmax = input[["geneSetExpression_projection_color_scale_range"]][2],
        reversescale = TRUE,
        line = list(
          color = "rgb(196,196,196)",
          width = 1
        ),
        size = input[["geneSetExpression_projection_dot_size"]]
      ),
      hoverinfo = "text",
      text = ~paste(
        "<b>Cell</b>: ", geneSetExpression_plot_data()$cell_barcode, "<br>",
        "<b>Sample</b>: ", geneSetExpression_plot_data()$sample, "<br>",
        "<b>Cluster</b>: ", geneSetExpression_plot_data()$cluster, "<br>",
        "<b>Transcripts</b>: ", formatC(geneSetExpression_plot_data()$nUMI, format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>Expressed genes</b>: ", formatC(geneSetExpression_plot_data()$nGene, format = "f", big.mark = ",", digits = 0), "<br>",
        "<b>Expression level</b>: ", formatC(geneSetExpression_plot_data()$level, format = "f", big.mark = ",", digits = 3)
      ),
      source = "geneSetExpression_projection"
    ) %>%
    plotly::layout(
      xaxis = list(
        title = colnames(geneSetExpression_plot_data())[1],
        mirror = TRUE,
        showline = TRUE,
        zeroline = FALSE,
        range = c(
          input[["geneSetExpression_projection_scale_x_manual_range"]][1],
          input[["geneSetExpression_projection_scale_x_manual_range"]][2]
        )
      ),
      yaxis = list(
        title = colnames(geneSetExpression_plot_data())[2],
        mirror = TRUE,
        showline = TRUE,
        zeroline = FALSE,
        range = c(
          input[["geneSetExpression_projection_scale_y_manual_range"]][1],
          input[["geneSetExpression_projection_scale_y_manual_range"]][2]
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
observeEvent(input[["geneSetExpression_projection_info"]], {
  showModal(
    modalDialog(
      geneSetExpression_projection_info[["text"]],
      title = geneSetExpression_projection_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Export function.
##----------------------------------------------------------------------------##
observeEvent(input[["geneSetExpression_projection_export"]], {
  req(
    input[["geneSetExpression_projection_to_display"]],
    input[["geneSetExpression_projection_plotting_order"]],
    input[["geneSetExpression_projection_dot_size"]],
    input[["geneSetExpression_projection_dot_opacity"]],
    input[["geneSetExpression_projection_color_scale"]],
    input[["geneSetExpression_projection_color_scale_range"]],
    input[["geneSetExpression_projection_scale_x_manual_range"]],
    input[["geneSetExpression_projection_scale_y_manual_range"]]
  )
  library("ggplot2")
  if ( exists("plot_export_path") ) {
    xlim <- c(
      input[["geneSetExpression_projection_scale_x_manual_range"]][1],
      input[["geneSetExpression_projection_scale_x_manual_range"]][2]
    )
    ylim <- c(
      input[["geneSetExpression_projection_scale_y_manual_range"]][1],
      input[["geneSetExpression_projection_scale_y_manual_range"]][2]
    )

    out_filename <- paste0(
      plot_export_path, "Cerebro_", sample_data()$experiment$experiment_name,
      "_gene_set_expression_", input[["geneSetExpression_select_geneSet"]], "_",
      input[["geneSetExpression_projection_to_display"]]
    )

    if ( input[["geneSetExpression_projection_plotting_order"]] == "Random" ) {
      out_filename <- paste0(out_filename, "_random_order.pdf")
    } else if ( input[["geneSetExpression_projection_plotting_order"]] == "Highest expression on top" ) {
      out_filename <- paste0(out_filename, "_highest_expression_on_top.pdf")
    }

    if ( ncol(sample_data()$projections[[ input[["geneSetExpression_projection_to_display"]] ]]) == 3 ) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Sorry!",
        text = "It's currently not possible to create PDF plots from 3D dimensional reductions. Please use the PNG export button in the panel or a 2D dimensional reduction instead.",
        type = "error"
      )
    } else {
      p <- ggplot(
          geneSetExpression_plot_data(),
          aes_q(
            x = as.name(colnames(geneSetExpression_plot_data())[1]),
            y = as.name(colnames(geneSetExpression_plot_data())[2]),
            fill = as.name("level")
          )
        ) +
        geom_point(
          shape = 21,
          size = input[["geneSetExpression_projection_dot_size"]]/3,
          stroke = 0.2,
          color = "#c4c4c4",
          alpha = input[["geneSetExpression_projection_dot_opacity"]]
        ) +
        lims(x = xlim, y = ylim) +
        theme_bw()

        if ( input[["geneSetExpression_projection_color_scale"]] == 'viridis' ) {
          p <- p + viridis::scale_fill_viridis(
            option = "viridis",
            limits = input[["geneSetExpression_projection_color_scale_range"]],
            oob = scales::squish,
            direction = -1,
            name = "Log-normalised\nexpression",
            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
          )
        } else {
          p <- p + scale_fill_distiller(
            palette = input[["geneSetExpression_projection_color_scale"]],
            limits = input[["geneSetExpression_projection_color_scale_range"]],
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

output[["geneSetExpression_details_selected_cells"]] <- DT::renderDataTable(server = FALSE, {
  if ( is.null(plotly::event_data("plotly_selected", source = "geneSetExpression_projection")) ) {
    table <- tibble(
      cell_barcode = character(),
      level = numeric(),
      sample = character(),
      cluster = character(),
      nUMI = numeric(),
      nGene = numeric()
    )
  } else if ( length(plotly::event_data("plotly_selected", source = "geneSetExpression_projection")) == 0 ) {
    table <- tibble(
      cell_barcode = character(),
      level = numeric(),
      sample = character(),
      cluster = character(),
      nUMI = numeric(),
      nGene = numeric()
    )
  } else {
    selected_cells <- plotly::event_data("plotly_selected", source = "geneSetExpression_projection") %>%
      dplyr::mutate(identifier = paste0(x, '-', y))
    table <- geneSetExpression_plot_data() %>%
      dplyr::rename(X1 = 1, X2 = 2) %>%
      dplyr::mutate(identifier = paste0(X1, '-', X2)) %>%
      dplyr::filter(identifier %in% selected_cells$identifier) %>%
      dplyr::select(cell_barcode, level, sample, cluster, nUMI, nGene)
    if ( nrow(table) == 0 ) {
      table <- tibble(
        cell_barcode = character(),
        level = numeric(),
        sample = character(),
        cluster = character(),
        nUMI = numeric(),
        nGene = numeric()
      )
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
    'Expression of selected gene set' = level,
    'Sample' = sample,
    'Cluster' = cluster,
    '# of transcripts' = nUMI,
    '# of expressed genes' = nGene
  ) %>%
  formattable::formattable(list(
    'Expression of selected gene set' = formattable::color_tile("white", "orange"),
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
              filename = "gene_set_expression_details_of_selected_cells",
              title = "Gene set expression details of selected cells"
            ),
            list(
              extend = "excel",
              filename = "gene_set_expression_details_of_selected_cells",
              title = "Gene set expression details of selected cells"
            )
          )
        )
      )
    )
  ) %>%
  DT::formatStyle(
    columns = c('Expression of selected gene set', '# of transcripts', '# of expressed genes'),
    textAlign = 'right'
  ) %>%
  DT::formatStyle(
    columns = c('Sample','Cluster'),
    textAlign = 'center'
  )
})

# info box
observeEvent(input[["geneSetExpression_details_selected_cells_info"]], {
  showModal(
    modalDialog(
      geneSetExpression_details_selected_cells_info$text,
      title = geneSetExpression_details_selected_cells_info$title,
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Expression in selected cells.
##----------------------------------------------------------------------------##

# violin + box plot
output[["geneSetExpression_in_selected_cells"]] <- plotly::renderPlotly({
  if (
    is.null(plotly::event_data("plotly_selected", source = "geneSetExpression_projection")) |
    length(plotly::event_data("plotly_selected", source = "geneSetExpression_projection")) == 0
  ) {
    data <- geneSetExpression_plot_data() %>% dplyr::mutate(group = 'not selected')
  } else {
    selected_cells <- plotly::event_data("plotly_selected", source = "geneSetExpression_projection") %>%
      dplyr::mutate(identifier = paste0(x, '-', y))
    data <- geneSetExpression_plot_data() %>%
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
observeEvent(input[["geneSetExpression_in_selected_cells_info"]], {
  showModal(
    modalDialog(
      geneSetExpression_in_selected_cells_info$text,
      title = geneSetExpression_in_selected_cells_info$title,
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Expression by sample.
##----------------------------------------------------------------------------##
# box plot
output[["geneSetExpression_by_sample"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    geneSetExpression_plot_data(),
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
      range = c(0, max(geneSetExpression_plot_data()$level) * 1.2),
      hoverformat=".2f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# info box
observeEvent(input[["geneSetExpression_by_sample_info"]], {
  showModal(
    modalDialog(
      geneSetExpression_by_sample_info[["text"]],
      title = geneSetExpression_by_sample_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Expression by cluster.
##----------------------------------------------------------------------------##
# box plot
output[["geneSetExpression_by_cluster"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    geneSetExpression_plot_data(),
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
      range = c(0, max(geneSetExpression_plot_data()$level) * 1.2),
      hoverformat = ".2f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# info box
observeEvent(input[["geneSetExpression_by_cluster_info"]], {
  showModal(
    modalDialog(
      geneSetExpression_by_cluster_info[["text"]],
      title = geneSetExpression_by_cluster_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})


##----------------------------------------------------------------------------##
## Expression by gene.
##----------------------------------------------------------------------------##

# bar plot
output[["geneSetExpression_by_gene"]] <- plotly::renderPlotly({
  req(
    input[["geneSetExpression_projection_color_scale"]]
  )
  if ( length(geneSetData()$genes_to_display_present) == 0 ) {
    expression_levels <- data.frame(
      "gene" = character(),
      "expression" = integer()
    )
  } else if ( length(geneSetData()$genes_to_display_present) == 1 ) {
    expression_levels <- data.frame(
      "gene" = geneSetData()$genes_to_display_present,
      "expression" = mean(sample_data()$expression[ geneSetData()$genes_to_display_present , ])
    )
  } else {
    expression_levels <- data.frame(
      "gene" = rownames(sample_data()$expression[ geneSetData()$genes_to_display_present , ]),
      "expression" = Matrix::rowMeans(sample_data()$expression[ geneSetData()$genes_to_display_present , ])
    ) %>%
    arrange(-expression) %>%
    top_n(50, expression)
  }
  # color scale
  if ( input[["geneSetExpression_projection_color_scale"]] == 'viridis' ) {
    color_scale <- 'Viridis'
  } else {
    color_scale <- input[["geneSetExpression_projection_color_scale"]]
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
observeEvent(input[["geneSetExpression_by_gene_info"]], {
  showModal(
    modalDialog(
      expression_by_gene_info[["text"]],
      title = expression_by_gene_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})
