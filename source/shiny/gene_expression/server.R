##--------------------------------------------------------------------------##
## Panel: Gene expression
##--------------------------------------------------------------------------##

# reactive data
genesToPlot <- reactive({
  genesToPlot <- list()
  if ( is.null(input$expression_genes_input) ) {
    genesToPlot$genes_to_display <- ""
  } else {
    genesToPlot$genes_to_display <- input$expression_genes_input %>%
      strsplit(",| |;|\n") %>%
      unlist() %>%
      gsub(pattern = " ", replacement = "", fixed = TRUE) %>%
      unique()
  }
  genesToPlot$genes_to_display_here <- rownames(sample_data()$expression)[ match(tolower(genesToPlot$genes_to_display), tolower(rownames(sample_data()$expression))) ]
  genesToPlot$genes_to_display_present <- genesToPlot$genes_to_display_here[ which(!is.na(genesToPlot$genes_to_display_here)) ]
  genesToPlot$genes_to_display_missing <- genesToPlot$genes_to_display[ which(is.na(genesToPlot$genes_to_display_here)) ]
  genesToPlot
})

# select genes to be displayed
output$expression_genes_displayed <- renderText({
  paste0(
    "<b>Showing expression for ",
    length(genesToPlot()$genes_to_display_present), " gene(s):</b><br>",
    paste0(genesToPlot()$genes_to_display_present, collapse = ", "),
    "<br><b>",
    length(genesToPlot()$genes_to_display_missing),
    " gene(s) are not in data set: </b><br>",
    paste0(genesToPlot()$genes_to_display_missing, collapse = ", ")
  )
})

##--------------------------------------------------------------------------##
# UI
output$expression_UI <- renderUI({
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
      choices = sample_data()$samples$overview$sample,
      selected = sample_data()$samples$overview$sample,
      options = list("actions-box" = TRUE),
      multiple = TRUE
    ),
    shinyWidgets::pickerInput(
      "expression_clusters_to_display",
      label = "Clusters to display:",
      choices = sample_data()$clusters$overview$cluster,
      selected = sample_data()$clusters$overview$cluster,
      options = list("actions-box" = TRUE),
      multiple = TRUE
    ),
    sliderInput(
      "expression_percentage_cells_to_show",
      label = "Show % of cells",
      min = 10,
      max = 100,
      step = 10,
      value = 100
    ),
    selectInput("expression_plotting_order", label = "Plotting order:",
      choices = c("Random", "Highest expression on top")
    ),
    sliderInput("expression_projection_dot_size", label = "Point size:",
      min = 0, max = 50, value = 25, step = 1
    ),
    sliderInput("expression_projection_opacity", label = "Point opacity:",
      min = 0, max = 1, value = 1, step = 0.05
    )
  )
})

output$expression_scales <- renderUI({
  projection_to_display <- if ( is.null(input$expression_projection_to_display) || is.na(input$expression_projection_to_display) ) names(sample_data()$projections)[1] else input$expression_projection_to_display
  range_x_min <- round(
      min(sample_data()$projections[[ projection_to_display ]][,1]) * 1.1)
  range_x_max <- round(
      max(sample_data()$projections[[ projection_to_display ]][,1]) * 1.1)
  range_y_min <- round(
      min(sample_data()$projections[[ projection_to_display ]][,2]) * 1.1)
  range_y_max <- round(
      max(sample_data()$projections[[ projection_to_display ]][,2]) * 1.1)
  tagList(
    sliderInput(
      "expression_projection_scale_x_manual_range",
      label = "X axis",
      min = range_x_min, max = range_x_max,
      value = c(range_x_min, range_x_max)
    ),
    sliderInput(
      "expression_projection_scale_y_manual_range",
      label = "Y axis",
      min = range_y_min, max = range_y_max,
      value = c(range_y_min, range_y_max)
    )
  )
})

##--------------------------------------------------------------------------##

# projection with scatterD3
output$expression_projection <- scatterD3::renderScatterD3({
  req(input$expression_projection_to_display)
  req(input$expression_samples_to_display)
  req(input$expression_clusters_to_display)
  req(input$expression_plotting_order)
  projection_to_display <- input$expression_projection_to_display
  samples_to_display <- input$expression_samples_to_display
  clusters_to_display <- input$expression_clusters_to_display
  cells_to_display <- which(
      grepl(sample_data()$cells$sample, pattern = paste0("^", samples_to_display, "$", collapse = "|")) == TRUE & 
      grepl(sample_data()$cells$cluster, pattern = paste0("^", clusters_to_display, "$", collapse = "|")) == TRUE
    )
  # randomly remove cells
  if ( input$expression_percentage_cells_to_show < 100 ) {
    number_of_cells_to_plot <- ceiling(
      input$expression_percentage_cells_to_show / 100 * length(cells_to_display)
    )
    cells_to_display <- cells_to_display[ sample(1:length(cells_to_display), number_of_cells_to_plot) ]
  }
  to_plot <- cbind(
      sample_data()$projections[[ projection_to_display ]][ cells_to_display , ],
      sample_data()$cells[ cells_to_display , ]
    )
  if ( length(genesToPlot()$genes_to_display_present) == 0 ) {
    to_plot$level <- 0
  } else if ( length(genesToPlot()$genes_to_display_present) == 1 ) {
    to_plot$level <- genesToPlot()$genes_to_display_present %>%
      sample_data()$expression[ . , cells_to_display ]
  } else {
    to_plot$level <- genesToPlot()$genes_to_display_present %>%
      sample_data()$expression[ . , cells_to_display ] %>%
      colMeans() %>%
      as.vector()
  }
  if ( input$expression_plotting_order == "Random" ) {
    to_plot <- sample(1:nrow(to_plot), nrow(to_plot)) %>% to_plot[ . , ]
  } else if ( input$expression_plotting_order == "Highest expression on top" ) {
    to_plot <- to_plot[ order(to_plot$level, decreasing = FALSE) , ]
  }
  scatterD3::scatterD3(
    x = to_plot[ , 1 ],
    y = to_plot[ , 2 ],
    xlab = colnames(to_plot)[ 1 ],
    ylab = colnames(to_plot)[ 2 ],
    xlim = c(
        input$expression_projection_scale_x_manual_range[1],
        input$expression_projection_scale_x_manual_range[2]
      ),
    ylim = c(
        input$expression_projection_scale_y_manual_range[1],
        input$expression_projection_scale_y_manual_range[2]
      ),
    point_size = input$expression_projection_dot_size,
    col_var = to_plot$level,
    col_lab = "Gene expression",
    col_continuous = TRUE,
    point_opacity = input$expression_projection_opacity,
    transitions = FALSE,
    legend_width = 0,
    menu = FALSE,
    tooltip_text = paste0(
      "<b>Sample</b>: ", to_plot[ , "sample" ], "<br/>",
      "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br/>",
      "<b>nUMI</b>: ", to_plot[ , "nUMI" ], "<br/>",
      "<b>nGene</b>: ", to_plot[ , "nGene" ], "<br/>",
      # "<b>Expr. MT</b>: ", format(to_plot[ , "percent_mt"   ]*100, digits = 1), "%<br/>",
      # "<b>Expr. ribo</b>: ", format(to_plot[ , "percent_ribo" ]*100, digits = 1), "%<br/>"
    )
  )
})

observeEvent(input$expression_projection_info, {
  showModal(
    modalDialog(
      expression_projection_info$text,
      title = expression_projection_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})

observeEvent(input$expression_projection_export, {
  library("ggplot2")

  projection_to_display <- input$expression_projection_to_display
  samples_to_display <- input$expression_samples_to_display
  clusters_to_display <- input$expression_clusters_to_display
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

  xlim <- c(
      input$expression_projection_scale_x_manual_range[1],
      input$expression_projection_scale_x_manual_range[2]
    )
  ylim <- c(
      input$expression_projection_scale_y_manual_range[1],
      input$expression_projection_scale_y_manual_range[2]
    )

  if ( length(genesToPlot()$genes_to_display_present) == 0 ) {
    to_plot$level <- 0
    out_filename <- paste0(
        plot_export_path, "Cerebro_",
        sample_data()$experiment$experiment_name, "_gene_expression_none"
      )
  } else if ( length(genesToPlot()$genes_to_display_present) == 1 ) {
    to_plot$level <- genesToPlot()$genes_to_display_present %>%
      sample_data()$expression[ . , cells_to_display ]
    out_filename <- paste0(
        plot_export_path, "Cerebro_",
        sample_data()$experiment$experiment_name, "_gene_expression_",
        genesToPlot()$genes_to_display_present, "_",
        input$expression_projection_to_display
      )
  } else {
    to_plot$level <- genesToPlot()$genes_to_display_present %>%
      sample_data()$expression[ . , cells_to_display ] %>%
      colMeans() %>%
      as.vector()
    out_filename <- paste0(
        plot_export_path, "Cerebro_",
        sample_data()$experiment$experiment_name, "_gene_expression_",
        genesToPlot()$genes_to_display_present[1],
        "_and_others_", input$expression_projection_to_display
      )
  }

  if ( input$expression_plotting_order == "Random" ) {
    to_plot <- sample(1:nrow(to_plot), nrow(to_plot)) %>% to_plot[ . , ]
    out_filename <- paste0(out_filename, "_random_order.pdf")
  } else if ( input$expression_plotting_order == "Highest expression on top" ) {
    to_plot <- to_plot[ order(to_plot$level, decreasing = FALSE) , ]
    out_filename <- paste0(out_filename, "_highest_expression_on_top.pdf")
  }

  p <- ggplot(
      to_plot,
      aes_q(
        x = as.name(colnames(to_plot)[1]),
        y = as.name(colnames(to_plot)[2]),
        colour = as.name("level")
      )
    ) +
    geom_point() +
    viridis::scale_colour_viridis(
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
})

##--------------------------------------------------------------------------##

# box plot by sample
output$expression_by_sample <- plotly::renderPlotly({
  if ( length(genesToPlot()$genes_to_display_present) == 0 ) {
    expression_levels <- 0
  } else if ( length(genesToPlot()$genes_to_display_present) == 1 ) {
    expression_levels <- genesToPlot()$genes_to_display_present %>%
      sample_data()$expression[ . , ]
  } else {
    expression_levels <- genesToPlot()$genes_to_display_present %>%
      sample_data()$expression[ . , ] %>%
      colMeans() %>%
      as.vector()
  }
  data.frame(
    "sample" = sample_data()$cells[ , "sample" ],
    "expression" = expression_levels
  ) %>%
  plotly::plot_ly(
    x = ~sample,
    y = ~expression,
    type = "box",
    color = ~sample,
    colors = sample_data()$samples$colors,
    source = "subset",
    showlegend = FALSE,
    hoverinfo = "y",
    marker = list(size = 5)
  ) %>%
  plotly::layout(
    title = "",
    xaxis = list(title = ""),
    yaxis = list(title = "Expression level", hoverformat = ".4f"),
    dragmode = "select",
    hovermode = "compare"
  )
})

observeEvent(input$expression_by_sample_info, {
  showModal(
    modalDialog(
      expression_by_sample_info$text,
      title = expression_by_sample_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})

##--------------------------------------------------------------------------##

# box plot by cluster
output$expression_by_cluster <- plotly::renderPlotly({
  if ( length(genesToPlot()$genes_to_display_present) == 0 ) {
    expression_levels <- 0
  } else if ( length(genesToPlot()$genes_to_display_present) == 1 ) {
    expression_levels <- genesToPlot()$genes_to_display_present %>%
      sample_data()$expression[ . , ]
  } else {
    expression_levels <- genesToPlot()$genes_to_display_present %>%
      sample_data()$expression[ . , ] %>%
      colMeans() %>%
      as.vector()
  }
  data.frame(
    "cluster" = sample_data()$cells[ , "cluster" ],
    "expression" = expression_levels
  ) %>%
  plotly::plot_ly(
    x = ~cluster,
    y = ~expression,
    type = "box",
    color = ~cluster,
    colors = sample_data()$clusters$colors,
    source = "subset",
    showlegend = FALSE,
    hoverinfo = "y",
    marker = list(size = 5)
  ) %>%
  plotly::layout(
    title = "",
    xaxis = list(title = ""),
    yaxis = list(title = "Expression level", hoverformat = ".4f"),
    dragmode =  "select",
    hovermode = "compare"
  )
})

observeEvent(input$expression_by_cluster_info, {
  showModal(
    modalDialog(
      expression_by_cluster_info$text,
      title = expression_by_cluster_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})
