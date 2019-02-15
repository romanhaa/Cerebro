
##--------------------------------------------------------------------------##
## Panel: Gene set expression
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##

# reactive data
geneSets <- reactive({
  if ( sample_data()$experiment$organism == "mm" ) {
    msigdbr::msigdbr(species = "Mus musculus")
  } else if ( sample_data()$experiment$organism == "hg" ) {
    msigdbr::msigdbr(species = "Homo sapiens")
  } else {
    msigdbr::msigdbr(species = "Mus musculus")
  }
})

# reactive data
geneSetData <- reactive({
  geneSetData <- list()
  if ( is.null(input$geneSetexpression_select_geneSet) || is.na(input$geneSetexpression_select_geneSet) || input$geneSetexpression_select_geneSet == "-" ) {
    geneSetData$genes_to_display_present <- NULL
  } else {
    geneSetData$genes_to_display <- geneSets()[ which(geneSets()$gs_name == input$geneSetexpression_select_geneSet) , "gene_symbol" ]$gene_symbol
    geneSetData$genes_to_display <- unique(geneSetData$genes_to_display)
    geneSetData$genes_to_display_here <- rownames(sample_data()$expression)[ match(tolower(geneSetData$genes_to_display), tolower(rownames(sample_data()$expression))) ]
    geneSetData$genes_to_display_present <- geneSetData$genes_to_display_here[ which(!is.na(geneSetData$genes_to_display_here)) ]
    geneSetData$genes_to_display_missing <- geneSetData$genes_to_display[ which(is.na(geneSetData$genes_to_display_here)) ]
  }
  geneSetData
})

# show which genes are in the data set and which aren"t
output$geneSetexpression_genes_displayed <- renderText({
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

##--------------------------------------------------------------------------##

# UI
output$geneSetexpression_UI <- renderUI({
  tagList(
    selectInput(
      "geneSetexpression_select_geneSet", label = "Gene set:",
      choices = c("-", unique(geneSets()$gs_name)), selected = "-"
    ),
    selectInput(
      "geneSetexpression_projection_to_display", label = "Projection:",
      choices = names(sample_data()$projections)
    ),
    shinyWidgets::pickerInput(
      "geneSetexpression_samples_to_display", label = "Samples to display:",
      choices = sample_data()$samples$overview$sample,
      selected = sample_data()$samples$overview$sample,
      options = list("actions-box"=TRUE), multiple = TRUE
    ),
    shinyWidgets::pickerInput(
      "geneSetexpression_clusters_to_display", label = "Clusters to display:",
      choices = sample_data()$clusters$overview$cluster,
      selected = sample_data()$clusters$overview$cluster,
      options = list("actions-box"=TRUE), multiple = TRUE
    ),
    sliderInput("geneSetexpression_percentage_cells_to_show",
      label = "Show % of cells",
      min = 10,
      max = 100,
      step = 10,
      value = 100
    ),
    selectInput(
      "geneSetexpression_plotting_order", label = "Plotting order:",
      choices = c("Random", "Highest expression on top")
    ),
    sliderInput(
      "geneSetexpression_projection_dot_size", label = "Point size:",
      min = 0, max = 50, value = 25, step = 1
    ),
    sliderInput(
      "geneSetexpression_projection_opacity", label = "Point opacity:",
      min = 0, max = 1, value = 1, step = 0.05
    )
  )
})

output$geneSetexpression_scales <- renderUI({
  projection_to_display <- if ( is.null(input$geneSetexpression_projection_to_display) || is.na(input$geneSetexpression_projection_to_display) ) names(sample_data()$projections)[1] else input$geneSetexpression_projection_to_display
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
      "geneSetexpression_projection_scale_x_manual_range", label = "X axis",
      min = range_x_min, max = range_x_max,
      value = c(range_x_min, range_x_max)
    ),
    sliderInput(
      "geneSetexpression_projection_scale_y_manual_range", label = "Y axis",
      min = range_y_min, max = range_y_max,
      value = c(range_y_min, range_y_max)
    )
  )
})

##--------------------------------------------------------------------------##

# projection with scatterD3
output$geneSetexpression_projection <- scatterD3::renderScatterD3({
  req(input$geneSetexpression_projection_to_display)
  req(input$geneSetexpression_samples_to_display)
  req(input$geneSetexpression_clusters_to_display)
  req(input$geneSetexpression_plotting_order)
  projection_to_display <- input$geneSetexpression_projection_to_display
  samples_to_display <- input$geneSetexpression_samples_to_display
  clusters_to_display <- input$geneSetexpression_clusters_to_display
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
  # randomly remove cells
  if ( input$geneSetexpression_percentage_cells_to_show < 100 ) {
    number_of_cells_to_plot <- ceiling(
      input$geneSetexpression_percentage_cells_to_show / 100 * length(cells_to_display)
    )
    cells_to_display <- cells_to_display[ sample(1:length(cells_to_display), number_of_cells_to_plot) ]
  }
  to_plot <- cbind(
      sample_data()$projections[[ projection_to_display ]][ cells_to_display , ],
      sample_data()$cells[ cells_to_display , ]
    )

  if ( length(geneSetData()$genes_to_display_present) == 0 ) {
    to_plot$level <- 0
  } else if ( length(geneSetData()$genes_to_display_present) == 1 ) {
    to_plot$level <- geneSetData()$genes_to_display_present %>%
      sample_data()$expression[ . , cells_to_display ]
  } else {
    to_plot$level <- geneSetData()$genes_to_display_present %>%
      sample_data()$expression[ . , cells_to_display ] %>%
      colMeans() %>%
      as.vector()
  }

  if ( input$geneSetexpression_plotting_order == "Random" ) {
    to_plot <- sample(1:nrow(to_plot), nrow(to_plot)) %>% to_plot[ . , ]
  } else if ( input$geneSetexpression_plotting_order == "Highest expression on top" ) {
    to_plot <- to_plot[ order(to_plot$level, decreasing=FALSE) , ]
  }
  scatterD3::scatterD3(
    x = to_plot[ , 1 ],
    y = to_plot[ , 2 ],
    xlab = colnames(to_plot)[ 1 ],
    ylab = colnames(to_plot)[ 2 ],
    xlim = c(
        input$geneSetexpression_projection_scale_x_manual_range[1],
        input$geneSetexpression_projection_scale_x_manual_range[2]
      ),
    ylim = c(
        input$geneSetexpression_projection_scale_y_manual_range[1],
        input$geneSetexpression_projection_scale_y_manual_range[2]
      ),
    point_size = input$geneSetexpression_projection_dot_size,
    col_var = to_plot$level,
    col_lab = "Gene expression",
    col_continuous = TRUE,
    point_opacity = input$geneSetexpression_projection_opacity,
    hover_size = 4,
    hover_opacity = 1,
    transitions = FALSE,
    legend_width = 0,
    menu = FALSE,
    tooltip_text = paste0(
      "<b>Cell</b>: ", to_plot[ , "cell_barcode" ], "<br/>",
      "<b>Sample</b>: ", to_plot[ , "sample" ], "<br/>",
      "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br/>",
      "<b>Expression level</b>: ", formatC(to_plot[ , "level" ], digits = 3), "<br/>"
    )
  )
})

observeEvent(input$geneSetexpression_projection_info, {
  showModal(
    modalDialog(
      geneSetexpression_projection_info$text,
      title = geneSetexpression_projection_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})

observeEvent(input$geneSetexpression_projection_export, {
  library("ggplot2")

  projection_to_display <- input$geneSetexpression_projection_to_display
  samples_to_display <- input$geneSetexpression_samples_to_display
  clusters_to_display <- input$geneSetexpression_clusters_to_display
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
      input$geneSetexpression_projection_scale_x_manual_range[1],
      input$geneSetexpression_projection_scale_x_manual_range[2]
    )
  ylim <- c(
      input$geneSetexpression_projection_scale_y_manual_range[1],
      input$geneSetexpression_projection_scale_y_manual_range[2]
    )

  if ( length(geneSetData()$genes_to_display_present) == 0 ) {
    to_plot$level <- 0
  } else if ( length(geneSetData()$genes_to_display_present) == 1 ) {
    to_plot$level <- geneSetData()$genes_to_display_present %>%
      sample_data()$expression[ . , cells_to_display ]
  } else {
    to_plot$level <- geneSetData()$genes_to_display_present %>%
      sample_data()$expression[ . , cells_to_display ] %>%
      colMeans() %>%
      as.vector()
  }

  out_filename <- paste0(
      plot_export_path, "Cerebro_", sample_data()$experiment$experiment_name,
      "_gene_set_expression_", input$geneSetexpression_select_geneSet, "_",
      input$geneSetexpression_projection_to_display
    )

  if ( input$geneSetexpression_plotting_order == "Random" ) {
    to_plot <- sample(1:nrow(to_plot), nrow(to_plot)) %>% to_plot[ . , ]
    out_filename <- paste0(out_filename, "_random_order.pdf")
  } else if ( input$geneSetexpression_plotting_order == "Highest expression on top" ) {
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
    scale_colour_distiller(
      palette = "YlGnBu",
      direction = 1,
      name = "Average\nlog-normalised\nexpression",
      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),
      limits = c(
        min(to_plot$level)-((max(to_plot$level)-min(to_plot$level))/10),
        max(to_plot$level)
      )
    ) +
    lims(x = xlim, y = ylim) +
    theme_bw()

  pdf(NULL)
  ggsave(out_filename, p, height = 8, width = 11)

  if (file.exists(out_filename)) {
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
output$geneSetexpression_by_sample <- plotly::renderPlotly({
  if ( length(geneSetData()$genes_to_display_present) == 0 ) {
    expression_levels <- 0
  } else if ( length(geneSetData()$genes_to_display_present) == 1 ) {
    expression_levels <- geneSetData()$genes_to_display_present %>%
      sample_data()$expression[ . , ]
  } else {
    expression_levels <- geneSetData()$genes_to_display_present %>%
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
    yaxis = list(title = "Expression level", hoverformat=".2f"),
    dragmode = "select",
    hovermode = "compare"
  )
})

observeEvent(input$geneSetexpression_by_sample_info, {
  showModal(
    modalDialog(
      geneSetexpression_by_sample_info$text,
      title = geneSetexpression_by_sample_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})

##--------------------------------------------------------------------------##

# box plot by cluster
output$geneSetexpression_by_cluster <- plotly::renderPlotly({
  if ( length(geneSetData()$genes_to_display_present) == 0 ) {
    expression_levels <- 0
  } else if ( length(geneSetData()$genes_to_display_present) == 1 ) {
    expression_levels <- geneSetData()$genes_to_display_present %>%
      sample_data()$expression[ . , ]
  } else {
    expression_levels <- geneSetData()$genes_to_display_present %>%
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
    yaxis = list(title = "Expression level", hoverformat = ".2f"),
    dragmode = "select",
    hovermode = "compare"
  )
})

observeEvent(input$geneSetexpression_by_cluster_info, {
  showModal(
    modalDialog(
      geneSetexpression_by_cluster_info$text,
      title = geneSetexpression_by_cluster_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})
