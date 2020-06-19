##----------------------------------------------------------------------------##
## Tab: Clusters.
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## Cluster tree.
##----------------------------------------------------------------------------##

# UI element
output[["clusters_tree_UI"]] <- renderUI({
  if ( !is.null(sample_data()$clusters$tree) ) {
    plotOutput("clusters_tree_plot")
  } else {
    textOutput("clusters_tree_text")
  }
})

# plot
output[["clusters_tree_plot"]] <- renderPlot({
  tree <- sample_data()$clusters$tree
  tree$tip.label <- paste0("Cluster ", tree$tip.label)
  colors_tree <- reactive_colors()$clusters
  ggtree::ggtree(tree, aes(x, y)) +
    scale_y_reverse() +
    ggtree::geom_tree() +
    ggtree::theme_tree() +
    ggtree::geom_tiplab(size = 5, hjust = -0.2) +
    ggtree::geom_tippoint(color = colors_tree, shape = 16, size = 6) +
    coord_cartesian(clip = 'off') +
    theme(plot.margin = unit(c(0,2.5,0,0), 'cm'))
})

# alternative text
output[["clusters_tree_text"]] <- renderText({ "Data not available." })

# info box
observeEvent(input[["clusters_tree_info"]], {
  showModal(
    modalDialog(
      clusters_tree_info[["text"]],
      title = clusters_tree_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Clusters by samples.
##----------------------------------------------------------------------------##

# UI element: buttons
output[["clusters_by_sample_UI_buttons"]] <- renderUI({
  tagList(
    shinyWidgets::materialSwitch(
      inputId = "clusters_by_sample_select_metric_for_bar_plot",
      label = "Show composition in percent [%]:",
      status = "primary",
      inline = TRUE
    ),
    shinyWidgets::materialSwitch(
      inputId = "clusters_by_sample_show_table",
      label = "Show table:",
      status = "primary",
      inline = TRUE
    )
  )
})

# UI element: rest
output[["clusters_by_sample_UI_rest"]] <- renderUI({
  tagList(
    plotly::plotlyOutput("clusters_by_sample_plot"),
    {
      if ( !is.null(input[["clusters_by_sample_show_table"]]) && input[["clusters_by_sample_show_table"]] == TRUE ) {
        DT::dataTableOutput("clusters_by_sample_table")
      }
    }
  )
})

# bar plot
output[["clusters_by_sample_plot"]] <- plotly::renderPlotly({
  # calculate table (must be merged later if user chooses to display in %)
  temp_table_original <- calculateTableAB('cluster','sample')
  # process table
  temp_table_to_plot <- temp_table_original %>%
    select(-total_cell_count) %>%
    reshape2::melt(id.vars = "cluster") %>%
    rename(sample = variable, cells = value)
  if ( input[['clusters_by_sample_select_metric_for_bar_plot']] != TRUE ) {
    # generate bar plot with actual cell counts
    temp_table_to_plot %>%
    plotly::plot_ly(
      x = ~cluster,
      y = ~cells,
      type = "bar",
      color = ~sample,
      colors = reactive_colors()$samples,
      hoverinfo = "text",
      text = ~paste0("<b>", .$sample, ": </b>", formatC(.$cells, big.mark = ','))
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
    # normalize counts to 100% and generate bar plot in percent
    temp_table_to_plot %>%
    left_join(
      .,
      temp_table_original[ , c("cluster", "total_cell_count") ],
      by = "cluster"
    ) %>%
    mutate(pct = cells / total_cell_count * 100) %>%
    plotly::plot_ly(
      x = ~cluster,
      y = ~pct,
      type = "bar",
      color = ~sample,
      colors = reactive_colors()$samples,
      hoverinfo = "text",
      text = ~paste0("<b>", .$sample, ": </b>", format(round(.$pct, 1), nsmall = 1), "%")
    ) %>%
    plotly::layout(
      xaxis = list(
        title = "",
        mirror = TRUE,
        showline = TRUE
      ),
      yaxis = list(
        title = "Percentage [%]",
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
output[["clusters_by_sample_table"]] <- DT::renderDataTable({
  # generate table
  temp_table <- calculateTableAB('cluster','sample')
  if ( input[["clusters_by_sample_select_metric_for_bar_plot"]] == TRUE ) {
    # normalize counts to 100% percent
    for ( i in 3:ncol(temp_table) ) {
      temp_table[,i] <- round(temp_table[,i] / temp_table$total_cell_count * 100, digits = 1)
    }
  }
  # process table and convert to DT
  temp_table %>%
  rename(
    Cluster = cluster,
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
      lengthMenu = c(20, 30, 50, 100),
      pageLength = 20
    )
  )
})

# info box
observeEvent(input[["clusters_by_sample_info"]], {
  showModal(
    modalDialog(
      clusters_by_sample_info[["text"]],
      title = clusters_by_sample_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## nUMI
##----------------------------------------------------------------------------##

# UI element
output[["clusters_nUMI_UI"]] <- renderUI({
  if ( "nUMI" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("clusters_nUMI_plot")
  } else {
    textOutput("clusters_nUMI_text")
  }
})

# box plot
output[["clusters_nUMI_plot"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~cluster,
    y = ~nUMI,
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
      title = "Number of transcripts",
      hoverformat = ".0f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# alternative text
output[["clusters_nUMI_text"]] <- renderText({
  "Data not available."
})

# info box
observeEvent(input[["clusters_nUMI_info"]], {
  showModal(
    modalDialog(
      clusters_nUMI_info[["text"]],
      title = clusters_nUMI_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## nGene
##----------------------------------------------------------------------------##

# UI element
output[["clusters_nGene_UI"]] <- renderUI({
  if ( "nGene" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("clusters_nGene_plot")
  } else {
    textOutput("clusters_nGene_text")
  }
})

# box plot
output[["clusters_nGene_plot"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~cluster,
    y = ~nGene,
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
      title ="",
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

# alternative text
output[["clusters_nGene_text"]] <- renderText({
  "Data not available."
})

# info box
observeEvent(input[["clusters_nGene_info"]], {
  showModal(
    modalDialog(
      clusters_nGene_info[["text"]],
      title = clusters_nGene_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## percent_mt
##----------------------------------------------------------------------------##

# UI element
output[["clusters_percent_mt_UI"]] <- renderUI({
  if ( "percent_mt" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("clusters_percent_mt_plot")
  } else {
    textOutput("clusters_percent_mt_text")
  }
})

# box plot
output[["clusters_percent_mt_plot"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~cluster,
    y = ~percent_mt*100,
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
      title = "Percentage of transcripts [%]",
      range = c(0, 100),
      hoverformat = ".1f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# alternative text
output[["clusters_percent_mt_text"]] <- renderText({
  "Data not available."
})

# info box
observeEvent(input[["clusters_percent_mt_info"]], {
  showModal(
    modalDialog(
      clusters_percent_mt_info[["text"]],
      title = clusters_percent_mt_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## percent_ribo
##----------------------------------------------------------------------------##

# UI element
output[["clusters_percent_ribo_UI"]] <- renderUI({
  if ( "percent_ribo" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("clusters_percent_ribo_plot")
  } else {
    textOutput("clusters_percent_ribo_text")
  }
})

# alternative text
output[["clusters_percent_ribo_text"]] <- renderText({
  "Data not available."
})

# box plot
output[["clusters_percent_ribo_plot"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~cluster,
    y = ~percent_ribo*100,
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
      title = "Percentage of transcripts [%]",
      range = c(0, 100),
      hoverformat = ".2f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# info box
observeEvent(input[["clusters_percent_ribo_info"]], {
  showModal(
    modalDialog(
      clusters_percent_ribo_info[["text"]],
      title = clusters_percent_ribo_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## cell cycle: Seurat
##----------------------------------------------------------------------------##

# UI element: buttons
output[["clusters_by_cell_cycle_seurat_UI_buttons"]] <- renderUI({
  if ( "cell_cycle_seurat" %in% colnames(sample_data()$cells) ) {
    tagList(
      shinyWidgets::materialSwitch(
        inputId = "clusters_by_cell_cycle_seurat_select_metric_for_bar_plot",
        label = "Show composition in percent [%]:",
        status = "primary",
        inline = TRUE
      ),
      shinyWidgets::materialSwitch(
        inputId = "clusters_by_cell_cycle_seurat_show_table",
        label = "Show table:",
        status = "primary",
        inline = TRUE
      )
    )
  } else {
    textOutput("clusters_by_cell_cycle_seurat_text")
  }
})

# UI element: rest
output[["clusters_by_cell_cycle_seurat_UI_rest"]] <- renderUI({
  if ( "cell_cycle_seurat" %in% colnames(sample_data()$cells) ) {
    tagList(
      plotly::plotlyOutput("clusters_by_cell_cycle_seurat_plot"),
      {
        if ( !is.null(input[["clusters_by_cell_cycle_seurat_show_table"]]) && input[["clusters_by_cell_cycle_seurat_show_table"]] == TRUE ) {
          DT::dataTableOutput("clusters_by_cell_cycle_seurat_table")
        }
      }
    )
  }
})

# bar plot
output[["clusters_by_cell_cycle_seurat_plot"]] <- plotly::renderPlotly({
  temp_table_original <- calculateTableAB('cluster','cell_cycle_seurat')
  temp_table_to_plot <- temp_table_original %>%
    select(-total_cell_count) %>%
    reshape2::melt(id.vars = "cluster") %>%
    rename(phase = variable, cells = value) %>%
    mutate(phase = factor(phase, levels = c("G1", "S", "G2M")))
  if ( input[['clusters_by_cell_cycle_seurat_select_metric_for_bar_plot']] != TRUE ) {
    temp_table_to_plot %>%
    plotly::plot_ly(
      x = ~cluster,
      y = ~cells,
      type = "bar",
      color = ~phase,
      colors = cell_cycle_colorset,
      hoverinfo = "text",
      text = ~paste0("<b>", .$phase, ": </b>", formatC(.$cells, big.mark = ','))
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
    temp_table_to_plot %>%
    left_join(
      .,
      temp_table_original[ , c("cluster", "total_cell_count") ],
      by = "cluster"
    ) %>%
    mutate(pct = cells / total_cell_count * 100) %>%
    plotly::plot_ly(
      x = ~cluster,
      y = ~pct,
      type = "bar",
      color = ~phase,
      colors = cell_cycle_colorset,
      hoverinfo = "text",
      text = ~paste0("<b>", .$phase, ": </b>", format(round(.$pct, 1), nsmall = 1), "%")
    ) %>%
    plotly::layout(
      xaxis = list(
        title = "",
        mirror = TRUE,
        showline = TRUE
      ),
      yaxis = list(
        title = "Percentage [%]",
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
output[["clusters_by_cell_cycle_seurat_table"]] <- DT::renderDataTable({
  temp_table <- calculateTableAB('cluster','cell_cycle_seurat')
  if ( input[["clusters_by_cell_cycle_seurat_select_metric_for_bar_plot"]] == TRUE ) {
    for ( i in 3:ncol(temp_table) ) {
      temp_table[,i] <- round(temp_table[,i] / temp_table$total_cell_count * 100, digits = 1)
    }
  }
  temp_table %>%
  rename(
    Cluster = cluster,
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
output[["clusters_by_cell_cycle_seurat_text"]] <- renderText({
  "Data not available."
})

# info box
observeEvent(input[["clusters_by_cell_cycle_seurat_info"]], {
  showModal(
    modalDialog(
      clusters_by_cell_cycle_seurat_info[["text"]],
      title = clusters_by_cell_cycle_seurat_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## cell cycle: Cyclone
##----------------------------------------------------------------------------##

# UI element: buttons
output[["clusters_by_cell_cycle_cyclone_UI_buttons"]] <- renderUI({
  if ( "cell_cycle_cyclone" %in% colnames(sample_data()$cells) ) {
    tagList(
      shinyWidgets::materialSwitch(
        inputId = "clusters_by_cell_cycle_cyclone_select_metric_for_bar_plot",
        label = "Show composition in percent [%]:",
        status = "primary",
        inline = TRUE
      ),
      shinyWidgets::materialSwitch(
        inputId = "clusters_by_cell_cycle_cyclone_show_table",
        label = "Show table:",
        status = "primary",
        inline = TRUE
      )
    )
  } else {
    textOutput("clusters_by_cell_cycle_cyclone_text")
  }
})

# UI element: rest
output[["clusters_by_cell_cycle_cyclone_UI_rest"]] <- renderUI({
  if ( "cell_cycle_cyclone" %in% colnames(sample_data()$cells) ) {
    tagList(
      plotly::plotlyOutput("clusters_by_cell_cycle_cyclone_plot"),
      {
        if ( !is.null(input[["clusters_by_cell_cycle_cyclone_show_table"]]) && input[["clusters_by_cell_cycle_cyclone_show_table"]] == TRUE ) {
          DT::dataTableOutput("clusters_by_cell_cycle_cyclone_table")
        }
      }
    )
  }
})

# bar plot
output[["clusters_by_cell_cycle_cyclone_plot"]] <- plotly::renderPlotly({
  temp_table_original <- calculateTableAB('cluster','cell_cycle_cyclone')
  temp_table_to_plot <- temp_table_original %>%
    select(-total_cell_count) %>%
    reshape2::melt(id.vars = "cluster") %>%
    rename(phase = variable, cells = value) %>%
    mutate(phase = factor(phase, levels = c("G1", "S", "G2M", "-")))
  if ( input[['clusters_by_cell_cycle_cyclone_select_metric_for_bar_plot']] != TRUE ) {
    temp_table_to_plot %>%
    plotly::plot_ly(
      x = ~cluster,
      y = ~cells,
      type = "bar",
      color = ~phase,
      colors = cell_cycle_colorset,
      hoverinfo = "text",
      text = ~paste0("<b>", .$phase, ": </b>", formatC(.$cells, big.mark = ','))
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
    temp_table_to_plot %>%
    left_join(
      .,
      temp_table_original[ , c("cluster", "total_cell_count") ],
      by = "cluster"
    ) %>%
    mutate(pct = cells / total_cell_count * 100) %>%
    plotly::plot_ly(
      x = ~cluster,
      y = ~pct,
      type = "bar",
      color = ~phase,
      colors = cell_cycle_colorset,
      hoverinfo = "text",
      text = ~paste0("<b>", .$phase, ": </b>", format(round(.$pct, 1), nsmall = 1), "%")
    ) %>%
    plotly::layout(
      xaxis = list(
        title = "",
        mirror = TRUE,
        showline = TRUE
      ),
      yaxis = list(
        title = "Percentage [%]",
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
output[["clusters_by_cell_cycle_cyclone_table"]] <- DT::renderDataTable({
  temp_table <- calculateTableAB('cluster','cell_cycle_cyclone')
  if ( input[["clusters_by_cell_cycle_cyclone_select_metric_for_bar_plot"]] == TRUE ) {
    for ( i in 3:ncol(temp_table) ) {
      temp_table[,i] <- round(temp_table[,i] / temp_table$total_cell_count * 100, digits = 1)
    }
  }
  temp_table %>%
  rename(
    Cluster = cluster,
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
output[["clusters_by_cell_cycle_cyclone_text"]] <- renderText({
  "Data not available."
})

# info box
observeEvent(input[["clusters_by_cell_cycle_cyclone_info"]], {
  showModal(
    modalDialog(
      clusters_by_cell_cycle_cyclone_info[["text"]],
      title = clusters_by_cell_cycle_cyclone_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})
