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
  colors_tree <- colors[1:length(tree$tip.label)]
  ggplot(tree, aes(x, y)) +
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

# UI element
output[["clusters_by_sample_UI"]] <- renderUI({
  if ( !is.null(sample_data()$clusters$by_sample) ) {
    tagList(
      DT::dataTableOutput("clusters_by_sample_table"),
      plotly::plotlyOutput("clusters_by_sample_plot")
    )
  } else {
    textOutput("clusters_by_sample_text")
  }
})

# table
output[["clusters_by_sample_table"]] <- DT::renderDataTable({
  sample_data()$clusters$by_sample %>%
  rename(
    Cluster = cluster,
    "# of cells" = total_cell_count
  ) %>%
  DT::datatable(
    filter = "none",
    selection = "multiple",
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

# bar plot
output[["clusters_by_sample_plot"]] <- plotly::renderPlotly({
  sample_data()$clusters$by_sample %>%
  select(-total_cell_count) %>%
  reshape2::melt(id.vars = "cluster") %>%
  rename(sample = variable, cells = value) %>%
  left_join(
    .,
    sample_data()$clusters$by_sample[ , c("cluster", "total_cell_count") ],
    by = "cluster"
  ) %>%
  mutate(pct = cells / total_cell_count * 100) %>%
  plotly::plot_ly(
    x = ~cluster,
    y = ~pct,
    type = "bar",
    color = ~sample,
    colors = sample_data()$samples$colors,
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
})

# alternative text
output[["clusters_by_sample_text"]] <- renderText({
    "Only 1 sample in this data set."
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

# UI element
output[["clusters_by_cell_cycle_seurat_UI"]] <- renderUI({
  if ( !is.null(sample_data()$clusters$by_cell_cycle_seurat) ) {
    plotly::plotlyOutput("clusters_by_cell_cycle_seurat_plot")
  } else {
    textOutput("clusters_by_cell_cycle_seurat_text")
  }
})

# bar plot
output[["clusters_by_cell_cycle_seurat_plot"]] <- plotly::renderPlotly({
  sample_data()$clusters$by_cell_cycle_seurat %>%
  select(-total_cell_count) %>%
  reshape2::melt(id.vars = "cluster") %>%
  rename(phase = variable, cells = value) %>%
  mutate(
    phase = factor(phase, levels = c("G1", "S", "G2M")),
  ) %>%
  left_join(
    .,
    sample_data()$clusters$by_cell_cycle_seurat[ , c("cluster", "total_cell_count") ],
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

# UI element
output[["clusters_by_cell_cycle_cyclone_UI"]] <- renderUI(
  if ( !is.null(sample_data()$cells$cell_cycle_cyclone) ) {
    plotly::plotlyOutput("clusters_by_cell_cycle_cyclone_plot")
  } else {
    textOutput("clusters_by_cell_cycle_cyclone_text")
  }
)

# bar plot
output[["clusters_by_cell_cycle_cyclone_plot"]] <- plotly::renderPlotly({
  sample_data()$clusters$by_cell_cycle_cyclone %>%
  select(-total_cell_count) %>%
  reshape2::melt(id.vars = "cluster") %>%
  rename(phase = variable, cells = value) %>%
  mutate(
    phase = factor(phase, levels = c("G1", "S", "G2M", "-")),
  ) %>%
  left_join(
    .,
    sample_data()$clusters$by_cell_cycle_cyclone[ , c("cluster", "total_cell_count") ],
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
