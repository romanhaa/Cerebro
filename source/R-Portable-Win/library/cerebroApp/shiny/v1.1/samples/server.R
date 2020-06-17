##----------------------------------------------------------------------------##
## Tab: Samples.
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## Samples by clusters: table + bar plot.
##----------------------------------------------------------------------------##

# UI element
output[["samples_by_cluster_UI"]] <- renderUI({
  if ( !is.null(sample_data()$samples$by_cluster) ) {
    tagList(
      DT::dataTableOutput("samples_by_cluster_table"),
      plotly::plotlyOutput("samples_by_cluster_plot")
    )
  } else {
    textOutput("samples_by_cluster_text")
  }
})

# table
output[["samples_by_cluster_table"]] <- DT::renderDataTable({
  sample_data()$samples$by_cluster %>%
  rename(
    Sample = sample,
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
      lengthMenu = c(15, 30, 50, 100),
      pageLength = 15
    )
  )
})

# bar plot
output[["samples_by_cluster_plot"]] <- plotly::renderPlotly({
  sample_data()$samples$by_cluster %>%
  select(-total_cell_count) %>%
  reshape2::melt(id.vars = "sample") %>%
  rename(cluster = variable, cells = value) %>%
  left_join(
    .,
    sample_data()$samples$by_cluster[ , c("sample", "total_cell_count") ],
    by = "sample"
  ) %>%
  mutate(pct = cells / total_cell_count * 100) %>%
  plotly::plot_ly(
    x = ~sample,
    y = ~pct,
    type = "bar",
    color = ~cluster,
    colors = sample_data()$clusters$colors,
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
      title = "Percentage [%]",
      range = c(0,100),
      hoverformat = ".2f",
      mirror = TRUE,
      zeroline = FALSE,
      showline = TRUE
    ),
    barmode = "stack",
    hovermode = "closest"
  )
})

# alternative text
output[["samples_by_cluster_text"]] <- renderText({
    "Only 1 cluster in this data set."
  })

# info button
observeEvent(input[["samples_by_cluster_info"]], {
  showModal(
    modalDialog(
      samples_by_cluster_info[["text"]],
      title = samples_by_cluster_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## nUMI
##----------------------------------------------------------------------------##

# UI element
output[["samples_nUMI_UI"]] <- renderUI({
  if ( "nUMI" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("samples_nUMI_plot")
  } else {
    textOutput("samples_nUMI_text")
  }
})

# box plot
output[["samples_nUMI_plot"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~sample,
    y = ~nUMI,
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
output[["samples_nUMI_text"]] <- renderText({
    "Column with number of transcript per cell not available."
  })

# info button
observeEvent(input[["samples_nUMI_info"]], {
  showModal(
    modalDialog(
      samples_nUMI_info[["text"]],
      title = samples_nUMI_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## nGene
##----------------------------------------------------------------------------##

# UI element
output[["samples_nGene_UI"]] <- renderUI({
  if ( "nGene" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("samples_nGene_plot")
  } else {
    textOutput("samples_nGene_text")
  }
})

# box plot
output[["samples_nGene_plot"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~sample,
    y = ~nGene,
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
output[["samples_nGene_text"]] <- renderText({
    "Column with number of expressed genes per cell not available."
  })

# info button
observeEvent(input[["samples_nGene_info"]], {
  showModal(
    modalDialog(
      samples_nGene_info[["text"]],
      title = samples_nGene_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## percent_mt
##----------------------------------------------------------------------------##

# UI element
output[["samples_percent_mt_UI"]] <- renderUI({
  if ( "percent_mt" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("samples_percent_mt_plot")
  } else {
    textOutput("samples_percent_mt_text")
  }
})

# box plot
output[["samples_percent_mt_plot"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~sample,
    y = ~percent_mt*100,
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
      title = "Percentage of transcripts [%]",
      range = c(0,100),
      hoverformat = ".1f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# alternative text
output[["samples_percent_mt_text"]] <- renderText({
    "Column with percentage of mitochondrial expression not available."
  })

# info button
observeEvent(input[["samples_percent_mt_info"]], {
  showModal(
    modalDialog(
      samples_percent_mt_info[["text"]],
      title = samples_percent_mt_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## percent_ribo
##----------------------------------------------------------------------------##

# UI element
output[["samples_percent_ribo_UI"]] <- renderUI({
  if ( "percent_ribo" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("samples_percent_ribo_plot")
  } else {
    textOutput("samples_percent_ribo_text")
  }
})

# box plot
output[["samples_percent_ribo_plot"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~sample,
    y = ~percent_ribo*100,
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
    marker = list(size = 5)
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
      range = c(0,100),
      hoverformat = ".1f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# alternative text
output[["samples_percent_ribo_text"]] <- renderText({
    "Column with percentage of ribosomal expression not available."
  })

# info button
observeEvent(input[["samples_percent_ribo_info"]], {
  showModal(
    modalDialog(
      samples_percent_ribo_info[["text"]],
      title = samples_percent_ribo_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## cell cycle: Seurat
##----------------------------------------------------------------------------##

# UI element
output[["samples_by_cell_cycle_seurat_UI"]] <- renderUI({
  if ( !is.null(sample_data()$samples$by_cell_cycle_seurat) ) {
    plotly::plotlyOutput("samples_by_cell_cycle_seurat_plot")
  } else {
    textOutput("samples_by_cell_cycle_seurat_text")
  }
})

# bar plot
output[["samples_by_cell_cycle_seurat_plot"]] <- plotly::renderPlotly({
  sample_data()$samples$by_cell_cycle_seurat %>%
  select(-total_cell_count) %>%
  reshape2::melt(id.vars = "sample") %>%
  rename(phase = variable, cells = value) %>%
  mutate(
    phase = factor(phase, levels = c("G1", "S", "G2M")),
  ) %>%
  left_join(
    .,
    sample_data()$samples$by_cell_cycle_seurat[ , c("sample", "total_cell_count") ],
    by = "sample"
  ) %>%
  mutate(pct = cells / total_cell_count * 100) %>%
  plotly::plot_ly(
    x = ~sample,
    y = ~pct,
    type = "bar",
    color = ~phase,
    colors = cell_cycle_colorset,
    hoverinfo = "text",
    text = ~paste0("<b>", .$phase, ": </b>", format(round(.$pct, 1), nsmall = 1), "%")
  ) %>%
  plotly::layout(
    xaxis = list(
      title ="",
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
output[["samples_by_cell_cycle_seurat_text"]] <- renderText({
    "Data not available."
  })

# info button
observeEvent(input[["samples_by_cell_cycle_seurat_info"]], {
  showModal(
    modalDialog(
      samples_by_cell_cycle_seurat_info[["text"]],
      title = samples_by_cell_cycle_seurat_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## cell cycle: Cyclone
##----------------------------------------------------------------------------##

# UI element
output[["samples_by_cell_cycle_cyclone_UI"]] <- renderUI({
  if ( !is.null(sample_data()$samples$by_cell_cycle_cyclone) ) {
    plotly::plotlyOutput("samples_by_cell_cycle_cyclone_plot")
  } else {
    textOutput("samples_by_cell_cycle_cyclone_text")
  }
})

# bar plot
output[["samples_by_cell_cycle_cyclone_plot"]] <- plotly::renderPlotly({
  sample_data()$samples$by_cell_cycle_cyclone %>%
  select(-total_cell_count) %>%
  reshape2::melt(id.vars = "sample") %>%
  rename(phase = variable, cells = value) %>%
  mutate(
    phase = factor(phase, levels = c("G1", "S", "G2M", "-")),
  ) %>%
  left_join(
    .,
    sample_data()$samples$by_cell_cycle_cyclone[ , c("sample", "total_cell_count") ],
    by = "sample"
  ) %>%
  mutate(pct = cells / total_cell_count * 100) %>%
  plotly::plot_ly(
    x = ~sample,
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
output[["samples_by_cell_cycle_cyclone_text"]] <- renderText({
    "Data not available."
  })

# info button
observeEvent(input[["samples_by_cell_cycle_cyclone_info"]], {
  showModal(
    modalDialog(
      samples_by_cell_cycle_cyclone_info[["text"]],
      title = samples_by_cell_cycle_cyclone_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})
