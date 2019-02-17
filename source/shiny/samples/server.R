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
  mutate(pct = cells / total_cell_count) %>%
  plotly::plot_ly(
    x = ~sample,
    y = ~pct*100,
    type = "bar",
    color = ~cluster,
    colors = sample_data()$clusters$colors,
    text = ~pct*100,
    hoverinfo = "name+y"
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
output[["samples_box_nUMI_UI"]] <- renderUI({
  if ( "nUMI" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("samples_box_nUMI_plot")
  } else {
    textOutput("samples_box_nUMI_text")
  }
})

# box plot
output[["samples_box_nUMI_plot"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~sample,
    y = ~nUMI,
    type = "box",
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
      title = "Number of UMIs",
      type = "log",
      hoverformat = ".2f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# alternative text
output[["samples_box_nUMI_text"]] <- renderText({
    "Column with number of transcript per cell not available."
  })

# info button
observeEvent(input[["samples_box_nUMI_info"]], {
  showModal(
    modalDialog(
      samples_box_nUMI_info[["text"]],
      title = samples_box_nUMI_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## nGene
##----------------------------------------------------------------------------##

# UI element
output[["samples_box_nGene_UI"]] <- renderUI({
  if ( "nGene" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("samples_box_nGene_plot")
  } else {
    textOutput("samples_box_nGene_text")
  }
})

# box plot
output[["samples_box_nGene_plot"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~sample,
    y = ~nGene,
    type = "box",
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
      type = "log",
      hoverformat = ".2f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# alternative text
output[["samples_box_nGene_text"]] <- renderText({
    "Column with number of expressed genes per cell not available."
  })

# info button
observeEvent(input[["samples_box_nGene_info"]], {
  showModal(
    modalDialog(
      samples_box_nGene_info[["text"]],
      title = samples_box_nGene_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## percent_mt
##----------------------------------------------------------------------------##

# UI element
output[["samples_box_percent_mt_UI"]] <- renderUI({
  if ( "percent_mt" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("samples_box_percent_mt_plot")
  } else {
    textOutput("samples_box_percent_mt_text")
  }
})

# box plot
output[["samples_box_percent_mt_plot"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~sample,
    y = ~percent_mt*100,
    type = "box",
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
      title = "Percentage of mitochondrial gene expression",
      range = c(0,100),
      hoverformat = ".2f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# alternative text
output[["samples_box_percent_mt_text"]] <- renderText({
    "Column with percentage of mitochondrial expression not available."
  })

# info button
observeEvent(input[["samples_box_percent_mt_info"]], {
  showModal(
    modalDialog(
      samples_box_percent_mt_info[["text"]],
      title = samples_box_percent_mt_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## percent_ribo
##----------------------------------------------------------------------------##

# UI element
output[["samples_box_percent_ribo_UI"]] <- renderUI({
  if ( "percent_ribo" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("samples_box_percent_ribo_plot")
  } else {
    textOutput("samples_box_percent_ribo_text")
  }
})

# box plot
output[["samples_box_percent_ribo_plot"]] <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~sample,
    y = ~percent_ribo*100,
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
    xaxis = list(
      title = "",
      mirror = TRUE,
      showline = TRUE
    ),
    yaxis = list(
      title = "Percentage of ribosomal gene expression",
      range = c(0,100),
      hoverformat = ".2f",
      mirror = TRUE,
      showline = TRUE
    ),
    dragmode = "select",
    hovermode = "compare"
  )
})

# alternative text
output[["samples_box_percent_ribo_text"]] <- renderText({
    "Column with percentage of ribosomal expression not available."
  })

# info button
observeEvent(input[["samples_box_percent_ribo_info"]], {
  showModal(
    modalDialog(
      samples_box_percent_ribo_info[["text"]],
      title = samples_box_percent_ribo_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## cell cycle: Seurat
##----------------------------------------------------------------------------##

# UI element
output[["samples_by_cell_cycle_Regev_UI"]] <- renderUI({
  if ( !is.null(sample_data()$samples$by_cell_cycle_Regev) ) {
    plotly::plotlyOutput("samples_by_cell_cycle_Regev_plot")
  } else {
    textOutput("samples_by_cell_cycle_Regev_text")
  }
})

# bar plot
output[["samples_by_cell_cycle_Regev_plot"]] <- plotly::renderPlotly({
  sample_data()$samples$by_cell_cycle_Regev %>%
  select(-total_cell_count) %>%
  reshape2::melt(id.vars = "sample") %>%
  rename(phase = variable, cells = value) %>%
  mutate(
    phase = factor(phase, levels = c("G1", "S", "G2M")),
  ) %>%
  left_join(
    .,
    sample_data()$samples$by_cell_cycle_Regev[ , c("sample", "total_cell_count") ],
    by = "sample"
  ) %>%
  mutate(pct = cells / total_cell_count) %>%
  plotly::plot_ly(
    x = ~sample,
    y = ~pct*100,
    type = "bar",
    color = ~phase,
    colors = cell_cycle_colorset,
    text = ~pct*100,
    hoverinfo = "name+y"
  ) %>%
  plotly::layout(
    xaxis = list(
      title ="",
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
})

# alternative text
output[["samples_by_cell_cycle_Regev_text"]] <- renderText({
    "Data not available."
  })

# info button
observeEvent(input[["samples_by_cell_cycle_Regev_info"]], {
  showModal(
    modalDialog(
      samples_by_cell_cycle_Regev_info[["text"]],
      title = samples_by_cell_cycle_Regev_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## cell cycle: Cyclone
##----------------------------------------------------------------------------##

# UI element
output[["samples_by_cell_cycle_Cyclone_UI"]] <- renderUI({
  if ( !is.null(sample_data()$samples$by_cell_cycle_Cyclone) ) {
    plotly::plotlyOutput("samples_by_cell_cycle_Cyclone_plot")
  } else {
    textOutput("samples_by_cell_cycle_Cyclone_text")
  }
})

# bar plot
output[["samples_by_cell_cycle_Cyclone_plot"]] <- plotly::renderPlotly({
  sample_data()$samples$by_cell_cycle_Cyclone %>%
  select(-total_cell_count) %>%
  reshape2::melt(id.vars = "sample") %>%
  rename(phase = variable, cells = value) %>%
  mutate(
    phase = factor(phase, levels = c("G1", "S", "G2M", "-")),
  ) %>%
  left_join(
    .,
    sample_data()$samples$by_cell_cycle_Cyclone[ , c("sample", "total_cell_count") ],
    by = "sample"
  ) %>%
  mutate(pct = cells / total_cell_count) %>%
  plotly::plot_ly(
    x = ~sample,
    y = ~pct*100,
    type = "bar",
    color = ~phase,
    colors = cell_cycle_colorset,
    text = ~pct*100,
    hoverinfo = "name+y"
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
})

# alternative text
output[["samples_by_cell_cycle_Cyclone_text"]] <- renderText({
    "Data not available."
  })

# info button
observeEvent(input[["samples_by_cell_cycle_Cyclone_info"]], {
  showModal(
    modalDialog(
      samples_by_cell_cycle_Cyclone_info[["text"]],
      title = samples_by_cell_cycle_Cyclone_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
}) 
