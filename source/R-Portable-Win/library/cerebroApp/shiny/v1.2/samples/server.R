##----------------------------------------------------------------------------##
## Tab: Samples.
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## Samples by clusters: table + bar plot.
##----------------------------------------------------------------------------##

# UI element: buttons
output[["samples_by_cluster_UI_buttons"]] <- renderUI({
  tagList(
    shinyWidgets::materialSwitch(
      inputId = "samples_by_cluster_select_metric_for_bar_plot",
      label = "Show composition in percent [%]:",
      status = "primary",
      inline = TRUE
    ),
    shinyWidgets::materialSwitch(
      inputId = "samples_by_cluster_show_table",
      label = "Show table:",
      status = "primary",
      inline = TRUE
    )
  )
})

# UI element: rest
output[["samples_by_cluster_UI_rest"]] <- renderUI({
  tagList(
    plotly::plotlyOutput("samples_by_cluster_plot"),
    {
      if ( !is.null(input[["samples_by_cluster_show_table"]]) && input[["samples_by_cluster_show_table"]] == TRUE ) {
        DT::dataTableOutput("samples_by_cluster_table")
      }
    }
  )
})

# bar plot
output[["samples_by_cluster_plot"]] <- plotly::renderPlotly({
  # calculate table (must be merged later if user chooses to display in %)
  temp_table_original <- calculateTableAB('sample','cluster')
  # process table
  temp_table_to_plot <- temp_table_original %>%
    select(-total_cell_count) %>%
    reshape2::melt(id.vars = "sample") %>%
    rename(cluster = variable, cells = value)
  if ( input[["samples_by_cluster_select_metric_for_bar_plot"]] != TRUE ) {
    # generate bar plot with actual cell counts
    temp_table_to_plot %>%
    plotly::plot_ly(
      x = ~sample,
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
    # normalize counts to 100% and generate bar plot in percent
    temp_table_to_plot %>%
    left_join(
      .,
      temp_table_original[ , c("sample", "total_cell_count") ],
      by = "sample"
    ) %>%
    mutate(pct = cells / total_cell_count * 100) %>%
    plotly::plot_ly(
      x = ~sample,
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
        title = "Percent of cells [%]",
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
output[["samples_by_cluster_table"]] <- DT::renderDataTable({
  # generate table
  temp_table <- calculateTableAB('sample','cluster')
  if ( input[["samples_by_cluster_select_metric_for_bar_plot"]] == TRUE ) {
    # normalize counts to 100% percent
    for ( i in 3:ncol(temp_table) ) {
      temp_table[,i] <- round(temp_table[,i] / temp_table$total_cell_count * 100, digits = 1)
    }
  }
  # process table and convert to DT
  temp_table %>%
  rename(
    Sample = sample,
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
    colors = reactive_colors()$samples,
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

# UI element: buttons
output[["samples_by_cell_cycle_seurat_UI_buttons"]] <- renderUI({
  if ( "cell_cycle_seurat" %in% colnames(sample_data()$cells) ) {
    tagList(
      shinyWidgets::materialSwitch(
        inputId = "samples_by_cell_cycle_seurat_select_metric_for_bar_plot",
        label = "Show composition in percent [%]:",
        status = "primary",
        inline = TRUE
      ),
      shinyWidgets::materialSwitch(
        inputId = "samples_by_cell_cycle_seurat_show_table",
        label = "Show table:",
        status = "primary",
        inline = TRUE
      )
    )
  } else {
    textOutput("samples_by_cell_cycle_seurat_text")
  }
})

# UI element: rest
output[["samples_by_cell_cycle_seurat_UI_rest"]] <- renderUI({
  if ( "cell_cycle_seurat" %in% colnames(sample_data()$cells) ) {
    tagList(
      plotly::plotlyOutput("samples_by_cell_cycle_seurat_plot"),
      {
        if ( !is.null(input[["samples_by_cell_cycle_seurat_show_table"]]) && input[["samples_by_cell_cycle_seurat_show_table"]] == TRUE ) {
          DT::dataTableOutput("samples_by_cell_cycle_seurat_table")
        }
      }
    )
  }
})

# bar plot
output[["samples_by_cell_cycle_seurat_plot"]] <- plotly::renderPlotly({
  temp_table_original <- calculateTableAB('sample','cell_cycle_seurat')
  temp_table_to_plot <- temp_table_original %>%
    select(-total_cell_count) %>%
    reshape2::melt(id.vars = "sample") %>%
    rename(phase = variable, cells = value) %>%
    mutate(phase = factor(phase, levels = c("G1", "S", "G2M")))
  if ( input[["samples_by_cell_cycle_seurat_select_metric_for_bar_plot"]] != TRUE ) {
    temp_table_to_plot %>%
    plotly::plot_ly(
      x = ~sample,
      y = ~cells,
      type = "bar",
      color = ~phase,
      colors = cell_cycle_colorset,
      hoverinfo = "text",
      text = ~paste0("<b>", .$phase, ": </b>", formatC(.$cells, big.mark = ','))
    ) %>%
    plotly::layout(
      xaxis = list(
        title ="",
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
      temp_table_original[ , c("sample", "total_cell_count") ],
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
        title = "Percent of cells [%]",
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
output[["samples_by_cell_cycle_seurat_table"]] <- DT::renderDataTable({
  temp_table <- calculateTableAB('sample','cell_cycle_seurat')
  if ( input[["samples_by_cell_cycle_seurat_select_metric_for_bar_plot"]] == TRUE ) {
    for ( i in 3:ncol(temp_table) ) {
      temp_table[,i] <- round(temp_table[,i] / temp_table$total_cell_count * 100, digits = 1)
    }
  }
  temp_table %>%
  rename(
    Sample = sample,
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

# UI element: buttons
output[["samples_by_cell_cycle_cyclone_UI_buttons"]] <- renderUI({
  if ( "cell_cycle_cyclone" %in% colnames(sample_data()$cells) ) {
    tagList(
      shinyWidgets::materialSwitch(
        inputId = "samples_by_cell_cycle_cyclone_select_metric_for_bar_plot",
        label = "Show composition in percent [%]:",
        status = "primary",
        inline = TRUE
      ),
      shinyWidgets::materialSwitch(
        inputId = "samples_by_cell_cycle_cyclone_show_table",
        label = "Show table:",
        status = "primary",
        inline = TRUE
      )
    )
  } else {
    textOutput("samples_by_cell_cycle_cyclone_text")
  }
})

# UI element: rest
output[["samples_by_cell_cycle_cyclone_UI_rest"]] <- renderUI({
  if ( "cell_cycle_cyclone" %in% colnames(sample_data()$cells) ) {
    tagList(
      plotly::plotlyOutput("samples_by_cell_cycle_cyclone_plot"),
      {
        if ( !is.null(input[["samples_by_cell_cycle_cyclone_show_table"]]) && input[["samples_by_cell_cycle_cyclone_show_table"]] == TRUE ) {
          DT::dataTableOutput("samples_by_cell_cycle_cyclone_table")
        }
      }
    )
  }
})

# bar plot
output[["samples_by_cell_cycle_cyclone_plot"]] <- plotly::renderPlotly({
  temp_table_original <- calculateTableAB('sample','cell_cycle_cyclone')
  temp_table_to_plot <- temp_table_original %>%
    select(-total_cell_count) %>%
    reshape2::melt(id.vars = "sample") %>%
    rename(phase = variable, cells = value) %>%
    mutate(phase = factor(phase, levels = c("G1", "S", "G2M", "-")))
  if ( input[["samples_by_cell_cycle_cyclone_select_metric_for_bar_plot"]] != TRUE ) {
    temp_table_to_plot %>%
    plotly::plot_ly(
      x = ~sample,
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
      temp_table_original[ , c("sample", "total_cell_count") ],
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
        title = "Percent of cells [%]",
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
output[["samples_by_cell_cycle_cyclone_table"]] <- DT::renderDataTable({
  temp_table <- calculateTableAB('sample','cell_cycle_cyclone')
  if ( input[["samples_by_cell_cycle_cyclone_select_metric_for_bar_plot"]] == TRUE ) {
    for ( i in 3:ncol(temp_table) ) {
      temp_table[,i] <- round(temp_table[,i] / temp_table$total_cell_count * 100, digits = 1)
    }
  }
  temp_table %>%
  rename(
    Sample = sample,
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
