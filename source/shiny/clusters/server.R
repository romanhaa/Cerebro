##--------------------------------------------------------------------------##
## Tab: Clusters.
##--------------------------------------------------------------------------##

# UI element for cluster tree
output$clusters_tree_UI <- renderUI({
  if ( !is.null(sample_data()$clusters$tree) ) {
    plotOutput("clusters_tree_plot")
  } else {
    textOutput("clusters_tree_text")
  }
})

output$clusters_tree_plot <- renderPlot({
  library("ggtree")
  tree <- sample_data()$clusters$tree
  tree$tip.label <- paste0("Cluster ", tree$tip.label)
  colors_tree <- colors[1:length(tree$tip.label)]
  ggplot(tree, aes(x, y)) + 
    ggplot2::scale_y_reverse() +
    xlim(0, max(tree$edge.length * 1.1)) +
    geom_tree() +
    theme_tree() +
    geom_tiplab(size = 5, hjust = -0.2) +
    geom_tippoint(color = colors_tree, shape = 16, size = 6)
})

output$clusters_tree_text <- renderText({ "Data not available." })

observeEvent(input$clusters_tree_info, {
  showModal(
    modalDialog(
      clusters_tree_info$text,
      title = clusters_tree_info$title, easyClose = TRUE, footer = NULL
    )
  )
})

##--------------------------------------------------------------------------##

output$clusters_by_sample_UI <- renderUI({
  if ( !is.null(sample_data()$clusters$by_sample) ) {
    tagList(
      DT::dataTableOutput("clusters_by_sample_table"),
      plotly::plotlyOutput("clusters_by_sample_plot")
    )
  } else {
    textOutput("clusters_by_sample_text")
  }
})

# cell counts by cluster and sample
output$clusters_by_sample_table <- DT::renderDataTable({
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

# bar plot of clusters by sample
output$clusters_by_sample_plot <- plotly::renderPlotly({
  sample_data()$clusters$by_sample %>%
  select(-total_cell_count) %>%
  reshape2::melt(id.vars = "cluster") %>%
  rename(sample = variable, cells = value) %>%
  left_join(
    .,
    sample_data()$clusters$by_sample[ , c("cluster", "total_cell_count") ],
    by = "cluster"
  ) %>%
  mutate(pct = cells / total_cell_count) %>%
  plotly::plot_ly(
    x = ~cluster,
    y = ~pct*100,
    type = "bar",
    color = ~sample,
    colors = sample_data()$samples$colors,
    text = ~pct*100,
    hoverinfo = "name+y"
  ) %>%
  plotly::layout(
    xaxis = list(title = ""),
    yaxis = list(title = "Percentage (%)", hoverformat = ".2f"),
    barmode = "stack",
    hovermode = "compare"
  )
})

# alternative text for bar plot of clusters by sample
output$clusters_by_sample_text <- renderText({
    "Only 1 sample in this data set."
  })

observeEvent(input$clusters_by_sample_info, {
  showModal(
    modalDialog(
      clusters_by_sample_info$text,
      title = clusters_by_sample_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})

##--------------------------------------------------------------------------##

output$clusters_box_nUMI_UI <- renderUI({
  if ( "nUMI" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("clusters_box_nUMI_plot")
  } else {
    textOutput("clusters_box_nUMI_text")
  }
})

# box plot of number of transcripts per cluster
output$clusters_box_nUMI_plot <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~cluster,
    y = ~nUMI,
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
    yaxis = list(title = "Number of UMIs", type = "log", hoverformat = ".2f"),
    dragmode = "select",
    hovermode = "compare"
  )
})

output$clusters_box_nUMI_text <- renderText({
  "Data not available."
})

observeEvent(input$clusters_box_nUMI_info, {
  showModal(
    modalDialog(
      clusters_box_nUMI_info$text,
      title = clusters_box_nUMI_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})

##--------------------------------------------------------------------------##

output$clusters_box_nGene_UI <- renderUI({
  if ( "nGene" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("clusters_box_nGene_plot")
  } else {
    textOutput("clusters_box_nGene_text")
  }
})

# box plot of number of expressed genes per cluster
output$clusters_box_nGene_plot <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~cluster, 
    y = ~nGene,
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
    xaxis = list(title =""),
    yaxis = list(title = "Number of expressed genes", type = "log",
        hoverformat = ".2f"),
    dragmode = "select",
    hovermode = "compare"
  )
})

output$clusters_box_nGene_text <- renderText({
  "Data not available."
})

observeEvent(input$clusters_box_nGene_info, {
  showModal(
    modalDialog(
      clusters_box_nGene_info$text,
      title = clusters_box_nGene_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})

##--------------------------------------------------------------------------##

output$clusters_box_percent_mt_UI <- renderUI({
  if ( "percent_mito" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("clusters_box_percent_mt_plot")
  } else {
    textOutput("clusters_box_percent_mt_text")
  }
})

# box plot of percentage of mitochondrial gene expression per cluster
output$clusters_box_percent_mt_plot <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~cluster,
    y = ~percent_mt*100,
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
    yaxis = list(title = "Percentage of mitochondrial gene expression",
      range = c(0, 100), hoverformat = ".2f"),
    dragmode = "select",
    hovermode = "compare"
  )
})

output$clusters_box_percent_mt_text <- renderText({
  "Data not available."
})

observeEvent(input$clusters_box_percent_mt_info, {
  showModal(
    modalDialog(
      clusters_box_percent_mt_info$text,
      title = clusters_box_percent_mt_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})

##--------------------------------------------------------------------------##

output$clusters_box_percent_ribo_UI <- renderUI({
  if ( "percent_ribo" %in% names(sample_data()$cells) ) {
    plotly::plotlyOutput("clusters_box_percent_ribo_plot")
  } else {
    textOutput("clusters_box_percent_ribo_text")
  }
})

output$clusters_box_percent_ribo_text <- renderText({
  "Data not available."
})

# box plot of percentage of ribosomal gene expression per cluster
output$clusters_box_percent_ribo_plot <- plotly::renderPlotly({
  plotly::plot_ly(
    sample_data()$cells,
    x = ~cluster,
    y = ~percent_ribo*100,
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
    yaxis = list(title = "Percentage of ribosomal gene expression",
      range = c(0, 100), hoverformat = ".2f"),
    dragmode = "select",
    hovermode = "compare"
  )
})

observeEvent(input$clusters_box_percent_ribo_info, {
  showModal(
    modalDialog(
      clusters_box_percent_ribo_info$text,
      title = clusters_box_percent_ribo_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})

##--------------------------------------------------------------------------##

# UI element for bar plot of clusters by cell cycle (Regev)
output$clusters_by_cell_cycle_Regev_UI <- renderUI({
  if ( !is.null(sample_data()$clusters$by_cell_cycle_Regev) ) {
    plotly::plotlyOutput("clusters_by_cell_cycle_Regev_plot")
  } else {
    textOutput("clusters_by_cell_cycle_Regev_text")
  }
})

# bar plot of clusters by cell cycle (Regev)
output$clusters_by_cell_cycle_Regev_plot <- plotly::renderPlotly({
  sample_data()$clusters$by_cell_cycle_Regev %>%
  select(-total_cell_count) %>%
  reshape2::melt(id.vars = "cluster") %>%
  rename(phase = variable, cells = value) %>%
  mutate(
    phase = factor(phase, levels = c("G1", "S", "G2M")),
  ) %>%
  left_join(
    .,
    sample_data()$clusters$by_cell_cycle_Regev[ , c("cluster", "total_cell_count") ],
    by = "cluster"
  ) %>%
  mutate(pct = cells / total_cell_count) %>%
  plotly::plot_ly(
    x = ~cluster,
    y = ~pct*100,
    type = "bar",
    color = ~phase,
    colors = cell_cycle_colorset,
    text = ~pct*100,
    hoverinfo = "name+y"
  ) %>%
  plotly::layout(
    xaxis = list(title = ""),
    yaxis = list(title = "Percentage (%)", hoverformat = ".2f"),
    barmode = "stack",
    hovermode = "compare"
  ) 
})

# alternative text for bar plot of clusters by cell cycle (Regev)
output$clusters_by_cell_cycle_Regev_text <- renderText({
    "Data not available."
  })

observeEvent(input$clusters_by_cell_cycle_Regev_info, {
  showModal(
    modalDialog(
      clusters_by_cell_cycle_Regev_info$text,
      title = clusters_by_cell_cycle_Regev_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})

##--------------------------------------------------------------------------##

# UI element for bar plot of clusters by cell cycle (Cyclone)
output$clusters_by_cell_cycle_Cyclone_UI <- renderUI(
  if ( !is.null(sample_data()$cells$cell_cycle_Cyclone) ) {
    plotly::plotlyOutput("clusters_by_cell_cycle_Cyclone_plot")
  } else {
    textOutput("clusters_by_cell_cycle_Cyclone_text")
  }
)

# bar plot of clusters by cell cycle (Cyclone)
output$clusters_by_cell_cycle_Cyclone_plot <- plotly::renderPlotly({
  sample_data()$clusters$by_cell_cycle_Cyclone %>%
  select(-total_cell_count) %>%
  reshape2::melt(id.vars = "cluster") %>%
  rename(phase = variable, cells = value) %>%
  mutate(
    phase = factor(phase, levels = c("G1", "S", "G2M", "-")),
  ) %>%
  left_join(
    .,
    sample_data()$clusters$by_cell_cycle_Cyclone[ , c("cluster", "total_cell_count") ],
    by = "cluster"
  ) %>%
  mutate(pct = cells / total_cell_count) %>%
  plotly::plot_ly( 
    x = ~cluster,
    y = ~pct*100,
    type = "bar",
    color = ~phase,
    colors = cell_cycle_colorset,
    text = ~pct*100,
    hoverinfo = "name+y"
  ) %>%
  plotly::layout(
    xaxis = list(title = ""),
    yaxis = list(title = "Percentage (%)", hoverformat = ".2f"), 
    barmode = "stack",
    hovermode = "compare"
  )
})

# alternative text for bar plot of clusters by cell cycle (Cyclone)
output$clusters_by_cell_cycle_Cyclone_text <- renderText({
    "Data not available."
  })

observeEvent(input$clusters_by_cell_cycle_Cyclone_info, {
  showModal(
    modalDialog(
      clusters_by_cell_cycle_Cyclone_info$text,
      title = clusters_by_cell_cycle_Cyclone_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})
