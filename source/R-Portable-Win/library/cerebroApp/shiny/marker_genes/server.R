##----------------------------------------------------------------------------##
## Tab: Marker genes.
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## Sample.
##----------------------------------------------------------------------------##

# UI element
output[["marker_genes_by_sample_UI"]] <- renderUI({
  if ( !is.null(sample_data()$marker_genes$by_sample) ) {
    if ( is.data.frame(sample_data()$marker_genes$by_sample) ) {
      fluidRow(
        column(12,
          selectInput(
            "marker_genes_by_sample_input",
            label = NULL,
            choices = unique(sample_data()$marker_genes$by_sample$sample)
          ),
          DT::dataTableOutput("marker_genes_by_sample_table_present")
        )
      )
    } else if ( sample_data()$marker_genes$by_sample == "no_markers_found" ) {
      textOutput("marker_genes_by_sample_table_no_markers_found")
    }
  } else {
    textOutput("marker_genes_by_sample_table_missing")
  }
})

# table
output[["marker_genes_by_sample_table_present"]] <- DT::renderDataTable(server = FALSE, {
  req(input[["marker_genes_by_sample_input"]])
  if ( "on_cell_surface" %in% colnames(sample_data()$marker_genes$by_sample) ) {
    table <- sample_data()$marker_genes$by_sample[ which(sample_data()$marker_genes$by_sample$sample == input[["marker_genes_by_sample_input"]]) , c(2,4,5,6,7,8) ]
    colnames(table) <- c("Gene", "avg. logFC", "% cells in this sample", "% cells in other samples", "adj. p-value", "present on cell surface")
    table$"avg. logFC" <- round(table$"avg. logFC", digits=3)
    table$"% cells in this sample" <- formattable::percent(table$"% cells in this sample")
    table$"% cells in other samples" <- formattable::percent(table$"% cells in other samples")
    table <- table %>%
    mutate("adj. p-value"=formatC(table$"adj. p-value", format="e", digits=3)) %>%
    formattable::formattable(list(
      "avg. logFC" = formattable::color_tile("white", "orange"),
      "% cells in this sample" = formattable::color_bar("pink"),
      "% cells in other samples" = formattable::color_bar("pink"),
      "present on cell surface" = formattable::formatter("span", style=x~style(color=ifelse(x, "green", "red")))
    ))
  } else {
    table <- sample_data()$marker_genes$by_sample[ which(sample_data()$marker_genes$by_sample$sample == input[["marker_genes_by_sample_input"]]) , c(2,4,5,6,7) ]
    colnames(table) <- c("Gene", "avg. logFC", "% cells in this sample", "% cells in other samples", "adj. p-value")
    table$"avg. logFC" <- round(table$"avg. logFC", digits=3)
    table$"% cells in this sample" <- formattable::percent(table$"% cells in this sample")
    table$"% cells in other samples" <- formattable::percent(table$"% cells in other samples")
    table <- table %>%
    mutate("adj. p-value"=formatC(table$"adj. p-value", format="e", digits=3)) %>%
    formattable::formattable(list(
      "avg. logFC" = formattable::color_tile("white", "orange"),
      "% cells in this sample" = formattable::color_bar("pink"),
      "% cells in other samples" = formattable::color_bar("pink")
    ))
  }
  formattable::as.datatable(
    table,
    filter = "top",
    selection = "multiple",
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
              filename = "marker_genes_by_sample",
              title = "Marker genes by sample"
            ),
            list(
              extend = "excel",
              filename = "marker_genes_by_sample",
              title = "Marker genes by sample"
            )
          )
        )
      )
    )
  ) %>% 
  DT::formatStyle(
    columns = c("avg. logFC", "% cells in this sample", "% cells in other samples", "adj. p-value"),
    textAlign = "right"
  )
})

# alternative text
output[["marker_genes_by_sample_table_no_markers_found"]] <- renderText({
  "No marker genes identified for any of the samples."
})

# alternative text
output[["marker_genes_by_sample_table_missing"]] <- renderText({
  "Data not available. Possible reasons: Only 1 sample in this data set or data not generated."
})

# info box
observeEvent(input[["marker_genes_by_sample_info"]], {
  showModal(
    modalDialog(
      marker_genes_by_sample_info[["text"]],
      title = marker_genes_by_sample_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Cluster.
##----------------------------------------------------------------------------##

# UI element
output[["marker_genes_by_cluster_UI"]] <- renderUI({
  if ( !is.null(sample_data()$marker_genes$by_cluster) ) {
    if ( is.data.frame(sample_data()$marker_genes$by_cluster) ) {
      fluidRow(
        column(12,
          selectInput(
            "marker_genes_by_cluster_input",
            label = NULL,
            choices = unique(sample_data()$marker_genes$by_cluster$cluster)
          ),
          DT::dataTableOutput("marker_genes_by_cluster_table_present")
        )
      )
    } else if ( sample_data()$marker_genes$by_cluster == "no_markers_found" ) {
      textOutput("marker_genes_by_cluster_table_no_markers_found")
    }
  } else {
    textOutput("marker_genes_by_cluster_table_missing")
  }
})

# table
output[["marker_genes_by_cluster_table_present"]] <- DT::renderDataTable(server = FALSE, {
  req(input[["marker_genes_by_cluster_input"]])
  if ("on_cell_surface" %in% colnames(sample_data()$marker_genes$by_cluster)) {
    table <- sample_data()$marker_genes$by_cluster[ which(sample_data()$marker_genes$by_cluster$cluster == input[["marker_genes_by_cluster_input"]]) , c(2,4,5,6,7,8) ]
    colnames(table) <- c("Gene", "avg. logFC", "% cells in this cluster", "% cells in other clusters", "adj. p-value", "present on cell surface")
    table$"avg. logFC" <- round(table$"avg. logFC", digits=3)
    table$"% cells in this cluster" <- formattable::percent(table$"% cells in this cluster")
    table$"% cells in other clusters" <- formattable::percent(table$"% cells in other clusters")
    table <- table %>% mutate("adj. p-value"=formatC(table$"adj. p-value", format="e", digits=3)) %>%
    formattable::formattable(list(
      "avg. logFC" = formattable::color_tile("white", "orange"),
      "% cells in this cluster" = formattable::color_bar("pink"),
      "% cells in other clusters" = formattable::color_bar("pink"),
      "present on cell surface" = formattable::formatter("span", style=x~style(color=ifelse(x, "green", "red")))
    ))
  } else {
    table <- sample_data()$marker_genes$by_cluster[ which(sample_data()$marker_genes$by_cluster$cluster == input[["marker_genes_by_cluster_input"]]) , c(2,4,5,6,7) ]
    colnames(table) <- c("Gene", "avg. logFC", "% cells in this cluster", "% cells in other clusters", "adj. p-value")
    table$"avg. logFC" <- round(table$"avg. logFC", digits=3)
    table$"% cells in this cluster" <- formattable::percent(table$"% cells in this cluster")
    table$"% cells in other clusters" <- formattable::percent(table$"% cells in other clusters")
    table <- table %>% mutate("adj. p-value"=formatC(table$"adj. p-value", format="e", digits=3)) %>%
    formattable::formattable(list(
      "avg. logFC" = formattable::color_tile("white", "orange"),
      "% cells in this cluster" = formattable::color_bar("pink"),
      "% cells in other clusters" = formattable::color_bar("pink")
    ))
  }
  formattable::as.datatable(
    table,
    filter = "top",
    selection = "multiple",
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
              filename = "marker_genes_by_cluster",
              title = "Marker genes by cluster"
            ),
            list(
              extend = "excel",
              filename = "marker_genes_by_cluster",
              title = "Marker genes by cluster"
            )
          )
        )
      )
    )
  ) %>% 
  DT::formatStyle(
    columns = c("avg. logFC", "% cells in this cluster", "% cells in other clusters", "adj. p-value"),
    textAlign = "right"
  )
})

# alternative text
output[["marker_genes_by_cluster_table_no_markers_found"]] <- renderText({
  "No marker genes identified for any of the clusters."
})

# alternative text
output[["marker_genes_by_cluster_table_missing"]] <- renderText({
  "Data not available. Possible reasons: Only 1 cluster in this data set or data not generated."
})

# info box
observeEvent(input[["marker_genes_by_cluster_info"]], {
  showModal(
    modalDialog(
      marker_genes_by_cluster_info[["text"]],
      title = marker_genes_by_cluster_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})
