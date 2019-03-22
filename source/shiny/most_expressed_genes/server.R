##----------------------------------------------------------------------------##
## Tab: Most expressed genes.
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## Sample.
##----------------------------------------------------------------------------##

# UI element
output[["most_expressed_genes_by_sample_UI"]] <- renderUI({
  if ( !is.null(sample_data()$most_expressed_genes$by_sample) ) {
    fluidRow(
      column(12,
        selectInput(
          "most_expressed_genes_by_sample_input",
          label = NULL,
          choices = sample_data()$sample_names
        ),
        DT::dataTableOutput("most_expressed_genes_by_sample_table_present")
      )
    )
  } else {
    textOutput("most_expressed_genes_by_sample_table_missing")
  }
})

# table
output[["most_expressed_genes_by_sample_table_present"]] <- DT::renderDataTable(server = FALSE, {
  req(input[["most_expressed_genes_by_sample_input"]])
  sample_data()$most_expressed_genes$by_sample %>%
  filter(sample == input[["most_expressed_genes_by_sample_input"]]) %>%
  mutate(pct = formattable::percent(round(pct/100, digits = 4))) %>%
  rename(
    Sample = sample,
    Gene = gene,
    "% of total expression" = pct
  ) %>%
  formattable::formattable(
    list(
      "Sample" = formattable::color_tile(
          colors[ which(sample_data()$samples$overview$sample == input[["most_expressed_genes_by_sample_input"]]) ],
          colors[ which(sample_data()$samples$overview$sample == input[["most_expressed_genes_by_sample_input"]]) ]
        ),
      "% of total expression" = formattable::color_bar("pink")
    )
  ) %>%
  formattable::as.datatable(
    filter = "none",
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
              filename = "most_expressed_genes_by_sample",
              title = "Most expressed genes by sample"
            ),
            list(
              extend = "excel",
              filename = "most_expressed_genes_by_sample",
              title = "Most expressed genes by sample"
            )
          )
        )
      )
    )
  ) %>%
  DT::formatStyle("% of total expression", textAlign = "right")
})

# alternative text
output[["most_expressed_genes_by_sample_table_missing"]] <- renderText({
    "Data not available."
  })

# info box
observeEvent(input[["most_expressed_genes_by_sample_info"]], {
  showModal(
    modalDialog(
      most_expressed_genes_by_sample_info[["text"]],
      title = most_expressed_genes_by_sample_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Cluster.
##----------------------------------------------------------------------------##

# UI element
output[["most_expressed_genes_by_cluster_UI"]] <- renderUI({
  if ( !is.null(sample_data()$most_expressed_genes$by_cluster) ) {
    fluidRow(
      column(12,
        selectInput(
          "most_expressed_genes_by_cluster_input",
          label = NULL,
          choices = sample_data()$cluster_names
        ),
        DT::dataTableOutput("most_expressed_genes_by_cluster_table_present")
      )
    )
  } else {
    textOutput("most_expressed_genes_by_cluster_table_missing")
  }
})

# table
output[["most_expressed_genes_by_cluster_table_present"]] <- DT::renderDataTable(server = FALSE, {
  req(input[["most_expressed_genes_by_cluster_input"]])
  sample_data()$most_expressed_genes$by_cluster %>%
  filter(cluster == input[["most_expressed_genes_by_cluster_input"]]) %>%
  mutate(pct = formattable::percent(round(pct/100, digits = 4))) %>%
  rename(
    Cluster = cluster,
    Gene = gene,
    "% of total expression" = pct
  ) %>%
  formattable::formattable(
    list(
      "Cluster" = formattable::color_tile(
          colors[ which(sample_data()$clusters$overview$cluster == input[["most_expressed_genes_by_cluster_input"]]) ],
          colors[ which(sample_data()$clusters$overview$cluster == input[["most_expressed_genes_by_cluster_input"]]) ]
        ),
      "% of total expression" = formattable::color_bar("pink")
    )
  ) %>%
  formattable::as.datatable(
    filter = "none",
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
              filename = "most_expressed_genes_by_cluster",
              title = "Most expressed genes by cluster"
            ),
            list(
              extend = "excel",
              filename = "most_expressed_genes_by_cluster",
              title = "Most expressed genes by cluster"
            )
          )
        )
      )
    )
  ) %>%
  DT::formatStyle("% of total expression", textAlign = "right")
})

# alternative text
output[["most_expressed_genes_by_cluster_table_missing"]] <- renderText({
  "Data not available."
})

# info box
observeEvent(input[["most_expressed_genes_by_cluster_info"]], {
  showModal(
    modalDialog(
      most_expressed_genes_by_cluster_info[["text"]],
      title = most_expressed_genes_by_cluster_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})
