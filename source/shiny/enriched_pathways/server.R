##----------------------------------------------------------------------------##
## Tab: Enriched pathways
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## Samples.
##----------------------------------------------------------------------------##

# UI element: choose sample
output[["enriched_pathways_by_sample_select_sample_UI"]] <- renderUI({
  selectInput(
    "enriched_pathways_by_sample_select_sample",
    label = NULL,
    choices = sample_data()$sample_names
  )
})

# UI element: choose database
output[["enriched_pathways_by_sample_select_db_UI"]] <- renderUI({
  req(input[["enriched_pathways_by_sample_select_sample"]])
  selectInput(
    "enriched_pathways_by_sample_select_db",
    label = NULL,
    choices = names(sample_data()$marker_genes$by_sample_annotation[[ input[["enriched_pathways_by_sample_select_sample"]] ]])
  )
})

# UI element: display results
output[["enriched_pathways_by_sample_UI"]] <- renderUI({
  req(input[["enriched_pathways_by_sample_select_sample"]])
  req(input[["enriched_pathways_by_sample_select_db"]])
  if ( !is.null(sample_data()$marker_genes$by_sample_annotation) ) {
    DT::dataTableOutput("enriched_pathways_by_sample_table_present")
  } else {
    textOutput("enriched_pathways_by_sample_table_missing")
  }
})

# table
output[["enriched_pathways_by_sample_table_present"]] <- DT::renderDataTable(server = FALSE, {
  req(input[["enriched_pathways_by_sample_select_sample"]])
  req(input[["enriched_pathways_by_sample_select_db"]])
  sample_data()$marker_genes$by_sample_annotation[[ input[["enriched_pathways_by_sample_select_sample"]] ]][[ input[["enriched_pathways_by_sample_select_db"]] ]] %>%
  select(c(1,2,3,4,8,9)) %>%
  mutate(
    P.value = formatC(P.value, format = "e", digits = 3),
    Adjusted.P.value = formatC(Adjusted.P.value, format = "e", digits = 3),
    Combined.Score = formatC(Combined.Score, format = "f", digits = 2)
  ) %>%
  rename(
    "p-value" = P.value,
    "adj. p-value" = Adjusted.P.value,
    "combined score" = Combined.Score,
  ) %>%
  formattable::formattable(
    list("combined score" = formattable::color_bar("pink"))
  ) %>%
  formattable::as.datatable(
    filter = "top",
    selection = "multiple",
    escape = FALSE,
    autoHideNavigation = TRUE,
    rownames = FALSE,
    extensions = c("Buttons"),
    class = "cell-border stripe",
    options = list(
      columnDefs = list(list(visible = FALSE, targets = c(2,5))),
      scrollX = TRUE,
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
              filename = "enriched_pathways_by_sample",
              title = "Enriched pathways by sample"
            ),
            list(
              extend = "excel",
              filename = "enriched_pathways_by_sample",
              title = "Enriched pathways by sample"
            ),
            list(
              extend = "pdf",
              filename = "enriched_pathways_by_sample",
              title = "Enriched pathways by sample"
            )
          )
        )
      )
    )
  ) %>% 
  DT::formatStyle(columns = c("combined score"), textAlign = "right")
})

# alternative text
output[["enriched_pathways_by_sample_table_missing"]] <- renderText({
  "Only 1 sample in this data set or data not available."
})

# info box
observeEvent(input[["enriched_pathways_by_sample_info"]], {
  showModal(
    modalDialog(
      enriched_pathways_by_sample_info[["text"]],
      title = enriched_pathways_by_sample_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

##----------------------------------------------------------------------------##
## Clusters.
##----------------------------------------------------------------------------##

# UI element: choose cluster
output[["enriched_pathways_by_cluster_select_cluster_UI"]] <- renderUI({
  selectInput(
    "enriched_pathways_by_cluster_select_cluster",
    label = NULL,
    choices = sample_data()$cluster_names
  )
})

# UI element: choose database
output[["enriched_pathways_by_cluster_select_db_UI"]] <- renderUI({
  req(input[["enriched_pathways_by_cluster_select_cluster"]])
  selectInput(
    "enriched_pathways_by_cluster_select_db",
    label = NULL,
    choices = names(sample_data()$marker_genes$by_cluster_annotation[[ input[["enriched_pathways_by_cluster_select_cluster"]] ]])
  )
})

# UI element: display results
output$enriched_pathways_by_cluster_UI <- renderUI({
  req(input[["enriched_pathways_by_cluster_select_cluster"]])
  req(input[["enriched_pathways_by_cluster_select_db"]])
  if ( !is.null(sample_data()$marker_genes$by_cluster_annotation) ) {
    DT::dataTableOutput("enriched_pathways_by_cluster_table_present")
  } else {
    textOutput("enriched_pathways_by_cluster_table_missing")
  }
})

# table
output[["enriched_pathways_by_cluster_table_present"]] <- DT::renderDataTable(server = FALSE, {
  req(input[["enriched_pathways_by_cluster_select_cluster"]])
  req(input[["enriched_pathways_by_cluster_select_db"]])
  sample_data()$marker_genes$by_cluster_annotation[[ input[["enriched_pathways_by_cluster_select_cluster"]] ]][[ input[["enriched_pathways_by_cluster_select_db"]] ]] %>%
  select(c(1,2,3,4,8,9)) %>%
  mutate(
    P.value = formatC(P.value, format = "e", digits = 3),
    Adjusted.P.value = formatC(Adjusted.P.value, format = "e", digits = 3),
    Combined.Score = formatC(Combined.Score, format = "f", digits = 2)
  ) %>%
  rename(
    "p-value" = P.value,
    "adj. p-value" = Adjusted.P.value,
    "combined score" = Combined.Score,
  ) %>%
  formattable::formattable(
    list("combined score" = formattable::color_bar("pink"))
  ) %>%
  formattable::as.datatable(
    filter = "top",
    selection = "multiple",
    escape = FALSE,
    autoHideNavigation = TRUE,
    rownames = FALSE,
    extensions = c("Buttons"),
    class = "cell-border stripe",
    options = list(
      columnDefs = list(list(visible = FALSE, targets = c(2,5))),
      scrollX = TRUE,
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
              filename = "enriched_pathways_by_cluster",
              title = "Enriched pathways by cluster"
            ),
            list(
              extend = "excel",
              filename = "enriched_pathways_by_cluster",
              title = "Enriched pathways by cluster"
            ),
            list(
              extend = "pdf",
              filename = "enriched_pathways_by_cluster",
              title = "Enriched pathways by cluster"
            )
          )
        )
      )
    )
  ) %>% 
  DT::formatStyle(columns = c("combined score"), textAlign = "right")
})

# alternative text
output[["enriched_pathways_by_cluster_table_missing"]] <- renderText({
  "Only 1 cluster in this data set or data not available."
})

# info box
observeEvent(input[["enriched_pathways_by_cluster_info"]], {
  showModal(
    modalDialog(
      enriched_pathways_by_cluster_info[["text"]],
      title = enriched_pathways_by_cluster_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})
