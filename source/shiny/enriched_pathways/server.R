##--------------------------------------------------------------------------##
## Panel: Enriched pathways
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##

# by sample
output$enriched_pathways_by_sample_UI <- renderUI({
  if ( nrow(sample_data()$samples$overview) > 1 & !is.null(sample_data()$marker_genes$by_sample_annotation) ) {    
    tagList(
      fluidRow(
        column(4,
          selectInput("enriched_pathways_select_sample", label = NULL,
            choices = sample_data()$samples$overview$sample)
        ),
        column(8,
          selectInput("enriched_pathways_select_db_for_sample", label = NULL,
            choices = sample_data()$parameters$enrichr_dbs)
        )
      ),
      fluidRow(
        column(12,
          DT::dataTableOutput("enriched_pathways_by_sample_table_present")
        )
      )
    )
  } else {
    textOutput("enriched_pathways_by_sample_table_missing")
  }
})

output$enriched_pathways_by_sample_table_present <- DT::renderDataTable(server = FALSE, {
  req(input$enriched_pathways_select_sample)
  req(input$enriched_pathways_select_db_for_sample)
  sample_data()$marker_genes$by_sample_annotation[[ input$enriched_pathways_select_sample ]][[ input$enriched_pathways_select_db_for_sample ]] %>%
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

output$enriched_pathways_by_sample_table_missing <- renderText({
    "Only 1 sample in this data set or data not available."
  })

observeEvent(input$enriched_pathways_by_sample_info, {
  showModal(
    modalDialog(
      enriched_pathways_by_sample_info$text,
      title = enriched_pathways_by_sample_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})

##--------------------------------------------------------------------------##

# by cluster
output$enriched_pathways_by_cluster_UI <- renderUI({
  if ( nrow(sample_data()$clusters$overview) > 1 & !is.null(sample_data()$marker_genes$by_cluster_annotation) ) {    
    tagList(
      fluidRow(
        column(4,
          selectInput("enriched_pathways_select_cluster", label = NULL,
            choices = sample_data()$clusters$overview$cluster)
        ),
        column(8,
          selectInput("enriched_pathways_select_db_for_cluster", label = NULL,
            choices = sample_data()$parameters$enrichr_dbs)
        )
      ),
      fluidRow(
        column(12,
          DT::dataTableOutput("enriched_pathways_by_cluster_table_present")
        )
      )
    )
  } else {
    textOutput("enriched_pathways_by_cluster_table_missing")
  }
})

output$enriched_pathways_by_cluster_table_present <- DT::renderDataTable(server = FALSE, {
  req(input$enriched_pathways_select_cluster)
  req(input$enriched_pathways_select_db_for_cluster)
  sample_data()$marker_genes$by_cluster_annotation[[ input$enriched_pathways_select_cluster ]][[ input$enriched_pathways_select_db_for_cluster ]][ , c(1,2,3,4,8,9) ] %>%
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

output$enriched_pathways_by_cluster_table_missing <- renderText({
    "Only 1 cluster in this data set or data not available."
  })

observeEvent(input$enriched_pathways_by_cluster_info, {
  showModal(
    modalDialog(
      enriched_pathways_by_cluster_info$text,
      title = enriched_pathways_by_cluster_info$title,
      easyClose = TRUE, footer = NULL
    )
  )
})
