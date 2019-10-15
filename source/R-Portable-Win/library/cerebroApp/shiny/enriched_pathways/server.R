##----------------------------------------------------------------------------##
## Tab: Enriched pathways
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## Samples.
##----------------------------------------------------------------------------##

# UI element: choose source for pathway enrichement results (currently Enrichr or GSVA)
output[["enriched_pathways_by_sample_select_source_UI"]] <- renderUI({
  if ( is.null(sample_data()$enriched_pathways) ) {
    textOutput("enriched_pathways_by_sample_table_missing")
  } else {
    selectInput(
      "enriched_pathways_by_sample_select_source",
      label = NULL,
      choices = names(sample_data()$enriched_pathways)
    )
  }
})

# UI element: display results or alternative text
output[["enriched_pathways_by_sample_UI"]] <- renderUI({
  req(input[["enriched_pathways_by_sample_select_source"]])
  if ( input[["enriched_pathways_by_sample_select_source"]] == "enrichr" ) {
    if ( !is.null(sample_data()$enriched_pathways$enrichr$by_sample) ) {
      if ( is.list(sample_data()$enriched_pathways$enrichr$by_sample) ) {
        tagList(
          fluidRow(
            column(4,
              uiOutput("enriched_pathways_by_sample_select_sample_UI")
            ),
            column(8,
              uiOutput("enriched_pathways_by_sample_select_db_UI")
            )
          ),
          DT::dataTableOutput("enriched_pathways_by_sample_table_present")
        )
      } else if ( sample_data()$enriched_pathways$enrichr$by_sample == "no_markers_found" ) {
        textOutput("enriched_pathways_by_sample_table_no_markers_found")
      }
    } else {
      textOutput("enriched_pathways_by_sample_table_missing_enrichr")
    }
  } else if ( input[["enriched_pathways_by_sample_select_source"]] == "GSVA" ) {
    if ( !is.null(sample_data()$enriched_pathways$GSVA$by_sample) ) {
      if ( is.data.frame(sample_data()$enriched_pathways$GSVA$by_sample) ) {
        tagList(
          fluidRow(
            column(4,
              uiOutput("enriched_pathways_by_sample_select_sample_UI")
            )
          ),
          DT::dataTableOutput("enriched_pathways_by_sample_table_present")
        )
      } else if ( sample_data()$enriched_pathways$GSVA$by_sample == "no_gene_sets_enriched" ) {
        textOutput("enriched_pathways_by_sample_table_no_gene_sets_enriched")
      } else if ( sample_data()$enriched_pathways$GSVA$by_sample == "only_one_sample_in_data_set" ) {
        textOutput("enriched_pathways_by_sample_table_only_one_sample_in_data_set")
      }
    } else {
      textOutput("enriched_pathways_by_sample_table_missing_gsva")
    }
  }
})

# UI element: choose sample
output[["enriched_pathways_by_sample_select_sample_UI"]] <- renderUI({
  req(input[["enriched_pathways_by_sample_select_source"]])
  if ( input[["enriched_pathways_by_sample_select_source"]] == 'enrichr' ) {
    choices <- levels(sample_data()$enriched_pathways$enrichr$by_sample$sample) %>%
      intersect(., unique(sample_data()$enriched_pathways$enrichr$by_sample$sample))
  } else if ( input[["enriched_pathways_by_sample_select_source"]] == 'GSVA' ) {
    choices <- levels(sample_data()$enriched_pathways$GSVA$by_sample$group) %>%
      intersect(., unique(sample_data()$enriched_pathways$GSVA$by_sample$group))
  }
  selectInput(
    "enriched_pathways_by_sample_select_sample",
    label = NULL,
    choices = choices
  )
})

# UI element: choose database
output[["enriched_pathways_by_sample_select_db_UI"]] <- renderUI({
  req(
    input[["enriched_pathways_by_sample_select_source"]],
    input[["enriched_pathways_by_sample_select_sample"]]
  )
  choices <- sample_data()$enriched_pathways$enrichr$by_sample %>%
    dplyr::filter(sample == input[["enriched_pathways_by_sample_select_sample"]]) %>%
    dplyr::pull(db) %>%
    intersect(., levels(.))
  selectInput(
    "enriched_pathways_by_sample_select_db",
    label = NULL,
    choices = choices
  )
})

# table
output[["enriched_pathways_by_sample_table_present"]] <- DT::renderDataTable(server = FALSE, {
  req(
    input[["enriched_pathways_by_sample_select_source"]],
    input[["enriched_pathways_by_sample_select_sample"]],
    input[["enriched_pathways_by_sample_select_db"]]
  )
  if ( input[["enriched_pathways_by_sample_select_source"]] == "enrichr" & is.data.frame(sample_data()$enriched_pathways$enrichr$by_sample) ) {
    sample_data()$enriched_pathways$enrichr$by_sample %>%
    dplyr::filter(
      sample == input[["enriched_pathways_by_sample_select_sample"]],
      db == input[["enriched_pathways_by_sample_select_db"]]
    ) %>%
    dplyr::select(3,4,5,6,10,11) %>%
    dplyr::arrange(-Combined.Score) %>%
    dplyr::mutate(
      P.value = formatC(P.value, format = "e", digits = 2),
      Adjusted.P.value = formatC(Adjusted.P.value, format = "e", digits = 2),
      Combined.Score = formatC(Combined.Score, format = "f", digits = 2)
    ) %>%
    dplyr::rename(
      "p-value" = P.value,
      "adj. p-value" = Adjusted.P.value,
      "combined score" = Combined.Score,
    ) %>%
    formattable::formattable(
      list("combined score" = formattable::color_bar("pink"))
    ) %>%
    formattable::as.datatable(
      filter = "top",
      selection = "none",
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
                filename = "enriched_pathways_from_enrichr_by_sample",
                title = "Enriched pathways from Enrichr by sample"
              ),
              list(
                extend = "excel",
                filename = "enriched_pathways_from_enrichr_by_sample",
                title = "Enriched pathways from Enrichr by sample"
              )
            )
          )
        )
      )
    ) %>%
    DT::formatStyle(columns = c("combined score"), textAlign = "right")
  } else if ( input[["enriched_pathways_by_sample_select_source"]] == "GSVA" & is.data.frame(sample_data()$enriched_pathways$GSVA$by_sample) ) {
    sample_data()$enriched_pathways$GSVA$by_sample %>%
    dplyr::filter(group == input[["enriched_pathways_by_sample_select_sample"]]) %>%
    dplyr::select(-group) %>%
    dplyr::arrange(q_value) %>%
    dplyr::mutate(
      p_value = formatC(p_value, format = "e", digits = 2),
      q_value = formatC(q_value, format = "e", digits = 2)
    ) %>%
    dplyr::rename(
      "Gene set name" = name,
      "Description" = description,
      "Number of genes" = length,
      "Genes" = genes,
      "Enrichment score" = enrichment_score,
      "p-value" = p_value,
      "adj. p-value" = q_value
    ) %>%
    formattable::formattable() %>%
    formattable::as.datatable(
      filter = "top",
      selection = "none",
      escape = FALSE,
      autoHideNavigation = TRUE,
      rownames = FALSE,
      extensions = c("Buttons"),
      class = "cell-border stripe",
      options = list(
        columnDefs = list(list(visible = FALSE, targets = c(3,4))),
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
                filename = "enriched_pathways_from_GSVA_by_sample",
                title = "Enriched pathways from GSVA by sample"
              ),
              list(
                extend = "excel",
                filename = "enriched_pathways_from_GSVA_by_sample",
                title = "Enriched pathways from GSVA by sample"
              )
            )
          )
        )
      )
    )
  }
})

# alternative text messages
output[["enriched_pathways_by_sample_table_missing"]] <- renderText({
  "Data not available. Possible reason: Data not generated."
})

output[["enriched_pathways_by_sample_table_no_markers_found"]] <- renderText({
  "No marker genes identified to perform pathway enrichment analysis with."
})

output[["enriched_pathways_by_sample_table_missing_enrichr"]] <- renderText({
  "Data not available. Possible reasons: Only 1 sample in this data set, no marker genes found or data not generated."
})

output[["enriched_pathways_by_sample_table_no_gene_sets_enriched"]] <- renderText({
  "No gene sets were found to be enriched (with the selected statistical thresholds) in any sample."
})

output[["enriched_pathways_by_sample_table_only_one_sample_in_data_set"]] <- renderText({
  "The loaded data set consists of a single sample which means GSVA cannot be applied."
})

output[["enriched_pathways_by_sample_table_missing_gsva"]] <- renderText({
  "Data not available. Possible reason: Data not generated."
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

# UI element: choose source for pathway enrichement results (currently Enrichr or GSVA)
output[["enriched_pathways_by_cluster_select_source_UI"]] <- renderUI({
  if ( is.null(sample_data()$enriched_pathways) ) {
    textOutput("enriched_pathways_by_cluster_table_missing")
  } else {
    selectInput(
      "enriched_pathways_by_cluster_select_source",
      label = NULL,
      choices = names(sample_data()$enriched_pathways)
    )
  }
})

# UI element: display results or alternative text
output[["enriched_pathways_by_cluster_UI"]] <- renderUI({
  req(input[["enriched_pathways_by_cluster_select_source"]])
  if ( input[["enriched_pathways_by_cluster_select_source"]] == "enrichr" ) {
    if ( !is.null(sample_data()$enriched_pathways$enrichr$by_cluster) ) {
      if ( is.list(sample_data()$enriched_pathways$enrichr$by_cluster) ) {
        tagList(
          fluidRow(
            column(4,
              uiOutput("enriched_pathways_by_cluster_select_cluster_UI")
            ),
            column(8,
              uiOutput("enriched_pathways_by_cluster_select_db_UI")
            )
          ),
          DT::dataTableOutput("enriched_pathways_by_cluster_table_present")
        )
      } else if ( sample_data()$enriched_pathways$enrichr$by_cluster == "no_markers_found" ) {
        textOutput("enriched_pathways_by_cluster_table_no_markers_found")
      }
    } else {
      textOutput("enriched_pathways_by_cluster_table_missing_enrichr")
    }
  } else if ( input[["enriched_pathways_by_cluster_select_source"]] == "GSVA" ) {
    if ( !is.null(sample_data()$enriched_pathways$GSVA$by_cluster) ) {
      if ( is.data.frame(sample_data()$enriched_pathways$GSVA$by_cluster) ) {
        tagList(
          fluidRow(
            column(4,
              uiOutput("enriched_pathways_by_cluster_select_cluster_UI")
            )
          ),
          DT::dataTableOutput("enriched_pathways_by_cluster_table_present")
        )
      } else if ( sample_data()$enriched_pathways$GSVA$by_cluster == "no_gene_sets_enriched" ) {
        textOutput("enriched_pathways_by_cluster_table_no_gene_sets_enriched")
      } else if ( sample_data()$enriched_pathways$GSVA$by_cluster == "only_one_cluster_in_data_set" ) {
        textOutput("enriched_pathways_by_cluster_table_only_one_cluster_in_data_set")
      }
    } else {
      textOutput("enriched_pathways_by_cluster_table_missing_gsva")
    }
  }
})

# UI element: choose cluster
output[["enriched_pathways_by_cluster_select_cluster_UI"]] <- renderUI({
  req(input[["enriched_pathways_by_cluster_select_source"]])
  if ( input[["enriched_pathways_by_cluster_select_source"]] == 'enrichr' ) {
    choices <- levels(sample_data()$enriched_pathways$enrichr$by_cluster$cluster) %>%
      intersect(., unique(sample_data()$enriched_pathways$enrichr$by_cluster$cluster))
  } else if ( input[["enriched_pathways_by_cluster_select_source"]] == 'GSVA' ) {
    choices <- levels(sample_data()$enriched_pathways$GSVA$by_cluster$group) %>%
      intersect(., unique(sample_data()$enriched_pathways$GSVA$by_cluster$group))
  }
  selectInput(
    "enriched_pathways_by_cluster_select_cluster",
    label = NULL,
    choices = choices
  )
})

# UI element: choose database
output[["enriched_pathways_by_cluster_select_db_UI"]] <- renderUI({
  req(
    input[["enriched_pathways_by_cluster_select_source"]],
    input[["enriched_pathways_by_cluster_select_cluster"]]
  )
  choices <- sample_data()$enriched_pathways$enrichr$by_cluster %>%
    dplyr::filter(cluster == input[["enriched_pathways_by_cluster_select_cluster"]]) %>%
    dplyr::pull(db) %>%
    intersect(., levels(.))
  selectInput(
    "enriched_pathways_by_cluster_select_db",
    label = NULL,
    choices = choices
  )
})

# table
output[["enriched_pathways_by_cluster_table_present"]] <- DT::renderDataTable(server = FALSE, {
  req(
    input[["enriched_pathways_by_cluster_select_source"]],
    input[["enriched_pathways_by_cluster_select_cluster"]],
    input[["enriched_pathways_by_cluster_select_db"]]
  )
  if ( input[["enriched_pathways_by_cluster_select_source"]] == "enrichr" & is.data.frame(sample_data()$enriched_pathways$enrichr$by_cluster) ) {
    sample_data()$enriched_pathways$enrichr$by_cluster %>%
    dplyr::filter(
      cluster == input[["enriched_pathways_by_cluster_select_cluster"]],
      db == input[["enriched_pathways_by_cluster_select_db"]]
    ) %>%
    dplyr::select(3,4,5,6,10,11) %>%
    dplyr::arrange(-Combined.Score) %>%
    dplyr::mutate(
      P.value = formatC(P.value, format = "e", digits = 2),
      Adjusted.P.value = formatC(Adjusted.P.value, format = "e", digits = 2),
      Combined.Score = formatC(Combined.Score, format = "f", digits = 2)
    ) %>%
    dplyr::rename(
      "p-value" = P.value,
      "adj. p-value" = Adjusted.P.value,
      "combined score" = Combined.Score,
    ) %>%
    formattable::formattable(
      list("combined score" = formattable::color_bar("pink"))
    ) %>%
    formattable::as.datatable(
      filter = "top",
      selection = "none",
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
                filename = "enriched_pathways_from_enrichr_by_cluster",
                title = "Enriched pathways from Enrichr by cluster"
              ),
              list(
                extend = "excel",
                filename = "enriched_pathways_from_enrichr_by_cluster",
                title = "Enriched pathways from Enrichr by cluster"
              )
            )
          )
        )
      )
    ) %>%
    DT::formatStyle(columns = c("combined score"), textAlign = "right")
  } else if ( input[["enriched_pathways_by_cluster_select_source"]] == "GSVA" & is.data.frame(sample_data()$enriched_pathways$GSVA$by_cluster) ) {
    sample_data()$enriched_pathways$GSVA$by_cluster %>%
    dplyr::filter(group == input[["enriched_pathways_by_cluster_select_cluster"]]) %>%
    dplyr::select(-group) %>%
    dplyr::arrange(q_value) %>%
    dplyr::mutate(
      p_value = formatC(p_value, format = "e", digits = 2),
      q_value = formatC(q_value, format = "e", digits = 2)
    ) %>%
    dplyr::rename(
      "Gene set name" = name,
      "Description" = description,
      "Number of genes" = length,
      "Genes" = genes,
      "Enrichment score" = enrichment_score,
      "p-value" = p_value,
      "adj. p-value" = q_value
    ) %>%
    formattable::formattable() %>%
    formattable::as.datatable(
      filter = "top",
      selection = "none",
      escape = FALSE,
      autoHideNavigation = TRUE,
      rownames = FALSE,
      extensions = c("Buttons"),
      class = "cell-border stripe",
      options = list(
        columnDefs = list(list(visible = FALSE, targets = c(3,4))),
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
                filename = "enriched_pathways_from_GSVA_by_cluster",
                title = "Enriched pathways from GSVA by cluster"
              ),
              list(
                extend = "excel",
                filename = "enriched_pathways_from_GSVA_by_cluster",
                title = "Enriched pathways from GSVA by cluster"
              )
            )
          )
        )
      )
    )
  }
})

# alternative text messages
output[["enriched_pathways_by_cluster_table_missing"]] <- renderText({
  "Data not available. Possible reason: Data not generated."
})

output[["enriched_pathways_by_cluster_table_no_markers_found"]] <- renderText({
  "No marker genes identified to perform pathway enrichment analysis with."
})

output[["enriched_pathways_by_cluster_table_missing_enrichr"]] <- renderText({
  "Data not available. Possible reasons: Only 1 cluster in this data set, no marker genes found or data not generated."
})

output[["enriched_pathways_by_cluster_table_no_gene_sets_enriched"]] <- renderText({
  "Either the loaded data set consists of a single cluster (in which case GSVA cannot be applied) or no gene sets were found to be enriched (with the selected statistical thresholds) in any cluster."
})

output[["enriched_pathways_by_cluster_table_only_one_cluster_in_data_set"]] <- renderText({
  "The loaded data set consists of a single cluster which means GSVA cannot be applied."
})

output[["enriched_pathways_by_cluster_table_missing_gsva"]] <- renderText({
  "Data not available. Possible reason: Data not generated."
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
