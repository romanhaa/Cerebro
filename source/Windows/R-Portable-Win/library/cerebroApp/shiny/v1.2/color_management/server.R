##----------------------------------------------------------------------------##
## Tab: Color Management.
##----------------------------------------------------------------------------##

# UI element
output[["colors_for_samples_UI"]] <- renderUI({
  tagList(
    {
      color_list <- list()
      for ( i in 1:length(sample_data()$sample_names) )
      {
        sample_name <- sample_data()$sample_names[i]
        color_list[[i]] <- colourpicker::colourInput(
          paste0('color_sample_', sample_name),
          sample_name,
          reactive_colors()$samples[i]
        )
      }
      color_list
    }
  )
})

# info box
observeEvent(input[["colors_for_samples_info"]],
{
  showModal(
    modalDialog(
      colors_for_samples_info[["text"]],
      title = colors_for_samples_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})

# UI element
output[["colors_for_clusters_UI"]] <- renderUI({
  tagList(
    {
      color_list <- list()
      for ( i in 1:length(sample_data()$cluster_names) )
      {
        cluster_name <- sample_data()$cluster_names[i]
        color_list[[i]] <- colourpicker::colourInput(
          paste0('color_cluster_', cluster_name),
          paste0('Cluster: ', cluster_name),
          reactive_colors()$clusters[i]
        )
      }
      color_list
    }
  )
})

# info box
observeEvent(input[["colors_for_clusters_info"]],
{
  showModal(
    modalDialog(
      colors_for_clusters_info[["text"]],
      title = colors_for_clusters_info[["title"]],
      easyClose = TRUE,
      footer = NULL
    )
  )
})
