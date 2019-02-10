##--------------------------------------------------------------------------##
## Tab: Load data.
##--------------------------------------------------------------------------##

output$load_data_experiment_name <- renderValueBox({
  # valueBox(
  #   value = if ( !is.null(sample_data()$experiment$experiment_name)) sample_data()$experiment$experiment_name else "-",
  #   subtitle = "Experiment",
  #   color = "green"
  # )
  box(
    title = "Experiment",
    width = 5,
    background = "light-blue",
    sample_data()$experiment$experiment_name
  )
})

output$load_data_organism <- renderValueBox({
  # valueBox(
  #   value = if ( !is.null(sample_data()$experiment$organism)) sample_data()$experiment$organism else "-",
  #   subtitle = "Organism",
  #   color = "orange",
  #   icon = icon("car")
  # )
  box(
    title = "Organism",
    width = 5,
    background = "light-blue",
    sample_data()$experiment$organism
  )
})

output$load_data_number_of_cells <- renderValueBox({
  valueBox(
    value = if ( !is.null(sample_data()) ) nrow(sample_data()$cells) else 0,
    subtitle = "cells",
    color = "light-blue"
  )
})

output$load_data_number_of_samples <- renderValueBox({
  valueBox(
    value = if ( !is.null(sample_data()) ) length(unique(sample_data()$cells$sample)) else 0,
    subtitle = "samples",
    color = "light-blue"
  )
})

output$load_data_number_of_clusters <- renderValueBox({
  valueBox(
    value = if ( !is.null(sample_data()) ) length(unique(sample_data()$cells$cluster)) else 0,
    subtitle = "clusters",
    color = "light-blue"
  )
})
