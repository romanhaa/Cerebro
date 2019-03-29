##----------------------------------------------------------------------------##
## Tab: Load data.
##----------------------------------------------------------------------------##

# number of cells
output$load_data_number_of_cells <- renderValueBox({
  valueBox(
    value = if ( !is.null(sample_data()) ) formatC(nrow(sample_data()$cells), format = "f", big.mark = ",", digits = 0) else 0,
    subtitle = "cells",
    color = "light-blue"
  )
})

# number of samples
output$load_data_number_of_samples <- renderValueBox({
  valueBox(
    value = ifelse(!is.null(sample_data()), length(sample_data()$sample_names), 0),
    subtitle = "samples",
    color = "light-blue"
  )
})

# number of clusters
output$load_data_number_of_clusters <- renderValueBox({
  valueBox(
    value = ifelse(!is.null(sample_data()), length(sample_data()$cluster_names), 0),
    subtitle = "clusters",
    color = "light-blue"
  )
})

# experiment name
output$load_data_experiment_name <- renderValueBox({
  box(
    title = "Experiment",
    width = 5,
    background = "light-blue",
    sample_data()$experiment$experiment_name
  )
})

# organism
output$load_data_organism <- renderValueBox({
  box(
    title = "Organism",
    width = 5,
    background = "light-blue",
    sample_data()$experiment$organism
  )
})


