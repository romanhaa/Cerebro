##--------------------------------------------------------------------------##
## Panel: Load data.
##--------------------------------------------------------------------------##

output$load_data_number_of_cells <- renderValueBox({
  valueBox(
    value = if ( !is.null(sample_data()) ) nrow(sample_data()$cells) else 0,
    subtitle = "Total number of cells",
    color = "yellow"
  )
})

output$load_data_number_of_samples <- renderValueBox({
  valueBox(
    value = if ( !is.null(sample_data()) ) length(unique(sample_data()$cells$sample)) else 0,
    subtitle = "Number of samples",
    color = "aqua"
  )
})
