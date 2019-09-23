##----------------------------------------------------------------------------##
## Tab: Load data.
##----------------------------------------------------------------------------##

tab_load_data <- tabItem(
    tabName = "loadData",
    fluidRow(
      column(12,
        titlePanel("Load data"),
        fileInput(
          inputId = "input_file",
          label = "Select input data (.crb or .rds file)...",
          multiple = FALSE,
          accept = c(".rds",".crb",".cerebro"),
          width = '350px',
          buttonLabel = "Browse...",
          placeholder = "No file selected"
        )
      )
    ),
    fluidRow(
      valueBoxOutput("load_data_number_of_cells"),
      valueBoxOutput("load_data_number_of_samples"),
      valueBoxOutput("load_data_number_of_clusters")
    ),
    fluidRow(
      valueBoxOutput("load_data_experiment_name"),
      valueBoxOutput("load_data_organism")
    )
  )