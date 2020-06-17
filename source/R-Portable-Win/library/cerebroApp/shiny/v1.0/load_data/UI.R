##----------------------------------------------------------------------------##
## Tab: Load data.
##----------------------------------------------------------------------------##

tab_load_data <- tabItem(
  tabName = "loadData",
  fluidRow(
    column(12,
      titlePanel("Load data"),
      fileInput(
        inputId = "RDS_file",
        label = "Choose RDS file...",
        multiple = FALSE,
        accept = c(".rds",".crb",".cerebro"),
        width = NULL,
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
