##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
server <- function(input, output, session) {

  ##--------------------------------------------------------------------------##
  ## Color management.
  ##--------------------------------------------------------------------------##
  # Dutch palette from flatuicolors.com
  colors_dutch <- c(
    "#FFC312","#C4E538","#12CBC4","#FDA7DF","#ED4C67",
    "#F79F1F","#A3CB38","#1289A7","#D980FA","#B53471",
    "#EE5A24","#009432","#0652DD","#9980FA","#833471",
    "#EA2027","#006266","#1B1464","#5758BB","#6F1E51"
  )

  # Spanish palette from flatuicolors.com
  colors_spanish <- c(
    "#40407a","#706fd3","#f7f1e3","#34ace0","#33d9b2",
    "#2c2c54","#474787","#aaa69d","#227093","#218c74",
    "#ff5252","#ff793f","#d1ccc0","#ffb142","#ffda79",
    "#b33939","#cd6133","#84817a","#cc8e35","#ccae62"
  )

  colors <- c(colors_dutch, colors_spanish)

  cell_cycle_colorset <- setNames(
    c("#45aaf2", "#f1c40f", "#e74c3c", "#7f8c8d"),
    c("G1",      "S",       "G2M",     "-")
  )

  ##--------------------------------------------------------------------------##
  ## Central parameters.
  ##--------------------------------------------------------------------------##
  scatter_plot_dot_size <- list(
    min = 1,
    max = 20,
    step = 1,
    default = 5
  )

  scatter_plot_dot_opacity <- list(
    min = 0.1,
    max = 1.0,
    step = 0.1,
    default = 1.0
  )

  scatter_plot_percentage_cells_to_show <- list(
    min = 10,
    max = 100,
    step = 10,
    default = 100
  )

  ##--------------------------------------------------------------------------##
  ## Sidebar menu.
  ##--------------------------------------------------------------------------##
  output[["sidebar_menu"]] <- renderMenu({
    sidebarMenu(id = "sidebar",
      menuItem(
        "Load data", tabName = "loadData",
        icon = icon("spinner"), selected = TRUE
      ),
      menuItem(
        "Overview", tabName = "overview",
        icon = icon("binoculars")
      ),
      menuItem(
        "Samples", tabName = "samples",
        icon = icon("star")
      ),
      menuItem(
        "Clusters", tabName = "clusters",
        icon = icon("braille")
      ),
      menuItem(
        "Most expressed genes", tabName = "mostExpressedGenes",
        icon = icon("bullhorn")
      ),
      menuItem(
        "Marker genes", tabName = "markerGenes",
        icon = icon("magnet")
      ),
      menuItem(
        "Enriched pathways", tabName = "enrichedPathways",
        icon = icon("sitemap")
      ),
      menuItem(
        "Gene expression", tabName = "geneExpression",
        icon = icon("signal")
      ),
      menuItem(
        "Gene set expression", tabName = "geneSetExpression",
        icon = icon("list")
      ),
      menuItem(
        "Gene ID conversion", tabName = "geneIdConversion",
        icon = icon("barcode")
      ),
      menuItem(
        "Analysis info", tabName = "info",
        icon = icon("info")
      ),
      menuItem(
        "About", tabName = "about",
        icon = icon("at")
      )
    )
  })

  ##--------------------------------------------------------------------------##
  ## Sample data.
  ##--------------------------------------------------------------------------##
  sample_data <- reactive({
    if ( is.null(input$RDS_file) || is.na(input$RDS_file) ) {
      sample_data <- readRDS("resources/example.rds")
    } else {
      req(input$RDS_file)
      sample_data <- readRDS(input$RDS_file$datapath)
    }
    sample_data$sample_names <- levels(sample_data$cells$sample)
    sample_data$cluster_names <- levels(sample_data$cells$cluster)
    sample_data
  })

  ##--------------------------------------------------------------------------##
  ## Tabs.
  ##--------------------------------------------------------------------------##
  source("shiny/load_data/server.R", local = TRUE)
  source("shiny/overview/server.R", local = TRUE)
  source("shiny/samples/server.R", local = TRUE)
  source("shiny/clusters/server.R", local = TRUE)
  source("shiny/most_expressed_genes/server.R", local = TRUE)
  source("shiny/marker_genes/server.R", local = TRUE)
  source("shiny/enriched_pathways/server.R", local = TRUE)
  source("shiny/gene_expression/server.R", local = TRUE)
  source("shiny/gene_set_expression/server.R", local = TRUE)
  source("shiny/gene_id_conversion/server.R", local = TRUE)
  source("shiny/analysis_info/server.R", local = TRUE)
  source("shiny/about/server.R", local = TRUE)

}