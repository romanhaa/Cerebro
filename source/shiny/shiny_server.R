##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
server <- function(input, output, session) {

  ##--------------------------------------------------------------------------##
  ## Sidebar menu.
  ##--------------------------------------------------------------------------##
  output$sidebar_menu <- renderMenu({
    if ( mode == "open" ) {
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
    } else {
      sidebarMenu(id = "sidebar",
        menuItem(
          "Load data", tabName = "loadData", icon = icon("spinner")
        ),
        menuItem(
          "Overview", tabName = "overview", icon = icon("binoculars"),
          selected = TRUE
        ),
        menuItem(
          "Samples", tabName = "samples", icon = icon("star")
        ),
        menuItem(
          "Clusters", tabName = "clusters", icon = icon("braille")
        ),
        menuItem(
          "Most expressed genes", tabName = "mostExpressedGenes",
          icon = icon("bullhorn")
        ),
        menuItem(
          "Marker genes", tabName = "markerGenes", icon = icon("magnet")
        ),
        menuItem(
          "Enriched pathways", tabName = "enrichedPathways",
          icon = icon("sitemap")
        ),
        menuItem(
          "Gene expression", tabName = "geneExpression", icon = icon("signal")
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
          "Analysis info", tabName = "info", icon = icon("info")
        ),
        menuItem(
          "About", tabName = "about", icon = icon("at")
        )
      )
    }
  })

  ##--------------------------------------------------------------------------##
  ## Sample data.
  ##--------------------------------------------------------------------------##
  sample_data <- reactive({

    if ( mode == "boxed" ) {
      sample_data <- readRDS("resources/data.rds")
    } else if ( is.null(input$RDS_file) || is.na(input$RDS_file) ) {
      sample_data <- readRDS("resources/example.rds")
    } else {
      req(input$RDS_file)
      sample_data <- readRDS(input$RDS_file$datapath)
    }

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