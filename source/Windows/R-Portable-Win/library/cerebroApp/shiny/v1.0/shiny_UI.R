##----------------------------------------------------------------------------##
## Custom functions.
##----------------------------------------------------------------------------##
cerebroBox <- function(title, content) {
  box(
    title = title,
    status = "primary",
    solidHeader = TRUE,
    width = 12,
    collapsible = TRUE,
    content
  )
}

cerebroInfoButton <- function(id) {
  actionButton(
    inputId = id,
    label = "info",
    icon = NULL,
    class = "btn-xs",
    title = "Show additional information for this panel."
  )
}

boxTitle <- function(title) {
  p(title, style = "padding-right: 5px; display: inline")
}

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
source(system.file("shiny/v1.0/load_data/UI.R", package = "cerebroApp"), local = TRUE)
source(system.file("shiny/v1.0/overview/UI.R", package = "cerebroApp"), local = TRUE)
source(system.file("shiny/v1.0/samples/UI.R", package = "cerebroApp"), local = TRUE)
source(system.file("shiny/v1.0/clusters/UI.R", package = "cerebroApp"), local = TRUE)
source(system.file("shiny/v1.0/most_expressed_genes/UI.R", package = "cerebroApp"), local = TRUE)
source(system.file("shiny/v1.0/marker_genes/UI.R", package = "cerebroApp"), local = TRUE)
source(system.file("shiny/v1.0/enriched_pathways/UI.R", package = "cerebroApp"), local = TRUE)
source(system.file("shiny/v1.0/gene_expression/UI.R", package = "cerebroApp"), local = TRUE)
source(system.file("shiny/v1.0/gene_set_expression/UI.R", package = "cerebroApp"), local = TRUE)
source(system.file("shiny/v1.0/gene_id_conversion/UI.R", package = "cerebroApp"), local = TRUE)
source(system.file("shiny/v1.0/analysis_info/UI.R", package = "cerebroApp"), local = TRUE)
source(system.file("shiny/v1.0/about/UI.R", package = "cerebroApp"), local = TRUE)

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
ui <- dashboardPage(
  dashboardHeader(
    title = span("Cerebro ", style = "color: white; font-size: 28px; font-weight: bold")
  ),
  dashboardSidebar(
    tags$head(tags$style(HTML(".content-wrapper {overflow-x: scroll;}"))),
    sidebarMenu(
      sidebarMenuOutput("sidebar_menu")
    )
  ),
  dashboardBody(
    tags$script(HTML('$("body").addClass("fixed");')),
    tabItems(
      tab_load_data,
      tab_overview,
      tab_samples,
      tab_clusters,
      tab_most_expressed_genes,
      tab_marker_genes,
      tab_enriched_pathways,
      tab_gene_expression,
      tab_gene_set_expression,
      tab_gene_id_conversion,
      tab_analysis_info,
      tab_about
    )
  )
)
