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
source("shiny/load_data/UI.R", local = TRUE)
source("shiny/overview/UI.R", local = TRUE)
source("shiny/samples/UI.R", local = TRUE)
source("shiny/clusters/UI.R", local = TRUE)
source("shiny/most_expressed_genes/UI.R", local = TRUE)
source("shiny/marker_genes/UI.R", local = TRUE)
source("shiny/enriched_pathways/UI.R", local = TRUE)
source("shiny/gene_expression/UI.R", local = TRUE)
source("shiny/gene_set_expression/UI.R", local = TRUE)
source("shiny/gene_id_conversion/UI.R", local = TRUE)
source("shiny/analysis_info/UI.R", local = TRUE)
source("shiny/about/UI.R", local = TRUE)

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

