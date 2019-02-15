##----------------------------------------------------------------------------##
## Cerebro
## version 1.0
##
## Author:    Roman Hillje
## Institute: IEO
## Lab:       PGP
## Date:      
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
if (grepl(tolower(Sys.info()['sysname']), pattern='^win')) {
    
  .libPaths(paste0(getwd(), "/R-Portable-Win/library"))

  plot_export_path <- paste0(Sys.getenv("USERPROFILE"), "\\Desktop\\")

} else {

  .libPaths(paste0(getwd(), "/R-Portable-Mac/library"))

  plot_export_path <- "~/Desktop/"

}

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
if ( !is.element(el = "BiocManager", set = rownames(installed.packages())) ) {
  install.packages(
    "BiocManager",
    repos = "http://cran.us.r-project.org",
    dependencies = TRUE
  )
}

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
required_packages_CRAN <- c(
    "DT",
    "formattable",
    "ggplot2",
    "ggtree",
    "Matrix",
    "msigdbr",
    "plotly",
    "RColorBrewer",
    "reshape2",
    "scales",
    "scatterD3",
    "shiny",
    "shinydashboard",
    "shinyWidgets"
  )

for ( package in required_packages_CRAN ) {
  if ( !is.element(el = package, set = rownames(installed.packages())) ) {
    BiocManager::install(package)
  }
}

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
require("dplyr")
require("formattable")
require("shiny")
require("shinydashboard")
#require("highcharter")

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
source("shiny/overview/info.txt")
source("shiny/samples/info.txt")
source("shiny/clusters/info.txt")
source("shiny/most_expressed_genes/info.txt")
source("shiny/marker_genes/info.txt")
source("shiny/enriched_pathways/info.txt")
source("shiny/gene_expression/info.txt")
source("shiny/gene_set_expression/info.txt")

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
# system("type R")
# print(sessionInfo())
# Sys.setenv("R_LIBS_USER" = "")
# print(Sys.getenv())

##----------------------------------------------------------------------------##
## Allow upload of files up to 400 MB.
##----------------------------------------------------------------------------##
options(shiny.maxRequestSize = 400*1024^2) 

##----------------------------------------------------------------------------##
## App.
##----------------------------------------------------------------------------##
source("shiny/shiny_UI.R", local = TRUE)
source("shiny/shiny_server.R", local = TRUE)

shinyApp(ui, server)


