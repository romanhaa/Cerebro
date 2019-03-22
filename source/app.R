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
required_packages <- c(
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
  "shiny",
  "shinydashboard",
  "shinyWidgets"
)

packages_not_present <- required_packages[which(required_packages %in% rownames(installed.packages()) == FALSE)]

BiocManager::install(
  packages_not_present,
  type = "binary",
  update = FALSE
)

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
require("dplyr")
require("formattable")
require("shiny")
require("shinydashboard")

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
source("shiny/overview/info.R")
source("shiny/samples/info.R")
source("shiny/clusters/info.R")
source("shiny/most_expressed_genes/info.R")
source("shiny/marker_genes/info.R")
source("shiny/enriched_pathways/info.R")
source("shiny/gene_expression/info.R")
source("shiny/gene_set_expression/info.R")

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
# system("type R")
# print(sessionInfo())
# Sys.setenv("R_LIBS_USER" = "")
# print(Sys.getenv())

##----------------------------------------------------------------------------##
## Allow upload of files up to 800 MB.
##----------------------------------------------------------------------------##
options(shiny.maxRequestSize = 800*1024^2) 

##----------------------------------------------------------------------------##
## App.
##----------------------------------------------------------------------------##
source("shiny/shiny_UI.R", local = TRUE)
source("shiny/shiny_server.R", local = TRUE)

shinyApp(ui, server)


