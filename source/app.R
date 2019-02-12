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
library("dplyr")
library("formattable")
library("shiny")
library("shinydashboard")

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
## Color management.
##----------------------------------------------------------------------------##
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

##----------------------------------------------------------------------------##
## App.
##----------------------------------------------------------------------------##
source("shiny/shiny_UI.R", local = TRUE)
source("shiny/shiny_server.R")

shinyApp(ui, server)


