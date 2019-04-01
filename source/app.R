##----------------------------------------------------------------------------##
## Cerebro
## version 1.0.0
##
## Author:    Roman Hillje
## Institute: IEO
## Lab:       PGP
## Date:      2019-04-01
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
if (grepl(tolower(Sys.info()["sysname"]), pattern="^win")) {
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
if ( !is.element(el = "cerebroApp", set = rownames(installed.packages())) ) {
  BiocManager::install(
    "romanhaa/cerebroApp"
  )
}

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
cerebroApp::launchApp()
