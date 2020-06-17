##----------------------------------------------------------------------------##
## set .libPaths() to load packages from inside Cerebro standalone version
##----------------------------------------------------------------------------##
if ( grepl(tolower(Sys.info()["sysname"]), pattern = "^win") )
{
  .libPaths(paste0(getwd(), "/R-Portable-Win/library"))
} else if ( grepl(tolower(Sys.info()["sysname"]), pattern = "darwin") )
{
  .libPaths(paste0(getwd(), "/R-Portable-Mac/library"))
}

##----------------------------------------------------------------------------##
## launch Cerebro interface
##----------------------------------------------------------------------------##
cerebroApp::launchCerebro()
