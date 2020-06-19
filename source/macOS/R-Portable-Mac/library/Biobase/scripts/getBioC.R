# This function gets and installs the required Bioconductor libraries.
#
# libName: a vector of character string for the name of the library to be
# installed. Valid names include "all" - all the released packages,
# "affy" - packages "affy" plus exprs, "CDNA" - packages "CDNA" plus
# exprs, and "exprs" - packages "Biobase", "annotate", "genefilter",
# "geneploter", "edd", "Roc", and "tkWidgets".
# destdir: a character string for the directory where the downloaded
# packages will be stored.
# isDevel: a boolean indicating whether the released (FALSE) or
# developer (TRUE) version will be downloaded and installed.
# verbose: a boolean indicating whether any error related to the
# downloading process will be (TRUE) printed. Error messages will
# still be returned but invisible if berbose is set to FALSE.
# bundle: a boolean indicating whether packages will be downloaded as
# bundles (TRUE) or individual packages (FALSE). Valid bundle
# (e. g. "all", "exprs") or package ("tkWidgets", "annotate") have to
# be used in each case.
#
getBioC <- function (libName = "exprs", destdir = NULL, isDevel = FALSE,
                     verbose = TRUE, bundle = TRUE){

    on.exit(options(show.error.messages = TRUE))

    PLATFORM <- .Platform$OS.type
    DESTDIR <- ifelse(is.null(destdir), getwd(), destdir)
    messages <- NULL
    packs <- NULL

    if(bundle){
        for(i in libName){
            packs <- c(packs, getPackNames(i))
        }
    }else{
        packs <- libName
    }

    repository <- getPkgDisc(isDevel)

    for(i in packs){
        sourceUrl <- getDLURL(i, repository, PLATFORM)
        fileName <- getFileName(sourceUrl)
        # check the connection instead of downloading directly which
        # will write files of 0 size in the directory even when
        # the connection is not there.
        options(show.error.messages = FALSE)
        tryMe <- try(url(sourceUrl, "r"))
        options(show.error.messages = TRUE)

        if(inherits(tryMe, "try-error")){
           messages <- c(messages, paste("Get", i, "failed"))
        }else{
            close(tryMe)
            download.file(sourceUrl, fileName,
                         mode = getMode(PLATFORM), quiet = TRUE)
            options(show.error.messages = FALSE)
            tryMe <- try(installPack(PLATFORM, fileName))
            options(show.error.messages = TRUE)

            if(inherits(tryMe, "try-error")){
                messages <- c(messages,
                                paste("Install", i, "failed"))
            }
        }
    }
    if(is.null(messages))
        messages <- "Download was successful"
    if(verbose)
        print(messages)
    return(invisible(messages))
}

getPackNames <- function (libName){
    error <- paste("The library is not valid. Must be:",
                      "all, exprs, affy, or CDNA", sep = "\n")
    AFFY <- "affy"
    CDNA <- c("marrayInput", "marrayClasses", "marrayNorm",
           "marrayPlots")
    EXPRS <-c("Biobase", "annotate", "genefilter", "geneplotter",
              "edd", "ROC", "tkWidgets")
    switch(libName,
           "all" = return(c(EXPRS, AFFY, CDNA)),
           "exprs" = return(EXPRS),
           "affy" = return(c(EXPRS, AFFY)),
           "cdna" =,
           "CDNA" = return(c(EXPRS, CDNA)),
           stop(error))
}

getMode <- function(platform){
    switch(platform,
           "unix" = return("w"),
           "windows" = return("wb"),
           stop("OS system not surported"))
}

installPack <- function(platform, fileName){
    if(platform == "unix"){
        system(paste("R CMD INSTALL ", fileName, sep = ""), TRUE)
    }else{
        if(platform == "windows"){
            install.packages(fileName, .libPaths()[1], CRAN = NULL)
        }else{
            stop("The OS system is not supported")
        }
    }
}

getDLURL <- function(pakName, rep, platform){
    sourceURL <- NULL
    isPkg <- FALSE
    version <- 0
    higherV <- FALSE

    for(i in rep){
        if(gsub("^Package: *(.*)", "\\1", i) == pakName){
            isPkg <- TRUE
        }
        if(isPkg && regexpr("^Version:", i)[1] > 0){
            if(gsub("^Version: *(.*)", "\\1", i) > version)
                higherV <- TRUE
            else
                higherV <- FALSE
        }
        if(platform == "windows"){
            if(isPkg && higherV && regexpr("Win32URL", i)[1] > 0){
                sourceURL <- gsub("^Win32URL: *(.*)", "\\1", i)
                isPkg = FALSE
            }
        }else{
            if(isPkg && higherV && regexpr("/Source/", i)[1] > 0){
                sourceURL <- gsub("^SourceURL: *(.*)", "\\1", i)
                isPkg <- FALSE
            }
        }
    }
    return(sourceURL)
}

getPkgDisc <- function (isDevel){
    on.exit(options(show.error.messages = TRUE))

    if(isDevel)
        URL <-
            "http://www.bioconductor.org/packages/devel/distrib/PACKAGES"
    else
        URL <-
            "http://www.bioconductor.org/packages/release/distrib/PACKAGES"

    con <- url(URL)
    options(show.error.messages = FALSE)
    tryMe <- try(readLines(con))
    options(show.error.messages = TRUE)

    if(inherits(tryMe, "try-error"))
       stop("The url for BioC PACKAGES is incorrect")

    close(con)
    return(tryMe)
}

getFileName <- function(url){
    temp <- unlist(strsplit(url, "/"))
    return(temp[length(temp)])
}













