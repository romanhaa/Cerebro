# This function queries the repository for packages and finds the urls
# for all the packages a given package depends on. The urls will be
# returned as list that contains the repository, url, and version
# number of the package and a list of the urls for all the packages the
# package depends on.
#
# pkgName: a character string for the name of the package whose
# dependency will be queried.
# pkgVersion: a integer for the version number of the packages on
# which a given package depends. Default is to get the largest versin
# number.
# type: a character string for the type of the package. Should be
# either "unix" or "windows".
# repositories: a list of urls for the respositories to be
# searched. The default is bioconductor and then CRAN.
#
query.packages <- function (pkgName, pkgVersion = NULL, type = "unix",
                            repositories = getDefaultRep(bioCOnly = FALSE)){

    on.exit(options(show.error.messages = TRUE))

    pkgFound <- FALSE
    collectedItem <- NULL
    pkgPat <- getPkgPat(pkgName)
    returnList <- list(repository = NULL,
                       packUrl = NULL,
                       version = ifelse(is.null(pkgVersion), -1, pkgVersion),
                       type = ifelse(is.null(type), .Platform$OS.type, type),
                       depends = NULL)
    tempRep <- NULL
    tempSUrl <- NULL
    tempZUrl <- NULL
    tempDep <- NULL
    tempVer <- -1
    listReady <- FALSE

    doList <- function(){

        if(is.null(pkgVersion)){
            if(tempRep == paste(getOption("repositories"),
               "/PACKAGES", sep = "")){
                tempSUrl <<- paste(getOption("repositories"), "/",
                                   pkgName, "_", tempVer, ".tar.gz",
                                   sep = "")
                tempZUrl <<- avlCheck(pkgName)
            }
            if(tempVer >= returnList$version){
                returnList$repository <<- tempRep
                if(!is.null(tempSUrl) && !is.null(tempZUrl))
                    returnList$packUrl <<- ifelse(type == "unix",
                                                  tempSUrl, tempZUrl)
                returnList$version <<- tempVer
                if(!is.null(tempDep))
                    returnList$depends <<- formatLine(tempDep)
                listReady <<- TRUE
                pkgFound <<- FALSE
            }
        }else{
            if(tempVer == returnList$version){
                returnList$repository <<- tempRep
                if(!is.null(tempSUrl) && !is.null(tempZUrl))
                    returnList$packUrl <<- ifelse(type == "unix",
                                                  tempSUrl, tempZUrl)
                returnList$version <<- tempVer
                if(!is.null(tempDep))
                    returnList$depends <<- formatLine(tempDep)
                listReady <<- TRUE
                pagFound <<- FALSE
            }
        }
    }

    for(i in repositories){
        pakRep <- getRep(i)
        tempRep <- i
        for(j in pakRep){
            if(isPak(j)){
                # doList whenever we see a match
                if(pkgFound)
                    doList()
                if(regexpr(pkgPat, j) == 1)
                    pkgFound <- TRUE
                else
                    pkgFound <- FALSE
            }
            if(pkgFound){
                switch(sub("(^.*): .*", "\\1", j),
                       "Version" = tempVer <- sub("^.*: *(.*)","\\1",j),
                       "SourceURL" = tempSUrl <-
                       sub("^.*(http://.*)","\\1",j),
                       "Win32URL" = tempZUrl <- sub("^.*(http://.*)","\\1",j),
                       "Depends" = tempDep <- getDepends(j))
            }
        }
        # do list again in case the match is the last one
        if(pkgFound)
            doList()
         # Taken out to be in line with the new format of PACKAGES
#        if(listReady)
#            break
    }

    if(listReady)
        return(returnList)
    else{
        toPut <- ifelse(is.null(pkgVersion), "",
                           paste(" Version", pkgVersion))
        print(paste("Can not find package ", pkgName,
                    toPut, " in repositories", sep = ""))
        return(NULL)
    }
}

getDefaultRep <- function (bioCOnly = FALSE){
    if(bioCOnly)
        return(list(BioC =
           "http://www.bioconductor.org/packages/release/distrib/PACKAGES"))
    else
        return(list(BioC =
           "http://www.bioconductor.org/packages/release/distrib/PACKAGES",
           CRAN = paste(getOption("repositories"), "/PACKAGES", sep = "")))
}

getRep <- function(rep){
    on.exit(options(show.error.messages = TRUE))

    con <- url(rep)
    options(show.error.messages = FALSE)
    tryMe <- try(readLines(con))
    #for testing using a local file only
    #tryMe <- try(readLines("PACKAGES"))
    options(show.error.messages = TRUE)

    if(inherits(tryMe, "try-error"))
       stop(paste("Invalid repository url", rep))

    close(con)
    return(tryMe)
}

getPkgPat <- function(pkgName){
    return(paste("Package: *", pkgName, sep = ""))
}

isPak <- function(aLine){
    return(ifelse(regexpr("Package:.*", aLine) == 1, TRUE, FALSE))
}

getDepends <- function(dep){
    on.exit(options(show.error.messages = TRUE))
    depList <- NULL

    options(show.error.messages = FALSE)
    tryMe <- try(sub("^.*: *(.*)","\\1",dep))
    options(show.error.messages = TRUE)

    if(tryMe == "" || is.null(tryMe) || is.na(tryMe)){
        depList <- NULL
    }else{
        dep <- sub("^.*: *(.*)","\\1",dep)
        deps <- unlist(strsplit(dep, ","))
        for(i in 1:length(deps)){
            depList[[i]] <- deps[i]
        }
    }
    return(depList)
}

formatLine <- function(aLine){
    return(gsub("^ *(.*)", "\\1", aLine))
}

avlCheck <- function(pkgName){
    on.exit(options(show.error.message = TRUE))

    options(show.error.messages = FALSE)
    tryMe <- try(url(paste("http://cran.r-project.org/bin/",
                                   "windows/contrib/", pkgName,
                                   ".zip", sep = ""), "r"))
    if(!inherits(tryMe, "try-error")){
        close(tryMe)
        return(paste("http://cran.r-project.org/bin/",
                                   "windows/contrib/", pkgName,
                                   ".zip", sep = ""))
    }
    tryMe <- try(url(paste("http://cran.r-project.org/bin/",
                                   "windows/base/", pkgName,
                                   ".zip", sep = ""), "r"))
    if(!inherits(tryMe, "try-error")){
        close(tryMe)
        return(paste("http://cran.r-project.org/bin/",
                                   "windows/base/", pkgName,
                                   ".zip", sep = ""))
    }else{
        return("Unavailable")
    }
}













