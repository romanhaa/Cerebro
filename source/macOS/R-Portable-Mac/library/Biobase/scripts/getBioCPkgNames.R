# This function provides a widget for users to select packages to
# download from Bioconductor.
#
getBioCWidget <- function(bundle = TRUE){

    nameList <- list()
    popList <- function(toPut){
        nameList[[toPut]] <<- "BioC Package"
    }

    if(bundle)
        pkgNames <- getBioCBundle()
    else
        pkgNames <- getBioCPkgNames()
    trash <- sapply(pkgNames, popList)

    selected <- listSelect(nameList, "Select packages from the list",
                           NULL, NULL)
    toGet <- NULL
    for(i in names(selected)){
        if(selected[[i]]){
            toGet <- c(toGet, i)
        }
    }
    getBioC(libName = toGet, destdir = NULL, isDevel = FALSE,
                     verbose = TRUE, bundle = bundle)
}

getBioCPkgNames <- function(){

    pkgNames <- NULL
    biocURL <- getDefaultRep(bioCOnly = TRUE)
    repository <- getRep(biocURL[[1]])

    for(i in repository){
        if(isPak(i)){
            pkgNames <- c(pkgNames, gsub("^Package: *(.*)", "\\1", i))
        }
    }
    return(unique(pkgNames))
}

getBioCBundle <- function(){
    return(c("exprs", "affy", "cdna"))
}







