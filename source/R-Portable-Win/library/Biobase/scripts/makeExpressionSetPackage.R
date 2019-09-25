 ##take a phenoData object and create a valid - sort of - format
 ## section for a man page
 pD2Rd <- function(pD) {
   if(!inherits(pD, "AnnotatedDataFrame") )
     stop("only works for AnnotatedDataFrames")

   fmt = "\\format{\n  The format is:\n  An \\code{ExpressionSetObject} with covariates:\n"
   covs = "\\itemize{"
   vMD = varMetadata(pD)
   vL = varLabels(pD)
   for(i in 1:length(vL) ) {
     item = paste("\\item \\code{", vL[i], "}: ", vMD[i,1], sep="")
     covs = paste(covs, item, sep="\n")
   }
   paste(fmt, covs, "\n}\n}\n", sep="")     
 }


 makeExpressionSetPackage = function(expS, author, filePath=tempdir(),
    version = "1.0.0", license, email, biocViews="ExperimentData",
    packageName  )
 {
   if( !inherits(expS, "ExpressionSet") )
     stop("only works for ExpressionSets")
 
   if( missing(email) || !(is.character(email) && (length(email) == 1)
           && grep("@", email) == 1 ) )
     stop("invalid email address")

   if( !is.package_version(version) )
     version = package_version(version)

   if(missing(license) ) 
      license= "The Artistic License, Version 2.0"

   if( missing(packageName) )
       pkgname = deparse(substitute(expS))

   sym = list(AUTHOR = author, VERSION=as.character(version), LICENSE=license,
        TITLE = paste("Experimental Data Package:",pkgname), 
        MAINTAINER = paste(author, ", <", email, ">", sep = ""),
        BVIEWS = biocViews, DESCRIPTION = "place holder 1",
        FORMAT = "An instance of the ExpressionSet class")

   
   res = createPackage(pkgname, destinationDir=filePath,
         originDir = system.file("ExpressionSet", package="Biobase"),
         symbolValues = sym, unlink=TRUE)

   ##save the data file
   assign(pkgname, expS)
   save(list=pkgname, file = file.path( res$pkgdir, "data", 
                          paste(pkgname, ".rda", sep=""))) 
   return(res)
 }


library(ALL)
data(ALL)
makeExpressionSetPackage(ALL, author="Robert Gentleman", email="rgentlem@foo")

 x=pD2Rd(phenoData(ALL))

