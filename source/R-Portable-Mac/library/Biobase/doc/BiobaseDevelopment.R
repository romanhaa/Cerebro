### R code from vignette source 'BiobaseDevelopment.Rnw'

###################################################
### code chunk number 1: Biobase
###################################################
options(width=69)
library(Biobase)


###################################################
### code chunk number 2: eSet-class
###################################################
getClass("eSet")


###################################################
### code chunk number 3: eSet-validity
###################################################
getValidity(getClass("eSet"))


###################################################
### code chunk number 4: ExpressionSet-initialize (eval = FALSE)
###################################################
## obj <- new("ExpressionSet", phenoData = new("AnnotatedDataFrame"), experimentData = new("MIAME"), annotation = character(), exprs = new("matrix")) 


###################################################
### code chunk number 5: newAssayData (eval = FALSE)
###################################################
## assayDataNew("environment", elt)


###################################################
### code chunk number 6: assayData-storageMode
###################################################
data(sample.ExpressionSet)
storageMode(sample.ExpressionSet)
tryCatch(assayData(sample.ExpressionSet)$exprs <- log(exprs(sample.ExpressionSet)),
         error=function(err) cat(conditionMessage(err)))
exprs(sample.ExpressionSet) <- log(exprs(sample.ExpressionSet))


###################################################
### code chunk number 7: ExpressionSet-class
###################################################
getClass("ExpressionSet")
getValidity(getClass("ExpressionSet"))


###################################################
### code chunk number 8: SwirlSet-class
###################################################
setClass("SwirlSet", contains="eSet")


###################################################
### code chunk number 9: SwirlSet-initialize
###################################################
setMethod("initialize", "SwirlSet",
          function(.Object,
                   R = new("matrix"),
                   G = new("matrix"),
                   Rb = new("matrix"),
                   Gb = new("matrix"),
                   ...) {
            callNextMethod(.Object,
                           R=R, G=G, Rb=Rb, Gb=Gb,
                           ...)
        })


###################################################
### code chunk number 10: SwirlSet-initialize-2
###################################################
setMethod("initialize", "SwirlSet",
          function(.Object,
                   assayData=assayDataNew(
                     R=R, G=G, Rb=Rb, Gb=Gb),
                   R = new("matrix"),
                   G = new("matrix"),
                   Rb = new("matrix"),
                   Gb = new("matrix"),
                   ...) {
            if (!missing(assayData) && 
                any(!missing(R), !missing(G), !missing(Rb), !missing(Gb))) {
                warning("using 'assayData'; ignoring 'R', 'G', 'Rb', 'Gb'")
            }
            callNextMethod(.Object, assayData=assayData, ...)
        })


###################################################
### code chunk number 11: SwirlSet-new
###################################################
new("SwirlSet")


###################################################
### code chunk number 12: initialize-.Object (eval = FALSE)
###################################################
## setMethod("initialize", "MySet",
##           function(.Object, ...) {
##               .Object <- callNextMethod(.Object, ...)
##           })
##               .


###################################################
### code chunk number 13: SwirlSet-validity
###################################################
setValidity("SwirlSet", function(object) {
  assayDataValidMembers(assayData(object), c("R", "G", "Rb", "Gb"))
})


###################################################
### code chunk number 14: validity-sometimes (eval = FALSE)
###################################################
## myFancyFunction <- function(obj) {
##   assayData(obj) <- fancyAssaydData # obj invalid...
##   phenoData(obj) <- justAsFancyPhenoData # but now valid
##   validObject(obj)
##   obj
## }


###################################################
### code chunk number 15: updateObject-eg
###################################################
data(sample.ExpressionSet)
classVersion(sample.ExpressionSet)
obj <- updateObject(sample.ExpressionSet)


###################################################
### code chunk number 16: isCurrent
###################################################
isCurrent(sample.ExpressionSet)[c("eSet", "ExpressionSet")]


###################################################
### code chunk number 17: MultiSet-obj
###################################################
setClass("MySet",
         contains = "eSet",
         prototype = prototype(
           new("VersionedBiobase",
               versions=c(classVersion("eSet"), MySet="1.0.0"))))
obj <- new("MySet")
classVersion(obj)


###################################################
### code chunk number 18: MultiSetRevised
###################################################
setClass("MySet",
         contains = "eSet",
         prototype = prototype(
           new("VersionedBiobase",
               versions=c(classVersion("eSet"), MySet="1.0.1"))))
isCurrent(obj)


###################################################
### code chunk number 19: updateObject-MultiSet
###################################################
setMethod("updateObject", signature(object="MySet"),
          function(object, ..., verbose=FALSE) {
              if (verbose) message("updateObject(object = 'MySet')")
              object <- callNextMethod()
              if (isCurrent(object)["MySet"]) return(object)
              ## Create an updated instance.
              if (!isVersioned(object))
                  ## Radical surgery -- create a new, up-to-date instance
                  new("MySet",
                      assayData = updateObject(assayData(object),
                        ...., verbose=verbose),
                      phenoData = updateObject(phenoData(object),
                        ..., verbose=verbose),
                      experimentData = updateObject(experimentData(object),
                        ..., verbose=verbose),
                      annotation = updateObject(annotation(object),
                        ..., verbose=verbose))
              else {
                  ## Make minor changes, and update version by consulting class definition
                  classVersion(object)["MySet"] <-
                      classVersion("MySet")["MySet"]
                  object
              }
          })


###################################################
### code chunk number 20: updateObject
###################################################
classVersion(updateObject(obj))


###################################################
### code chunk number 21: classVersion-AnnotatedDataFrame
###################################################
classVersion(new("AnnotatedDataFrame"))


###################################################
### code chunk number 22: SwirlSet-version
###################################################
setClass("SwirlSet", contains = "eSet",
         prototype = prototype(
           new("VersionedBiobase",
               versions=c(classVersion("eSet"), SwirlSet="1.0.0"))))
classVersion(new("SwirlSet"))


###################################################
### code chunk number 23: arbitraryClassVersions
###################################################
obj <- new("SwirlSet")
classVersion(obj)["MyID"] <- "0.0.1"
classVersion(obj)
classVersion(updateObject(obj))


###################################################
### code chunk number 24: BiobaseDevelopment.Rnw:680-681
###################################################
toLatex(sessionInfo())


