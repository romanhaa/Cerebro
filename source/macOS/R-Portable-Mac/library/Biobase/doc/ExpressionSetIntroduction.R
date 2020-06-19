### R code from vignette source 'ExpressionSetIntroduction.Rnw'

###################################################
### code chunk number 1: init
###################################################
options(width=65)


###################################################
### code chunk number 2: install-pkg (eval = FALSE)
###################################################
## if (!require("BiocManager"))
##     install.packages("BiocManager")
## BiocManager::install("Biobase")


###################################################
### code chunk number 3: loadlib
###################################################
library("Biobase")


###################################################
### code chunk number 4: convert (eval = FALSE)
###################################################
## library(convert)
## as(object, "ExpressionSet")


###################################################
### code chunk number 5: read-table-geneData
###################################################
dataDirectory <- system.file("extdata", package="Biobase")
exprsFile <- file.path(dataDirectory, "exprsData.txt")
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t",
                              row.names=1,
                              as.is=TRUE))


###################################################
### code chunk number 6: exprsFile (eval = FALSE)
###################################################
## exprsFile <- "c:/path/to/exprsData.txt"


###################################################
### code chunk number 7: geneData-peak
###################################################
class(exprs)
dim(exprs)
colnames(exprs)
head(exprs[,1:5])


###################################################
### code chunk number 8: ExpressionSet-basic
###################################################
minimalSet <- ExpressionSet(assayData=exprs)


###################################################
### code chunk number 9: pData
###################################################
pDataFile <- file.path(dataDirectory, "pData.txt")
pData <- read.table(pDataFile,
                    row.names=1, header=TRUE, sep="\t")
dim(pData)
rownames(pData)
summary(pData)


###################################################
### code chunk number 10: geneCovariate-geneData-name-match
###################################################
all(rownames(pData)==colnames(exprs))


###################################################
### code chunk number 11: colnames
###################################################
names(pData)


###################################################
### code chunk number 12: sapplyClasses
###################################################
sapply(pData, class)


###################################################
### code chunk number 13: simpleSubsetting
###################################################
pData[c(15, 20), c("gender", "type")]
pData[pData$score>0.8,]


###################################################
### code chunk number 14: metadata-create
###################################################
metadata <- data.frame(labelDescription=
                       c("Patient gender", 
                         "Case/control status", 
                         "Tumor progress on XYZ scale"),
                       row.names=c("gender", "type", "score"))


###################################################
### code chunk number 15: AnnotatedDataFrame
###################################################
phenoData <- new("AnnotatedDataFrame", 
                 data=pData, varMetadata=metadata)
phenoData


###################################################
### code chunk number 16: AnnotatedDataFrame-subset
###################################################
head(pData(phenoData))
phenoData[c("A","Z"),"gender"]
pData(phenoData[phenoData$score>0.8,])


###################################################
### code chunk number 17: annotation
###################################################
annotation <- "hgu95av2"


###################################################
### code chunk number 18: R.MIAME
###################################################
experimentData <- new("MIAME",
  name="Pierre Fermat",
  lab="Francis Galton Lab",
  contact="pfermat@lab.not.exist",
  title="Smoking-Cancer Experiment",
  abstract="An example ExpressionSet",
  url="www.lab.not.exist",
  other=list(
    notes="Created from text files"
  ))


###################################################
### code chunk number 19: ExpressionSetFinally
###################################################
exampleSet <- ExpressionSet(assayData=exprs, 
                  phenoData=phenoData, 
                  experimentData=experimentData,
                  annotation="hgu95av2")


###################################################
### code chunk number 20: ExpressionSet-minimal
###################################################
minimalSet <- ExpressionSet(assayData=exprs)


###################################################
### code chunk number 21: helpExpressionSet (eval = FALSE)
###################################################
## help("ExpressionSet-class")


###################################################
### code chunk number 22: showExpressionSet
###################################################
exampleSet


###################################################
### code chunk number 23: usingDollar
###################################################
exampleSet$gender[1:5]
exampleSet$gender[1:5] == "Female"


###################################################
### code chunk number 24: featureNames
###################################################
featureNames(exampleSet)[1:5]


###################################################
### code chunk number 25: sampleNames
###################################################
sampleNames(exampleSet)[1:5]
varLabels(exampleSet)


###################################################
### code chunk number 26: exprs
###################################################
mat <- exprs(exampleSet)
dim(mat)


###################################################
### code chunk number 27: first10
###################################################
vv <- exampleSet[1:5, 1:3]
dim(vv)
featureNames(vv)
sampleNames(vv)


###################################################
### code chunk number 28: males
###################################################
males <- exampleSet[ , exampleSet$gender == "Male"]
males


###################################################
### code chunk number 29: ExpressionSetIntroduction.Rnw:490-491
###################################################
toLatex(sessionInfo())


