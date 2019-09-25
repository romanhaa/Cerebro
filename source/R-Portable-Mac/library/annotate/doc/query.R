### R code from vignette source 'query.Rnw'

###################################################
### code chunk number 1: data
###################################################
library("annotate")
data(sample.ExpressionSet)
affys <- featureNames(sample.ExpressionSet)[490:500]
affys


###################################################
### code chunk number 2: annotation
###################################################
library("hgu95av2.db")
ids <- getPMID(affys,"hgu95av2")
ids <- unlist(ids,use.names=FALSE)
ids <- unique(ids[!is.na(as.numeric(ids))])
length(ids)
ids[1:10]


###################################################
### code chunk number 3: getabsts
###################################################
x <- pubmed(ids[1:10])
a <- xmlRoot(x)
numAbst <- length(xmlChildren(a))
numAbst


###################################################
### code chunk number 4: query.Rnw:170-179
###################################################
arts <- vector("list", length=numAbst)
absts <- rep(NA, numAbst)
for (i in 1:numAbst) {
   ## Generate the PubMedAbst object for this abstract
   arts[[i]] <- buildPubMedAbst(a[[i]])
   ## Retrieve the abstract text for this abstract
   absts[i] <- abstText(arts[[i]])
}
arts[[7]]


###################################################
### code chunk number 5: query.Rnw:203-206
###################################################
found <- grep("cDNA",absts)
goodAbsts <- arts[found]
length(goodAbsts)


###################################################
### code chunk number 6: query.Rnw:222-224
###################################################
y <- genbank(ids[1:10], type="uid")
b <- xmlRoot(y)


###################################################
### code chunk number 7: abst2HTML
###################################################
fname <- tempfile()
pmAbst2HTML(goodAbsts, filename=fname)

fnameBase <- tempfile()
pmAbst2HTML(goodAbsts, filename=fnameBase, frames=TRUE)


###################################################
### code chunk number 8: query.Rnw:267-268
###################################################
sessionInfo()


