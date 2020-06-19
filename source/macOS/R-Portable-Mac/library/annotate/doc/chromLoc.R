### R code from vignette source 'chromLoc.Rnw'

###################################################
### code chunk number 1: buildCL
###################################################
library("annotate")
z <- buildChromLocation("hgu95av2")
z


###################################################
### code chunk number 2: showBasicMethods
###################################################
organism(z)

dataSource(z)

## The chromLocs list is extremely large.  Let's only
## look at one of the elements.
names(chromLocs(z))
chromLocs(z)[["Y"]]

get("32972_at", probesToChrom(z))

chromInfo(z)

get("32972_at", geneSymbols(z))



###################################################
### code chunk number 3: nChrom
###################################################
nChrom(z)


