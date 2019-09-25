### R code from vignette source 'esApply.Rnw'

###################################################
### code chunk number 1: R.hide
###################################################
library(Biobase)
data(sample.ExpressionSet)


###################################################
### code chunk number 2: R
###################################################
print(sample.ExpressionSet)
print(exprs(sample.ExpressionSet)[1,])
print(pData(sample.ExpressionSet)[1:2,1:3])


###################################################
### code chunk number 3: R
###################################################
print(rbind(exprs(sample.ExpressionSet[1,]),
sex <- t(pData(sample.ExpressionSet))[1,]))


###################################################
### code chunk number 4: R
###################################################
medContr <- function( y, x ) {
 ys <- split(y,x)
 median(ys[[1]]) - median(ys[[2]])
}


###################################################
### code chunk number 5: R
###################################################
print(apply(exprs(sample.ExpressionSet[1,,drop=F]), 1,
  medContr, pData(sample.ExpressionSet)[["sex"]]))


###################################################
### code chunk number 6: R
###################################################
medContr1 <- function(y) {
 ys <- split(y,sex)
 median(ys[[1]]) - median(ys[[2]])
}

print(esApply( sample.ExpressionSet, 1, medContr1)[1])


###################################################
### code chunk number 7: esApply.Rnw:126-127
###################################################
sessionInfo()


