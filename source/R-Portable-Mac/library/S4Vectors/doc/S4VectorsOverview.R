### R code from vignette source 'S4VectorsOverview.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex(use.unsrturl=FALSE)


###################################################
### code chunk number 2: options
###################################################
options(width=72)


###################################################
### code chunk number 3: install (eval = FALSE)
###################################################
## if (!require("BiocManager"))
##     install.packages("BiocManager")
## BiocManager::install("S4Vectors")


###################################################
### code chunk number 4: initialize
###################################################
library(S4Vectors)


###################################################
### code chunk number 5: Rle-extends-Vector
###################################################
showClass("Rle")


###################################################
### code chunk number 6: initialize
###################################################
set.seed(0)
lambda <- c(rep(0.001, 4500), seq(0.001, 10, length=500),
            seq(10, 0.001, length=500))
xVector <- rpois(1e7, lambda)
yVector <- rpois(1e7, lambda[c(251:length(lambda), 1:250)])
xRle <- Rle(xVector)
yRle <- Rle(yVector)


###################################################
### code chunk number 7: basic-ops
###################################################
length(xRle)
xRle[1]
zRle <- c(xRle, yRle)


###################################################
### code chunk number 8: seq-extraction
###################################################
xSnippet <- window(xRle, 4751, 4760)
xSnippet
head(xSnippet)
tail(xSnippet)
rev(xSnippet)
rep(xSnippet, 2)
subset(xSnippet, xSnippet >= 5L)


###################################################
### code chunk number 9: seq-concatenate
###################################################
c(xSnippet, rev(xSnippet))
append(xSnippet, xSnippet, after=3)


###################################################
### code chunk number 10: aggregate
###################################################
xSnippet
aggregate(xSnippet, start=1:8, width=3, FUN=median)


###################################################
### code chunk number 11: shiftApply-cor
###################################################
cor(xRle, yRle)
shifts <- seq(235, 265, by=3)
corrs  <- shiftApply(shifts, yRle, xRle, FUN=cor)


###################################################
### code chunk number 12: figshiftcorrs
###################################################
plot(shifts, corrs)


###################################################
### code chunk number 13: Rle-vector-compare
###################################################
as.vector(object.size(xRle) / object.size(xVector))
identical(as.vector(xRle), xVector)


###################################################
### code chunk number 14: Rle-accessors
###################################################
head(runValue(xRle))
head(runLength(xRle))


###################################################
### code chunk number 15: Rle-ops
###################################################
xRle > 0
xRle + yRle
xRle > 0 | yRle > 0


###################################################
### code chunk number 16: Rle-summary
###################################################
range(xRle)
sum(xRle > 0 | yRle > 0)


###################################################
### code chunk number 17: Rle-math
###################################################
log1p(xRle)


###################################################
### code chunk number 18: Rle-cor
###################################################
cor(xRle, yRle)
shiftApply(249:251, yRle, xRle,
           FUN=function(x, y) {var(x, y) / (sd(x) * sd(y))})


###################################################
### code chunk number 19: DataFrame-extends-List
###################################################
showClass("DataFrame")


###################################################
### code chunk number 20: DataFrame
###################################################
df <- DataFrame(x=xRle, y=yRle)
sapply(df, class)
sapply(df, summary)
sapply(as.data.frame(df), summary)
endoapply(df, `+`, 0.5)


###################################################
### code chunk number 21: SessionInfo
###################################################
sessionInfo()


