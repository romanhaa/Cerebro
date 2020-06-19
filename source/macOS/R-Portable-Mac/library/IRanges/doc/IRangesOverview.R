### R code from vignette source 'IRangesOverview.Rnw'
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
## BiocManager::install("IRanges")


###################################################
### code chunk number 4: initialize
###################################################
library(IRanges)


###################################################
### code chunk number 5: iranges-constructor
###################################################
ir1 <- IRanges(start=1:10, width=10:1)
ir1
ir2 <- IRanges(start=1:10, end=11)
ir3 <- IRanges(end=11, width=10:1)
identical(ir1, ir2) && identical(ir1, ir3)
ir <- IRanges(c(1, 8, 14, 15, 19, 34, 40),
              width=c(12, 6, 6, 15, 6, 2, 7))
ir


###################################################
### code chunk number 6: iranges-start
###################################################
start(ir)


###################################################
### code chunk number 7: iranges-end
###################################################
end(ir)


###################################################
### code chunk number 8: iranges-width
###################################################
width(ir)


###################################################
### code chunk number 9: iranges-subset-numeric
###################################################
ir[1:4]


###################################################
### code chunk number 10: iranges-subset-logical
###################################################
ir[start(ir) <= 15]


###################################################
### code chunk number 11: plotRanges
###################################################
plotRanges <- function(x, xlim=x, main=deparse(substitute(x)),
                       col="black", sep=0.5, ...) 
{
  height <- 1
  if (is(xlim, "IntegerRanges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
  title(main)
  axis(1)
}


###################################################
### code chunk number 12: ir-plotRanges
###################################################
plotRanges(ir)


###################################################
### code chunk number 13: ranges-reduce
###################################################
reduce(ir)
plotRanges(reduce(ir))


###################################################
### code chunk number 14: rangeslist-contructor
###################################################
rl <- IRangesList(ir, rev(ir))


###################################################
### code chunk number 15: rangeslist-start
###################################################
start(rl)


###################################################
### code chunk number 16: bracket-ranges
###################################################
set.seed(0)
lambda <- c(rep(0.001, 4500), seq(0.001, 10, length=500),
            seq(10, 0.001, length=500))
xVector <- rpois(1e7, lambda)
yVector <- rpois(1e7, lambda[c(251:length(lambda), 1:250)])
xRle <- Rle(xVector)
yRle <- Rle(yVector)
irextract <- IRanges(start=c(4501, 4901) , width=100)
xRle[irextract]


###################################################
### code chunk number 17: overlap-ranges
###################################################
ol <- findOverlaps(ir, reduce(ir))
as.matrix(ol)


###################################################
### code chunk number 18: ranges-coverage
###################################################
cov <- coverage(ir)
plotRanges(ir)
cov <- as.vector(cov)
mat <- cbind(seq_along(cov)-0.5, cov)
d <- diff(cov) != 0
mat <- rbind(cbind(mat[d,1]+1, mat[d,2]), mat)
mat <- mat[order(mat[,1]),]
lines(mat, col="red", lwd=4)
axis(2)


###################################################
### code chunk number 19: ranges-shift
###################################################
shift(ir, 10)


###################################################
### code chunk number 20: ranges-narrow
###################################################
narrow(ir, start=1:5, width=2)


###################################################
### code chunk number 21: ranges-restrict
###################################################
restrict(ir, start=2, end=3)


###################################################
### code chunk number 22: ranges-threebands
###################################################
threebands(ir, start=1:5, width=2)


###################################################
### code chunk number 23: ranges-plus
###################################################
ir + seq_len(length(ir))


###################################################
### code chunk number 24: ranges-asterisk
###################################################
ir * -2 # double the width


###################################################
### code chunk number 25: ranges-disjoin
###################################################
disjoin(ir)
plotRanges(disjoin(ir))


###################################################
### code chunk number 26: ranges-disjointBins
###################################################
disjointBins(ir)


###################################################
### code chunk number 27: ranges-reflect
###################################################
reflect(ir, IRanges(start(ir), width=width(ir)*2))


###################################################
### code chunk number 28: ranges-flank
###################################################
flank(ir, width=seq_len(length(ir)))


###################################################
### code chunk number 29: ranges-gaps
###################################################
gaps(ir, start=1, end=50)
plotRanges(gaps(ir, start=1, end=50), c(1,50))


###################################################
### code chunk number 30: ranges-pgap
###################################################



###################################################
### code chunk number 31: ranges-union
###################################################



###################################################
### code chunk number 32: ranges-punion
###################################################



###################################################
### code chunk number 33: ranges-intersect
###################################################



###################################################
### code chunk number 34: ranges-pintersect
###################################################



###################################################
### code chunk number 35: ranges-setdiff
###################################################



###################################################
### code chunk number 36: ranges-psetdiff
###################################################



###################################################
### code chunk number 37: Views-constructors
###################################################
xViews <- Views(xRle, xRle >= 1)
xViews <- slice(xRle, 1)
xRleList <- RleList(xRle, 2L * rev(xRle))
xViewsList <- slice(xRleList, 1)


###################################################
### code chunk number 38: views-looping
###################################################
head(viewSums(xViews))
viewSums(xViewsList)
head(viewMaxs(xViews))
viewMaxs(xViewsList)


###################################################
### code chunk number 39: AtomicList-intro
###################################################
showClass("RleList")


###################################################
### code chunk number 40: list-construct
###################################################
args(IntegerList)
cIntList1 <- IntegerList(x=xVector, y=yVector)
cIntList1
sIntList2 <- IntegerList(x=xVector, y=yVector, compress=FALSE)
sIntList2
## sparse integer list
xExploded <- lapply(xVector[1:5000], function(x) seq_len(x))
cIntList2 <- IntegerList(xExploded)
sIntList2 <- IntegerList(xExploded, compress=FALSE)
object.size(cIntList2)
object.size(sIntList2)


###################################################
### code chunk number 41: list-length
###################################################
length(cIntList2)
Rle(lengths(cIntList2))


###################################################
### code chunk number 42: list-lapply
###################################################
system.time(sapply(xExploded, mean))
system.time(sapply(sIntList2, mean))
system.time(sapply(cIntList2, mean))
identical(sapply(xExploded, mean), sapply(sIntList2, mean))
identical(sapply(xExploded, mean), sapply(cIntList2, mean))


###################################################
### code chunk number 43: list-groupgenerics
###################################################
xRleList > 0
yRleList <- RleList(yRle, 2L * rev(yRle))
xRleList + yRleList
sum(xRleList > 0 | yRleList > 0)


###################################################
### code chunk number 44: list-endoapply
###################################################
safe.max <- function(x) { if(length(x)) max(x) else integer(0) }
endoapply(sIntList2, safe.max)
endoapply(cIntList2, safe.max)
endoapply(sIntList2, safe.max)[[1]]


###################################################
### code chunk number 45: SessionInfo
###################################################
sessionInfo()


