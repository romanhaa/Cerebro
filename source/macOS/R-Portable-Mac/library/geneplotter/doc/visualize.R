### R code from vignette source 'visualize.Rnw'

###################################################
### code chunk number 1: getl
###################################################
library(geneplotter)


###################################################
### code chunk number 2: start
###################################################
data(sample.ExpressionSet)
eset = sample.ExpressionSet # legacy naming
mytt <- function(y, cov2) {
      ys <- split( y, cov2 )
      t.test( ys[[1]], ys[[2]] )
      }

ttout <- esApply(eset, 1, mytt, eset$type)
s1means <- sapply(ttout, function(x) x$estimate[1])
s2means <- sapply(ttout, function(x) x$estimate[2])
deciles <- quantile(c(s1means, s2means), probs=seq(0,1,.1))
s1class <- cut(s1means, deciles)
names(s1class) <- names(s1means)
s2class <- cut(s2means, deciles)
names(s2class) <- names(s2means)


###################################################
### code chunk number 3: f11
###################################################
cols <- dChip.colors(10)
def.par <- par(no.readonly = TRUE)# save default, for resetting...
nf <- layout(matrix(1:3,nr=1), widths=c(5,5,2))
chrObj <- buildChromLocation("hgu95av2")
cPlot(chrObj)
cColor(featureNames(eset), cols[s1class], chrObj)
cPlot(chrObj)
cColor(featureNames(eset), cols[s2class], chrObj)
image(1,1:10,matrix(1:10,nc=10),col=cols, axes=FALSE,
         xlab="", ylab="")
axis(2, at=(1:10), labels=levels(s1class), las=1)
par(def.par)


###################################################
### code chunk number 4: f22
###################################################
 par(mfrow=c(1,1)) 
 mycols <- c("red", "darkgreen", "blue")[eset$cov3]
 alongChrom(eset, "1", chrObj, plotFormat="cumulative", 
        col=mycols)


