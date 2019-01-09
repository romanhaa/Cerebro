### R code from vignette source 'diptest-issues.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(SweaveHooks= list(fig=function() par(mar=c(5.1, 4.1, 1.1, 2.1))),
        width = 75,
        digits = 7, # <-- here, keep R's default!
        prompt = "R> ", # <- "yuck!" - required by JSS
        continue=" ")
set.seed(47)
Sys.setenv(LANGUAGE = "en")
if(.Platform$OS.type != "windows")
  Sys.setlocale("LC_MESSAGES","C")

## In order to  save() and load() expensive results
thisDir <- system.file('doc', package='diptest')
xtraDir <- if(Sys.getenv("USER") == "maechler")
    "~/R/Pkgs/diptest/stuff" else thisDir
res1.file <- file.path(thisDir, "aggr_results.Rdata")



###################################################
### code chunk number 2: diagnose-lib
###################################################
if(nzchar(Sys.getenv("R_MM_PKG_CHECKING"))) print( .libPaths() )


###################################################
### code chunk number 3: dip_n-is-5
###################################################
getOption("SweaveHooks")[["fig"]]()
require("diptest") # after installing it ..
D5 <- replicate(10000, dip(runif(5)))
hist(D5, breaks=128, main = "Histogram of  replicate(10'000, dip(runif(5))))")


###################################################
### code chunk number 4: dip_n-is-8
###################################################
getOption("SweaveHooks")[["fig"]]()
D8 <- replicate(10000, dip(runif(8)))
hist(D8, breaks=128, main = "Histogram of  replicate(10'000, dip(runif(8))))")


###################################################
### code chunk number 5: sim--n-eq-11 (eval = FALSE)
###################################################
##  set.seed(11)
##  n <- 11
##  B.s11 <- 500000
##  D11 <- replicate(B.s11, dip(runif(n)))


###################################################
### code chunk number 6: 2nd-small-sample-phenomen--n-eq-11
###################################################
if(file.exists(ff <- file.path(thisDir, "hist-D11.rda"))) {
  load(ff)
} else { ## takes a few minutes
 set.seed(11)
 n <- 11
 B.s11 <- 500000
 D11 <- replicate(B.s11, dip(runif(n)))
  hD11 <- hist(D11, breaks=1e-6+(63:298)/(2*11*64), plot=FALSE)
  save(hD11, n, B.s11, file= ff)
}


###################################################
### code chunk number 7: 2nd-small-sample-phenomen--n-eq-11
###################################################
getOption("SweaveHooks")[["fig"]]()
B.str <- format(B.s11, sci=FALSE, big.mark="'")
plot(hD11, main = "",
     ## main = sprintf("Histogram of  replicate(%s, dip(runif(%d)))", B.str, n),
     border=NA, col="dark gray",
     xlab = substitute("Dip" ~~ D[.N.](U(group("[",list(0,1),"]"))), list(.N. = n)))
title(xlab= substitute(B == .B.SIM. ~ "replicates", list(.B.SIM. = B.str)),
      adj = .88)
lcol <- adjustcolor("orange4", 0.4)
abline(v = (1:3)/(2*n), col=lcol, lty=3, lwd=2)
axis(1, pos=0, at = (1:3)/(2*n),
     labels = expression(1/22, 2/22, 3/22), col=lcol, col.axis=lcol)


###################################################
### code chunk number 8: sqrt-n-qdip
###################################################
getOption("SweaveHooks")[["fig"]]()
data(qDiptab)
dnqd <- dimnames(qDiptab)
(nn. <- as.integer(dnqd[["n"]]))
matplot(nn., qDiptab*sqrt(nn.), type ="o", pch=1, cex = 0.4,
        log="x", xlab="n   [log scaled]",
        ylab = expression(sqrt(n) %*% q[D[n]]))
## Note that  1/2n  is the first possible value (with finite mass),,
## clearly visible for (very) small n:
lines(nn., sqrt(nn.)/(2*nn.), col=adjustcolor("yellow2",0.5), lwd=3)

P.p <- as.numeric(print(noquote(dnqd[["Pr"]])))
## Now look at one well known data set:
D <- dip(x <- faithful$waiting)
n <- length(x)
points(n, sqrt(n)*D, pch=13, cex=2, col= adjustcolor("blue2",.5), lwd=2)
## a simulated (approximate) $p$-value for D  is
mean(D <= replicate(10000, dip(runif(n)))) ## ~ 0.002


###################################################
### code chunk number 9: interpolate-dip-table
###################################################
## We are in this interval:
n0 <- nn.[i.n <- findInterval(n, nn.)]
n1 <- nn.[i.n +1] ; c(n0, n1)
f.n <- (n - n0)/(n1 - n0)# in [0, 1]
## Now "find" y-interval:
y.0 <- sqrt(n0)* qDiptab[i.n  ,]
y.1 <- sqrt(n1)* qDiptab[i.n+1,]
(Pval <- 1 - approx(y.0 + f.n*(y.1 - y.0),
                    P.p,
                    xout = sqrt(n) * D)[["y"]])
## 0.018095


###################################################
### code chunk number 10: statfac-dip.test
###################################################
data(statfaculty)
dip.test(statfaculty)


###################################################
### code chunk number 11: sessionInfo
###################################################
toLatex(sessionInfo())


