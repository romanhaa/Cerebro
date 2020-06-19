## ----foo,cache=FALSE,include=FALSE,echo=FALSE----
library(qvalue)
options(keep.source = TRUE, width = 48)
foo <- packageDescription("qvalue")

## ----citingqvalue, eval=FALSE-----------------
#  citation("qvalue")

## ----help_qvalue------------------------------
help(package="qvalue")

## ----quick_p----------------------------------
library(qvalue)
data(hedenfalk)
pvalues <- hedenfalk$p
qobj <- qvalue(p = pvalues)

## ----quick_stat-------------------------------
library(qvalue)
data(hedenfalk)
obs_stats <- hedenfalk$stat
null_stats <- hedenfalk$stat0
pvalues <- empPvals(stat = obs_stats, stat0 = null_stats)
qobj <- qvalue(p = pvalues)

## ----quick_access-----------------------------
qvalues <- qobj$qvalues
pi0 <- qobj$pi0
lfdr <- qobj$lfdr

## ----quick_sumviz, eval=FALSE-----------------
#  summary(qobj)
#  hist(qobj)
#  plot(qobj)

## ----load_qvalue------------------------------
data(hedenfalk)
names(hedenfalk)

## ----obsnullstat, dependson="load_qvalue"-----
null_stats <- hedenfalk$stat0
obs_stats <- hedenfalk$stat
pvalues <- empPvals(stat = obs_stats, stat0 = null_stats, pool = FALSE)

## ----pvalue_hist2, dependson=c("load_qvalue", "quick_p"),  fig.height=3, fig.width=5----
hist(hedenfalk$p, nclass = 20)

## ----pvalue_histBad, dependson=c("load_qvalue", "quick_p"), echo=FALSE, fig.height=3, fig.width=5----
set.seed(478)
p2 = c(hedenfalk$p, (runif(450, min=0.7, max=1))^(0.33))
somethingsWrong = list(p=p2)
hist(somethingsWrong$p, nclass=20, main="Problematic p-values", xlab="intentionally bad, simulated p-values")

## ----run_qvalue, dependson="load_qvalue"------
qobj <- qvalue(p = hedenfalk$p)

## ----outNames, dependson="run_qvalue"---------
names(qobj)

## ----summary_qvalue, dependson="run_qvalue"----
summary(qobj)

## ----pi0, dependson="run_qvalue"--------------
pi0 <- qobj$pi0

## ----pi0est, dependson="load_qvalue"----------
pi0 <- pi0est(p = hedenfalk$p, lambda = seq(0.1, 0.9, 0.1), pi0.method = "smoother" )
names(pi0)

## ----qvalue_ext, dependson="run_qvalue", fig.height=4----
qvalues <- qobj$qvalues

## ----tmp, dependson=c("qvalue_ext", "run_qvalue")----
max(qvalues[qobj$pvalues <= 0.01])

## ----fdrlevel, dependson="run_qvalue", eval=FALSE----
#  qobj_fdrlevel <- qvalue(p = hedenfalk$p, fdr.level = 0.1)
#  qobj$significant

## ----lfdr, dependson="run_qvalue"-------------
localFDR <- qobj$lfdr

## ----lfdr_f, dependson="load_qvalue"----------
localFDR <- lfdr(p = hedenfalk$p)

## ----plot_qobj, dependson=c("load_qvalue", "run_qvalue"), fig.width='\textwidth'----
plot(qobj)

## ----hist_qobj, dependson=c("load_qvalue", "run_qvalue"), fig.width='\textwidth'----
hist(qobj)

