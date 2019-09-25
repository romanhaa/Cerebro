### R code from vignette source 'GSVA.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: options
###################################################
options(width=60)
pdf.options(useDingbats=FALSE)


###################################################
### code chunk number 2: GSVA.Rnw:163-176
###################################################
library(GSVA)

p <- 20000    ## number of genes
n <- 30       ## number of samples
nGS <- 100    ## number of gene sets
min.sz <- 10  ## minimum gene set size
max.sz <- 100 ## maximum gene set size
X <- matrix(rnorm(p*n), nrow=p, dimnames=list(1:p, 1:n))
dim(X)
gs <- as.list(sample(min.sz:max.sz, size=nGS, replace=TRUE)) ## sample gene set sizes
gs <- lapply(gs, function(n, p) sample(1:p, size=n, replace=FALSE), p) ## sample gene sets
es.max <- gsva(X, gs, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
es.dif <- gsva(X, gs, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)


###################################################
### code chunk number 3: maxvsdif
###################################################
par(mfrow=c(1,2), mar=c(4, 4, 4, 1))
plot(density(as.vector(es.max)), main="Maximum deviation from zero",
     xlab="GSVA score", lwd=2, las=1, xaxt="n", xlim=c(-0.75, 0.75), cex.axis=0.8)
axis(1, at=seq(-0.75, 0.75, by=0.25), labels=seq(-0.75, 0.75, by=0.25), cex.axis=0.8)
plot(density(as.vector(es.dif)), main="Difference between largest\npositive and negative deviations",
     xlab="GSVA score", lwd=2, las=1, xaxt="n", xlim=c(-0.75, 0.75), cex.axis=0.8)
axis(1, at=seq(-0.75, 0.75, by=0.25), labels=seq(-0.75, 0.75, by=0.25), cex.axis=0.8)


###################################################
### code chunk number 4: GSVA.Rnw:307-312
###################################################
library(GSEABase)
library(GSVAdata)

data(c2BroadSets)
c2BroadSets


###################################################
### code chunk number 5: GSVA.Rnw:317-322
###################################################
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)


###################################################
### code chunk number 6: GSVA.Rnw:328-330
###################################################
cacheDir <- system.file("extdata", package="GSVA")
cachePrefix <- "cache4vignette_"


###################################################
### code chunk number 7: GSVA.Rnw:336-337 (eval = FALSE)
###################################################
## file.remove(paste(cacheDir, list.files(cacheDir, pattern=cachePrefix), sep="/"))


###################################################
### code chunk number 8: GSVA.Rnw:359-363
###################################################
data(leukemia)
leukemia_eset
head(pData(leukemia_eset))
table(leukemia_eset$subtype)


###################################################
### code chunk number 9: figIQR
###################################################
png(filename="GSVA-figIQR.png", width=500, height=500, res=150)
IQRs <- esApply(leukemia_eset, 1, IQR)
plot.ecdf(IQRs, pch=".", xlab="Interquartile range (IQR)", main="Leukemia data")
abline(v=quantile(IQRs, prob=0.5), lwd=2, col="red")
dev.off()


###################################################
### code chunk number 10: GSVA.Rnw:391-396
###################################################
filtered_eset <- nsFilter(leukemia_eset, require.entrez=TRUE, remove.dupEntrez=TRUE,
                          var.func=IQR, var.filter=TRUE, var.cutoff=0.5, filterByQuantile=TRUE,
                          feature.exclude="^AFFX")
filtered_eset
leukemia_filtered_eset <- filtered_eset$eset


###################################################
### code chunk number 11: GSVA.Rnw:410-413
###################################################
cache(leukemia_es <- gsva(leukemia_filtered_eset, c2BroadSets,
                           min.sz=10, max.sz=500, verbose=TRUE),
                           dir=cacheDir, prefix=cachePrefix)


###################################################
### code chunk number 12: GSVA.Rnw:423-425
###################################################
adjPvalueCutoff <- 0.001
logFCcutoff <- log2(2)


###################################################
### code chunk number 13: GSVA.Rnw:430-439
###################################################
design <- model.matrix(~ factor(leukemia_es$subtype))
colnames(design) <- c("ALL", "MLLvsALL")
fit <- lmFit(leukemia_es, design)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef="MLLvsALL", number=Inf)
DEgeneSets <- topTable(fit, coef="MLLvsALL", number=Inf,
                       p.value=adjPvalueCutoff, adjust="BH")
res <- decideTests(fit, p.value=adjPvalueCutoff)
summary(res)


###################################################
### code chunk number 14: GSVA.Rnw:445-455
###################################################
logFCcutoff <- log2(2)
design <- model.matrix(~ factor(leukemia_eset$subtype))
colnames(design) <- c("ALL", "MLLvsALL")
fit <- lmFit(leukemia_filtered_eset, design)
fit <- eBayes(fit)
allGenes <- topTable(fit, coef="MLLvsALL", number=Inf)
DEgenes <- topTable(fit, coef="MLLvsALL", number=Inf,
                    p.value=adjPvalueCutoff, adjust="BH", lfc=logFCcutoff)
res <- decideTests(fit, p.value=adjPvalueCutoff, lfc=logFCcutoff)
summary(res)


###################################################
### code chunk number 15: leukemiaVolcano
###################################################
png(filename="GSVA-leukemiaVolcano.png", width=800, height=500)
par(mfrow=c(1,2))
plot(allGeneSets$logFC, -log10(allGeneSets$P.Value), pch=".", cex=4, col=grey(0.75),
     main="Gene sets", xlab="GSVA enrichment score difference", ylab=expression(-log[10]~~~Raw~P-value))
abline(h=-log10(max(allGeneSets$P.Value[allGeneSets$adj.P.Val <= adjPvalueCutoff])),
       col=grey(0.5), lwd=1, lty=2)
points(allGeneSets$logFC[match(rownames(DEgeneSets), rownames(allGeneSets))],
       -log10(allGeneSets$P.Value[match(rownames(DEgeneSets), rownames(allGeneSets))]), pch=".",
       cex=4, col="red")
text(max(allGeneSets$logFC)*0.85,
         -log10(max(allGeneSets$P.Value[allGeneSets$adj.P.Val <= adjPvalueCutoff])),
         sprintf("%.1f%% FDR", 100*adjPvalueCutoff), pos=1)

plot(allGenes$logFC, -log10(allGenes$P.Value), pch=".", cex=4, col=grey(0.75),
     main="Genes", xlab="Log fold-change", ylab=expression(-log[10]~~~Raw~P-value))
abline(h=-log10(max(allGenes$P.Value[allGenes$adj.P.Val <= adjPvalueCutoff])),
       col=grey(0.5), lwd=1, lty=2)
abline(v=c(-logFCcutoff, logFCcutoff), col=grey(0.5), lwd=1, lty=2)
points(allGenes$logFC[match(rownames(DEgenes), rownames(allGenes))],
       -log10(allGenes$P.Value[match(rownames(DEgenes), rownames(allGenes))]), pch=".",
       cex=4, col="red")
text(max(allGenes$logFC)*0.85,
         -log10(max(allGenes$P.Value[allGenes$adj.P.Val <= adjPvalueCutoff])),
         sprintf("%.1f%% FDR", 100*adjPvalueCutoff), pos=1)
dev.off()


###################################################
### code chunk number 16: leukemiaHeatmapGeneSets
###################################################
png(filename="GSVA-leukemiaHeatmapGeneSets.png", width=500, height=500)
GSVAsco <- exprs(leukemia_es[rownames(DEgeneSets), ])
colorLegend <- c("darkred", "darkblue")
names(colorLegend) <- c("ALL", "MLL")
sample.color.map <- colorLegend[pData(leukemia_es)[, "subtype"]]
names(sample.color.map) <- colnames(GSVAsco)
sampleClustering <- hclust(as.dist(1-cor(GSVAsco, method="spearman")), method="complete")
geneSetClustering <- hclust(as.dist(1-cor(t(GSVAsco), method="pearson")), method="complete")
heatmap(GSVAsco, ColSideColors=sample.color.map, xlab="samples",
        ylab="Gene sets and pathways", margins=c(2, 20),
        labRow=substr(gsub("_", " ", gsub("^KEGG_|^REACTOME_|^BIOCARTA_", "", rownames(GSVAsco))), 1, 35),
        labCol="", scale="row",
        Colv=as.dendrogram(sampleClustering), Rowv=as.dendrogram(geneSetClustering))
legend("topleft", names(colorLegend), fill=colorLegend, inset=0.01, bg="white")
dev.off()


###################################################
### code chunk number 17: leukemiaHeatmapGenes
###################################################
png(filename="GSVA-leukemiaHeatmapGenes.png", width=500, height=500)
exps <- exprs(leukemia_eset[rownames(DEgenes), ])
colorLegend <- c("darkred", "darkblue")
names(colorLegend) <- c("ALL", "MLL")
sample.color.map <- colorLegend[pData(leukemia_eset)[, "subtype"]]
names(sample.color.map) <- colnames(exps)
sampleClustering <- hclust(as.dist(1-cor(exps, method="spearman")), method="complete")
geneClustering <- hclust(as.dist(1-cor(t(exps), method="pearson")), method="complete")
heatmap(exps, ColSideColors=sample.color.map, xlab="samples", ylab="Genes",
        labRow="", labCol="", scale="row", Colv=as.dendrogram(sampleClustering),
         Rowv=as.dendrogram(geneClustering), margins=c(2,2))
legend("topleft", names(colorLegend), fill=colorLegend, inset=0.01, bg="white")
dev.off()


###################################################
### code chunk number 18: GSVA.Rnw:562-569
###################################################
data(gbm_VerhaakEtAl)
gbm_eset
head(featureNames(gbm_eset))
table(gbm_eset$subtype)
data(brainTxDbSets)
sapply(brainTxDbSets, length)
lapply(brainTxDbSets, head)


###################################################
### code chunk number 19: GSVA.Rnw:574-575
###################################################
gbm_es <- gsva(gbm_eset, brainTxDbSets, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)


###################################################
### code chunk number 20: gbmSignature
###################################################
png(filename="GSVA-gbmSignature.png", width=700, height=500)
subtypeOrder <- c("Proneural", "Neural", "Classical", "Mesenchymal")
sampleOrderBySubtype <- sort(match(gbm_es$subtype, subtypeOrder), index.return=TRUE)$ix
subtypeXtable <- table(gbm_es$subtype)
subtypeColorLegend <- c(Proneural="red", Neural="green", Classical="blue", Mesenchymal="orange")
geneSetOrder <- c("astroglia_up", "astrocytic_up", "neuronal_up", "oligodendrocytic_up")
geneSetLabels <- gsub("_", " ", geneSetOrder)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
hmcol <- hmcol[length(hmcol):1]

heatmap(exprs(gbm_es)[geneSetOrder, sampleOrderBySubtype], Rowv=NA, Colv=NA,
        scale="row", margins=c(3,5), col=hmcol,
		    ColSideColors=rep(subtypeColorLegend[subtypeOrder], times=subtypeXtable[subtypeOrder]),
				labCol="", gbm_es$subtype[sampleOrderBySubtype],
        labRow=paste(toupper(substring(geneSetLabels, 1,1)), substring(geneSetLabels, 2), sep=""),
        cexRow=2, main=" \n ")
par(xpd=TRUE)
text(0.22,1.11, "Proneural", col="red", cex=1.2)
text(0.36,1.11, "Neural", col="green", cex=1.2)
text(0.48,1.11, "Classical", col="blue", cex=1.2)
text(0.66,1.11, "Mesenchymal", col="orange", cex=1.2)
mtext("Gene sets", side=4, line=0, cex=1.5)
mtext("Samples          ", side=1, line=4, cex=1.5)
dev.off()


###################################################
### code chunk number 21: GSVA.Rnw:649-684
###################################################
runSim <- function(p, n, gs.sz, S2N, fracDEgs) {
  sizeDEgs <- round(fracDEgs * gs.sz)
  group.n <- round(n / 2)

  sampleEffect <- rnorm(n, mean=0, sd=1)
  sampleEffectDE <- rnorm(n, mean=S2N, sd=0.5)
  probeEffect <- rnorm(p, mean=0, sd=1)
  noise <- matrix(rnorm(p*n, mean=0, sd=1), nrow=p, ncol=n)
  noiseDE <- matrix(rnorm(p*n, mean=0, sd=1), nrow=p, ncol=n)
  M <- outer(probeEffect, sampleEffect, "+") + noise
  M2 <- outer(probeEffect, sampleEffectDE, "+") + noiseDE
  M[1:sizeDEgs, 1:group.n] <- M2[1:sizeDEgs, 1:group.n]

  rownames(M) <- paste0("g", 1:nrow(M))
  geneSets <- list(H1GeneSet=paste0("g", 1:(gs.sz)),
                   H0GeneSet=paste0("g", (gs.sz+1):(2*gs.sz)))

  es.gsva <- gsva(M, geneSets, verbose=FALSE, parallel.sz=1)
  es.ss <- gsva(M, geneSets, method="ssgsea", verbose=FALSE, parallel.sz=1)
  es.z <- gsva(M, geneSets, method="zscore", verbose=FALSE, parallel.sz=1)
  es.plage <- gsva(M, geneSets, method="plage", verbose=FALSE, parallel.sz=1)

  h1.gsva.pval <- t.test(es.gsva["H1GeneSet", 1:group.n],es.gsva["H1GeneSet", (group.n+1):n])$p.value
  h1.ssgsea.pval <- t.test(es.ss["H1GeneSet", 1:group.n],es.ss["H1GeneSet", (group.n+1):n])$p.value
  h1.zscore.pval <- t.test(es.z["H1GeneSet", 1:group.n],es.z["H1GeneSet", (group.n+1):n])$p.value
  h1.plage.pval <- t.test(es.plage["H1GeneSet", 1:group.n],es.plage["H1GeneSet", (group.n+1):n])$p.value

  h0.gsva.pval <- t.test(es.gsva["H0GeneSet", 1:group.n],es.gsva["H0GeneSet", (group.n+1):n])$p.value
  h0.ssgsea.pval <- t.test(es.ss["H0GeneSet", 1:group.n],es.ss["H0GeneSet", (group.n+1):n])$p.value
  h0.zscore.pval <- t.test(es.z["H0GeneSet", 1:group.n],es.z["H0GeneSet", (group.n+1):n])$p.value
  h0.plage.pval <- t.test(es.plage["H0GeneSet", 1:group.n],es.plage["H0GeneSet", (group.n+1):n])$p.value

  c(h1.gsva.pval, h1.ssgsea.pval, h1.zscore.pval, h1.plage.pval,
    h0.gsva.pval, h0.ssgsea.pval, h0.zscore.pval, h0.plage.pval)
}


###################################################
### code chunk number 22: GSVA.Rnw:691-696
###################################################
estPwrTypIerr <- function(pvals, alpha=0.05) {
  N <- ncol(pvals)
  c(1 - sum(pvals[1, ] > alpha)/N, 1 - sum(pvals[2, ] > alpha)/N,1 - sum(pvals[3, ] > alpha)/N, 1 - sum(pvals[4, ] > alpha)/N,
        sum(pvals[5, ] <= alpha)/N, sum(pvals[6, ] <= alpha)/N, sum(pvals[7, ] <= alpha)/N, sum(pvals[8, ] <= alpha)/N)
}


###################################################
### code chunk number 23: GSVA.Rnw:705-726
###################################################
set.seed(1234)

exp1 <- cbind(estPwrTypIerr(replicate(60, runSim(1000, 10, gs.sz=30, S2N=0.5, fracDEgs=0.5))),
              estPwrTypIerr(replicate(60, runSim(1000, 20, gs.sz=30, S2N=0.5, fracDEgs=0.5))),
              estPwrTypIerr(replicate(60, runSim(1000, 40, gs.sz=30, S2N=0.5, fracDEgs=0.5))),
              estPwrTypIerr(replicate(60, runSim(1000, 60, gs.sz=30, S2N=0.5, fracDEgs=0.5))))

exp2 <- cbind(estPwrTypIerr(replicate(60, runSim(1000, 10, gs.sz=30, S2N=1.0, fracDEgs=0.5))),
              estPwrTypIerr(replicate(60, runSim(1000, 20, gs.sz=30, S2N=1.0, fracDEgs=0.5))),
              estPwrTypIerr(replicate(60, runSim(1000, 40, gs.sz=30, S2N=1.0, fracDEgs=0.5))),
              estPwrTypIerr(replicate(60, runSim(1000, 60, gs.sz=30, S2N=1.0, fracDEgs=0.5))))

exp3 <- cbind(estPwrTypIerr(replicate(60, runSim(1000, 10, gs.sz=30, S2N=0.5, fracDEgs=0.8))),
              estPwrTypIerr(replicate(60, runSim(1000, 20, gs.sz=30, S2N=0.5, fracDEgs=0.8))),
              estPwrTypIerr(replicate(60, runSim(1000, 40, gs.sz=30, S2N=0.5, fracDEgs=0.8))),
              estPwrTypIerr(replicate(60, runSim(1000, 60, gs.sz=30, S2N=0.5, fracDEgs=0.8))))

exp4 <- cbind(estPwrTypIerr(replicate(60, runSim(1000, 10, gs.sz=30, S2N=1.0, fracDEgs=0.8))),
              estPwrTypIerr(replicate(60, runSim(1000, 20, gs.sz=30, S2N=1.0, fracDEgs=0.8))),
              estPwrTypIerr(replicate(60, runSim(1000, 40, gs.sz=30, S2N=1.0, fracDEgs=0.8))),
              estPwrTypIerr(replicate(60, runSim(1000, 60, gs.sz=30, S2N=1.0, fracDEgs=0.8))))


###################################################
### code chunk number 24: powertype1errsim
###################################################
plotPower <- function(statpower, main, legendposition="bottomright", ...) {
  plot(statpower[1,], ylim=c(0, 1.0), type="b", lwd=2, pch=1, main=main,
       col="blue", ylab="Statistcal Power", xlab="Sample Size", xaxt="n")
  lines(statpower[2,], col="red", type="b", lwd=2, pch=2)
  lines(statpower[3,], col="darkgreen", type="b", lwd=2, pch=3)
  lines(statpower[4,], col="lightgreen", type="b", lwd=2, pch=4)
  if (!is.null(legendposition))
    legend(legendposition, c("GSVA","ssGSEA","z-score","PLAGE"),
           col=c("blue","red","darkgreen","lightgreen"),pch=1:4,lty=1,lwd=2,inset=0.02)
  axis(1,at=1:4, labels=c("10","20","40","60"))
}

plotType1Error <- function(tmp, title, legendposition="bottomright", alpha=0.05, ...){
  plot(tmp[5,],ylim=c(0, 0.2),type="b",lwd=2,pch=1,
       col="blue",ylab="Empirical Type-I Error",xlab="Sample Size",xaxt="n",main=title, ...)
  lines(tmp[6,],col="red",type="b",lwd=2,pch=2)
  lines(tmp[7,],col="darkgreen",type="b",lwd=2,pch=3)
  lines(tmp[8,],col="lightgreen",type="b",lwd=2,pch=4)
  if (!is.null(legendposition))
    legend(legendposition,c("GSVA","ssGSEA","z-score","PLAGE"),col=c("blue","red","darkgreen","lightgreen"),pch=1:4,lty=1,lwd=2,inset=0.02)
  axis(1,at=c(1:dim(tmp)[2]), labels=c("10","20","40","60"))
  abline(h=alpha, lty=2)
}

labelPlot <- function(lab, font, cex, offsetx=0.05, offsety=0.05) {
  par(xpd=TRUE)
  w <- par("usr")[2] - par("usr")[1]
  h <- par("usr")[4] - par("usr")[3]
  text(par("usr")[1]-w*offsetx, par("usr")[4]+h*offsety, lab, font=font, cex=cex)
  par(xpd=FALSE)
}

par(mfrow=c(4,2), mar=c(4, 4, 2, 1))
plotPower(exp1, main="", legendposition=NULL, las=1)
labelPlot("A", 2, 2, 0.2, 0.15)
plotType1Error(exp1,"",legendposition="topright", las=1)
labelPlot("B", 2, 2, 0.2, 0.15)
plotPower(exp2, main="", legendposition=NULL, las=1)
labelPlot("C", 2, 2, 0.2, 0.15)
plotType1Error(exp2,"",legendposition="topright", las=1)
labelPlot("D", 2, 2, 0.2, 0.15)
plotPower(exp3, main="", legendposition=NULL, las=1)
labelPlot("E", 2, 2, 0.2, 0.15)
plotType1Error(exp3,"",legendposition="topright", las=1)
labelPlot("F", 2, 2, 0.2, 0.15)
plotPower(exp4, main="", legendposition=NULL, las=1)
labelPlot("G", 2, 2, 0.2, 0.15)
plotType1Error(exp4,"",legendposition="topright", las=1)
labelPlot("H", 2, 2, 0.2, 0.15)


###################################################
### code chunk number 25: GSVA.Rnw:808-814
###################################################
data(commonPickrellHuang)

stopifnot(identical(featureNames(huangArrayRMAnoBatchCommon_eset),
                    featureNames(pickrellCountsArgonneCQNcommon_eset)))
stopifnot(identical(sampleNames(huangArrayRMAnoBatchCommon_eset),
                    sampleNames(pickrellCountsArgonneCQNcommon_eset)))


###################################################
### code chunk number 26: GSVA.Rnw:820-824
###################################################
canonicalC2BroadSets <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
                                      grep("^REACTOME", names(c2BroadSets)),
                                      grep("^BIOCARTA", names(c2BroadSets)))]
canonicalC2BroadSets


###################################################
### code chunk number 27: <
###################################################
data(genderGenesEntrez)

MSY <- GeneSet(msYgenesEntrez, geneIdType=EntrezIdentifier(),
               collectionType=BroadCollection(category="c2"), setName="MSY")
MSY
XiE <- GeneSet(XiEgenesEntrez, geneIdType=EntrezIdentifier(),
               collectionType=BroadCollection(category="c2"), setName="XiE")
XiE


canonicalC2BroadSets <- GeneSetCollection(c(canonicalC2BroadSets, MSY, XiE))
canonicalC2BroadSets


###################################################
### code chunk number 28: <
###################################################
esmicro <- gsva(huangArrayRMAnoBatchCommon_eset, canonicalC2BroadSets, min.sz=5, max.sz=500,
                mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
dim(esmicro)
esrnaseq <- gsva(pickrellCountsArgonneCQNcommon_eset, canonicalC2BroadSets, min.sz=5, max.sz=500,
                 kcdf="Poisson", mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
dim(esrnaseq)


###################################################
### code chunk number 29: GSVA.Rnw:864-879
###################################################
library(edgeR)

data(annotEntrez220212)
head(annotEntrez220212)

cpm <- cpm(exprs(pickrellCountsArgonneCQNcommon_eset))
dim(cpm)

common <- intersect(rownames(cpm), rownames(annotEntrez220212))
length(common)

rpkm <- sweep(cpm[common, ], 1, annotEntrez220212[common, "Length"] / 10^3, FUN="/")
dim(rpkm)

dim(huangArrayRMAnoBatchCommon_eset[rownames(rpkm), ])


###################################################
### code chunk number 30: GSVA.Rnw:885-894
###################################################
corsrowsgene <- sapply(1:nrow(huangArrayRMAnoBatchCommon_eset[rownames(rpkm), ]),
                       function(i, expmicro, exprnaseq) cor(expmicro[i, ], exprnaseq[i, ], method="pearson"),
                       exprs(huangArrayRMAnoBatchCommon_eset[rownames(rpkm), ]), log2(rpkm+0.1))
names(corsrowsgene) <- rownames(rpkm)

corsrowsgs <- sapply(1:nrow(esmicro),
                     function(i, esmicro, esrnaseq) cor(esmicro[i, ], esrnaseq[i, ], method="spearman"),
                     exprs(esmicro), exprs(esrnaseq))
names(corsrowsgs) <- rownames(esmicro)


###################################################
### code chunk number 31: RNAseqComp
###################################################
png(filename="GSVA-RNAseqComp.png", width=1100, height=1100, res=150)
par(mfrow=c(2,2), mar=c(4, 5, 3, 2))
hist(corsrowsgene, xlab="Spearman correlation", main="Gene level\n(RNA-seq RPKM vs Microarray RMA)",
     xlim=c(-1, 1), col="grey", las=1)
par(xpd=TRUE)
text(par("usr")[1]*1.5, par("usr")[4]*1.1, "A", font=2, cex=2)
par(xpd=FALSE)
hist(corsrowsgs, xlab="Spearman correlation", main="Gene set level\n(GSVA enrichment scores)",
     xlim=c(-1, 1), col="grey", las=1)
par(xpd=TRUE)
text(par("usr")[1]*1.5, par("usr")[4]*1.1, "B", font=2, cex=2)
par(xpd=FALSE)
plot(exprs(esrnaseq)["MSY", ], exprs(esmicro)["MSY", ], xlab="GSVA scores RNA-seq", ylab="GSVA scores microarray",
     main=sprintf("MSY R=%.2f", cor(exprs(esrnaseq)["MSY", ], exprs(esmicro)["MSY", ])), las=1, type="n")
     sprintf("MSY R=%.2f", cor(exprs(esrnaseq)["MSY", ], exprs(esmicro)["MSY", ]))
abline(lm(exprs(esmicro)["MSY", ] ~ exprs(esrnaseq)["MSY", ]), lwd=2, lty=2, col="grey")
points(exprs(esrnaseq)["MSY", pickrellCountsArgonneCQNcommon_eset$Gender == "Female"],
       exprs(esmicro)["MSY", huangArrayRMAnoBatchCommon_eset$Gender == "Female"], col="red", pch=21, bg="red", cex=1)
points(exprs(esrnaseq)["MSY", pickrellCountsArgonneCQNcommon_eset$Gender == "Male"],
       exprs(esmicro)["MSY", huangArrayRMAnoBatchCommon_eset$Gender == "Male"], col="blue", pch=21, bg="blue", cex=1)
par(xpd=TRUE)
text(par("usr")[1]*1.5, par("usr")[4]*1.1, "C", font=2, cex=2)
par(xpd=FALSE)
plot(exprs(esrnaseq)["XiE", ], exprs(esmicro)["XiE", ], xlab="GSVA scores RNA-seq", ylab="GSVA scores microarray",
     main=sprintf("XiE R=%.2f", cor(exprs(esrnaseq)["XiE", ], exprs(esmicro)["XiE", ])), las=1, type="n")
abline(lm(exprs(esmicro)["XiE", ] ~ exprs(esrnaseq)["XiE", ]), lwd=2, lty=2, col="grey")
points(exprs(esrnaseq["XiE", pickrellCountsArgonneCQNcommon_eset$Gender == "Female"]),
       exprs(esmicro)["XiE", huangArrayRMAnoBatchCommon_eset$Gender == "Female"], col="red", pch=21, bg="red", cex=1)
points(exprs(esrnaseq)["XiE", pickrellCountsArgonneCQNcommon_eset$Gender == "Male"],
       exprs(esmicro)["XiE", huangArrayRMAnoBatchCommon_eset$Gender == "Male"], col="blue", pch=21, bg="blue", cex=1)
par(xpd=TRUE)
text(par("usr")[1]*1.5, par("usr")[4]*1.1, "D", font=2, cex=2)
par(xpd=FALSE)
dev.off()


###################################################
### code chunk number 32: info
###################################################
toLatex(sessionInfo())


