### R code from vignette source 'prettyOutput.Rnw'

###################################################
### code chunk number 1: prettyOutput.Rnw:78-79
###################################################
options(width=70)


###################################################
### code chunk number 2: prettyOutput.Rnw:82-85
###################################################
library("annotate")
data(sample.ExpressionSet)
igenes <- featureNames(sample.ExpressionSet)[246:260]


###################################################
### code chunk number 3: prettyOutput.Rnw:89-98
###################################################
ug <- c("Hs.169284 // ---", "Hs.268515 // full length", "Hs.103419 // full length", "Hs.380429 // ---" ,"--- // ---",
        "Hs.169331 // full length", "Hs.381231 // full length", "Hs.283781 // full length", "--- // ---", "--- // ---",
        "Hs.3195 // full length", "--- // ---", "Hs.176660 // full length", "Hs.272484 // full length", "Hs.372679 // full length")
ll <- c("221823", "4330", "9637", "---", "---", "6331", "841", "27335", "---", "---", "6375", "---", "2543", "2578", "2215")
gb <- c("M57423", "Z70218", "L17328", "S81916", "U63332", "M77235", "X98175", "AB019392", "J03071", "D25272", "D63789",
        "D63789", "U19142", "U19147", "X16863")
sp <- c("P21108", "Q10571", "Q9UHY8", "Q16444", "---", "Q14524 /// Q8IZC9 /// Q8WTQ6 /// Q8WWN5 /// Q96J69", "Q14790", "Q9UBQ5",
        "---", "---", "P47992", "---", "Q13065 /// Q8IYC5", "Q13070", "O75015")



###################################################
### code chunk number 4: prettyOutput.Rnw:117-119
###################################################
gb
ll


###################################################
### code chunk number 5: prettyOutput.Rnw:127-130
###################################################
ug
ug <- sub(" //.*$", "", ug)
ug


###################################################
### code chunk number 6: prettyOutput.Rnw:141-144
###################################################
sp
sp <- strsplit(sub("---","&nbsp;",as.character(sp)), "///")
sp


###################################################
### code chunk number 7: expDat
###################################################
dat <- exprs(sample.ExpressionSet)[igenes,1:10]
FC <- rowMeans(dat[igenes,1:5]) - rowMeans(dat[igenes,6:10])
pval <- esApply(sample.ExpressionSet[igenes,1:10], 1, function(x) t.test(x[1:5], x[6:10])$p.value)
tstat <- esApply(sample.ExpressionSet[igenes,1:10], 1, function(x) t.test(x[1:5], x[6:10])$statistic)


###################################################
### code chunk number 8: prettyOutput.Rnw:170-177
###################################################
name <- c("hypothetical protein LOC221823",
          "meningioma (disrupted in balanced translocation) 1",
          "fasciculation and elongation protein zeta 2 (zygin II)",
          "Phosphoglycerate kinase {alternatively spliced}",
          "---","sodium channel, voltage-gated, type V, alpha polypeptide",
          "caspase 8, apoptosis-related cysteine protease","muscle specific gene","---","---","chemokine (C motif) ligand 1",
          "---","G antigen 1","G antigen 6","Fc fragment of IgG, low affinity IIIb, receptor for (CD16)")


###################################################
### code chunk number 9: prettyOutput.Rnw:179-182
###################################################
name
name <- gsub("---", "&nbsp;", name)
name


###################################################
### code chunk number 10: buildTable
###################################################
genelist <- list(igenes, ug, ll, gb, sp)
filename <- "Interesting_genes.html"
title <- "An Artificial Set of Interesting Genes"
othernames <- list(name, round(tstat, 2), round(pval, 3), round(FC, 1), round(dat, 2))
head <- c("Probe ID", "UniGene", "LocusLink", "GenBank", "SwissProt", "Gene Name", "t-statistic", "p-value",
          "Fold Change", "Sample 1", "Sample 2", "Sample 3", "Sample 4", "Sample 5", "Sample 6",
          "Sample 7", "Sample 8", "Sample 9", "Sample 10")
repository <- list("affy", "ug", "en", "gb", "sp")
htmlpage(genelist, filename, title, othernames, head, repository = repository)


###################################################
### code chunk number 11: prettyOutput.Rnw:208-209
###################################################
sessionInfo()


