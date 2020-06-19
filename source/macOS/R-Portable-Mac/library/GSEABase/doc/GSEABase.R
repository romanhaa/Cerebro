## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"--------------------
BiocStyle::latex()

## ----preliminaries, echo=FALSE---------------------------------------------
## FIXME: <adjMat> adjacency matrix -- color w. +/- 1
## FIXME: limma topTable --> GeneColorSet
## w. verbose=TRUE
suppressPackageStartupMessages({
    library(GSEABase)
    library(hgu95av2.db)
    library(GO.db)
})

## ----GeneSet---------------------------------------------------------------
data(sample.ExpressionSet) # from Biobase
egs <- GeneSet(sample.ExpressionSet[201:250,], setName="Sample")
egs

## ----geneIds---------------------------------------------------------------
head(geneIds(egs))

## ----details---------------------------------------------------------------
details(egs)

## --------------------------------------------------------------------------
  ## FIXME: GeneSet(AnnotationIdentifier("hgu95av2")) --> non-empty
  ## FIXME: GeneSet(AnnotationIdentifier("hgu95av2"),
  ## collectionType=GOCollection()) filters on GOCollection (or KEGG)

## ----GeneSet-methods-------------------------------------------------------
showMethods("GeneSet", inherited=FALSE)

## ----GeneIdentifierTypes---------------------------------------------------
names(slot(getClass("GeneIdentifierType"), "subclasses"))

## ----mapIdentifiers--------------------------------------------------------
mapIdentifiers(egs, EntrezIdentifier())

## ----GeneSet_Identifiers---------------------------------------------------
library(annotate)                       # getEG
eids <- unique(getEG(geneIds(egs), "hgu95av2"))
eids <- eids[!is.na(eids)]
GeneSet(EntrezIdentifier(), geneIds=as.character(eids))

## ----CollectionType--------------------------------------------------------
names(slot(getClass("CollectionType"), "subclasses"))

## ----GOCollection----------------------------------------------------------
GeneSet(GOCollection(c("GO:0005488", "GO:0019825"),
                     evidenceCode="IDA"),
        geneIdType=EntrezIdentifier("org.Hs.eg.db"),
        setName="Sample GO Collection")

## ----Broad-----------------------------------------------------------------
fl <- system.file("extdata", "Broad1.xml", package="GSEABase")
bgs <- GeneSet(BroadCollection(), urls=fl)
bgs

## ----Broad-to-annotation---------------------------------------------------
bgs1 <- mapIdentifiers(bgs, AnnotationIdentifier("hgu95av2"))
bgs1

## ----subset----------------------------------------------------------------
bgs[1:5]
bgs[c("GALNS", "LOC646365")]

## ----egs-bgs---------------------------------------------------------------
egs & bgs1

## ----subset-ExpressionSet--------------------------------------------------
sample.ExpressionSet[bgs,]

## ----GeneColorSet-setup, echo=FALSE, results="hide"------------------------
conn <- textConnection("
Entrez ID, Gene Symbol, Expression level, Phenotype response
##used to be MRP2
1244, ABCC2, Increase, Resistant
538, ATP7A, Increase, Resistant
540, ATP7B, Increase, Resistant
9961, MVP, Increase, Resistant
##the LRP below must be MVP
##LRP, Increase, Resistant - need to know which one
7507,XPA, Increase, Resistant
2067, ERCC1, Increase, Resistant
##TOP, Increase, Resistant  - need to know which one, notes say II
672, BRCA1, Increase, Resistant
3725, JUN, Increase, Resistant
#GCS, Increase, Resistant  - my notes say alpha-GCS - so which one?
##I only found gamma at PubMed as being related
2730, GCLM, Increase, Resistant")
tbl <- read.csv(conn, strip.white=TRUE, comment.char="#")
close(conn) 
unlink(conn)

## ----GeneColorSet-phenotype------------------------------------------------
tbl

## ----GeneColorSet-constructor----------------------------------------------
gcs <- GeneColorSet(EntrezIdentifier(),
                    setName="A color set",
                    geneIds=as.character(tbl$Entrez.ID),
                    phenotype="Cisplatin resistance",
                    geneColor=tbl$Expression.level,
                    phenotypeColor=tbl$Phenotype.response)
gcs

## ----GeneSetCollection-----------------------------------------------------
gsc <- GeneSetCollection(sample.ExpressionSet[201:250,], setType=GOCollection())
gsc
gsc[["GO:0005737"]]

## ----GeneSetCollection-GOCollection----------------------------------------
GeneSetCollection(sample.ExpressionSet[201:300,],
                  setType=GOCollection(evidenceCode="IMP"))

## ----GeneSetCollection-BroadCollection-------------------------------------
  ## FIXME: BroadCollection default to paste("c", 1:4, sep="")
  ## FIXME: GeneSetCollection(BroadCollection(), urls=fl); filters on bcCategory
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
gss <- getBroadSets(fl)
gss
names(gss)

## ----mapIds-GeneSetCollection----------------------------------------------
gsc <- mapIdentifiers(gsc, EntrezIdentifier())
gsc
gsc[["GO:0005737"]]

## ----ReportingTools--------------------------------------------------------
## 'interesting' gene sets
idx <- sapply(gsc, function(x) length(geneIds(x))) > 2

library(ReportingTools)
gscReport <- HTMLReport(
    shortName="gsc_example",
    title="GSEABase Vignette GeneSetCollection",
    basePath=tempdir())
publish(gsc[idx], gscReport, annotation.db="org.Hs.eg")
url <- finish(gscReport)

## ----ReportingTools-view, eval=FALSE---------------------------------------
#  browseURL(url)

