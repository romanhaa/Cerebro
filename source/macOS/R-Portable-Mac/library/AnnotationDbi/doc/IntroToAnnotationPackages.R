## ----style, eval=TRUE, echo=FALSE, results='asis'--------------------------
BiocStyle::latex()

## ----include=FALSE---------------------------------------------------------
library(knitr)
opts_chunk$set(tidy=FALSE)

## ----loadChip--------------------------------------------------------------
suppressPackageStartupMessages({
    library(hgu95av2.db)
})

## ----listContents----------------------------------------------------------
ls("package:hgu95av2.db")

## ----show------------------------------------------------------------------
hgu95av2.db

## ----columns---------------------------------------------------------------
columns(hgu95av2.db)

## ----help, eval=FALSE------------------------------------------------------
#  help("SYMBOL")

## ----keytypes--------------------------------------------------------------
keytypes(hgu95av2.db)

## ----keys------------------------------------------------------------------
head(keys(hgu95av2.db, keytype="SYMBOL"))

## ----selectChip------------------------------------------------------------
#1st get some example keys
k <- head(keys(hgu95av2.db,keytype="PROBEID"))
# then call select
select(hgu95av2.db, keys=k, columns=c("SYMBOL","GENENAME"), keytype="PROBEID")

## ----mapIdsChip------------------------------------------------------------
#1st get some example keys
k <- head(keys(hgu95av2.db,keytype="PROBEID"))
# then call mapIds
mapIds(hgu95av2.db, keys=k, column=c("GENENAME"), keytype="PROBEID")

## ----selectOrg1------------------------------------------------------------
library(org.Hs.eg.db)
columns(org.Hs.eg.db)

## ----selectOrg2, eval=FALSE------------------------------------------------
#  help("SYMBOL") ## for explanation of these columns and keytypes values

## ----selectOrg3------------------------------------------------------------
keytypes(org.Hs.eg.db)
uniKeys <- head(keys(org.Hs.eg.db, keytype="UNIPROT"))
cols <- c("SYMBOL", "PATH")
select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="UNIPROT")

## ----selectData------------------------------------------------------------
load(system.file("extdata", "resultTable.Rda", package="AnnotationDbi"))
head(resultTable)

## ----selectOrgData---------------------------------------------------------
annots <- select(org.Hs.eg.db, keys=rownames(resultTable),
                 columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
resultTable <- merge(resultTable, annots, by.x=0, by.y="ENTREZID")
head(resultTable)

## ----selectGO--------------------------------------------------------------
library(GO.db)
GOIDs <- c("GO:0042254","GO:0044183")
select(GO.db, keys=GOIDs, columns="DEFINITION", keytype="GOID")

## ----selectTxDb------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
columns(txdb)
keytypes(txdb)
keys <- head(keys(txdb, keytype="GENEID"))
cols <- c("TXID", "TXSTART")
select(txdb, keys=keys, columns=cols, keytype="GENEID")


## ----selectEnsDb-----------------------------------------------------------
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
edb

## List all columns
columns(edb)

## List all keytypes
keytypes(edb)

## Get the first
keys <- head(keys(edb, keytype="GENEID"))

## Get the data
select(edb, keys=keys, columns=c("TXID", "TXSEQSTART", "TXBIOTYPE"), 
       keytype="GENEID")

## ----selectEnsDb.Y---------------------------------------------------------
## Retrieve all gene IDs of all lincRNAs encoded on chromosome Y
linkY <- keys(edb,
              filter=list(GeneBiotypeFilter("lincRNA"), SeqNameFilter("Y")))
length(linkY)

## We get now all transcripts for these genes.
txs <- select(edb, keys=linkY, columns=c("TXID", "TXSEQSTART", "TXBIOTYPE"),
              keytype="GENEID")
nrow(txs)

## Alternatively, we could specify/pass the filters with the keys argument.
txs <- select(edb, keys=list(GeneBiotypeFilter("lincRNA"), SeqNameFilter("Y")),
              columns=c("TXID", "TXSEQSTART", "TXBIOTYPE"))
nrow(txs)

## ----SessionInfo, echo=FALSE-----------------------------------------------
sessionInfo()

