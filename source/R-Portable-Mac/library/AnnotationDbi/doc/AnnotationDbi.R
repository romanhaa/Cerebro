## ----style, eval=TRUE, echo=FALSE, results='asis'--------------------------
BiocStyle::latex()

## ----include=FALSE---------------------------------------------------------
library(knitr)
opts_chunk$set(tidy=FALSE)

## ----<available schemas, results='hide'------------------------------------
library(DBI)
library(org.Hs.eg.db)
library(AnnotationForge)
available.dbschemas()

## ----setup0, results='hide', echo=FALSE-------------------------------
options(continue=" ", prompt="R> ", width=72L)

## ----setup, results='hide'--------------------------------------------
library(hgu95av2.db)

## ----objects----------------------------------------------------------
ls("package:hgu95av2.db")

## ----Question #1, echo=FALSE, results='hide'--------------------------
library(hgu95av2.db)
search()
hgu95av2_dbschema()
org.Hs.eg_dbschema()

## ----QAlisting--------------------------------------------------------
qcdata = capture.output(hgu95av2())
head(qcdata, 20)

## ----mapcounts, eval=FALSE--------------------------------------------
#  hgu95av2MAPCOUNTS

## ----envApiDemo1------------------------------------------------------
all_probes <- ls(hgu95av2ENTREZID)
length(all_probes)

set.seed(0xa1beef)
probes <- sample(all_probes, 5)
probes

## ----envApiDemo2------------------------------------------------------
hgu95av2ENTREZID[[probes[1]]]
hgu95av2ENTREZID$"31882_at"

syms <- unlist(mget(probes, hgu95av2SYMBOL))
syms

## ----helpDemo, eval= FALSE, results='hide'----------------------------
#  ?hgu95av2CHRLOC

## ----Question #2, echo=FALSE, results='hide'--------------------------
mget(probes, hgu95av2CHRLOC, ifnotfound=NA)[1:2]

## ----as.list, eval=FALSE----------------------------------------------
#  system.time(as.list(hgu95av2SYMBOL)[1:10])
#  
#  ## vs:
#  
#  system.time(as.list(hgu95av2SYMBOL[1:10]))

## ----show.revmap------------------------------------------------------
unlist(mget(syms, revmap(hgu95av2SYMBOL)))

## ----thisworks--------------------------------------------------------
as.list(revmap(hgu95av2PATH)["00300"])

## ----revmap2----------------------------------------------------------
x <- hgu95av2PATH
## except for the name, this is exactly revmap(x)
revx <- hgu95av2PATH2PROBE
revx2 <- revmap(x, objName="PATH2PROBE")
revx2
identical(revx, revx2)

as.list(revx["00300"])

## ----revmap2b---------------------------------------------------------
Term("GO:0000018")
Definition("GO:0000018")

## ----Question #3, echo=FALSE, results='hide'--------------------------
rs = ls(revmap(org.Hs.egREFSEQ))[4:6]
EGs = mget(rs, revmap(org.Hs.egREFSEQ), ifnotfound=NA)
##Then get the GO terms.
GOs = mget(unlist(EGs), org.Hs.egGO, ifnotfound=NA)
GOs
##Extract the GOIDs from this list:
GOIDs = as.character(unique(sapply(GOs, names)))
##Then look up what these terms are:
Term(GOIDs)

## ----toTable----------------------------------------------------------
head(toTable(hgu95av2GO[probes]))

## ----undirectedMethod-------------------------------------------------
toTable(x)[1:6, ]
toTable(revx)[1:6, ]

## ----directedMethods--------------------------------------------------
length(x)
length(revx)
allProbeSetIds <- keys(x)
allKEGGIds <- keys(revx)

## ----moreUndirectedMethods--------------------------------------------
junk <- Lkeys(x)        # same for all maps in hgu95av2.db (except pseudo-map
                        # MAPCOUNTS)
Llength(x)              # nb of Lkeys
junk <- Rkeys(x)        # KEGG ids for PATH/PATH2PROBE maps, GO ids for
                        # GO/GO2PROBE/GO2ALLPROBES maps, etc...
Rlength(x)              # nb of Rkeys

## ----moreKeysMethods--------------------------------------------------
x = hgu95av2ENTREZID[1:10]
## Directed methods
mappedkeys(x)           # mapped keys
count.mappedkeys(x)     # nb of mapped keys
## Undirected methods
mappedLkeys(x)          # mapped left keys
count.mappedLkeys(x)    # nb of mapped Lkeys

## ----isNA-------------------------------------------------------------
y = hgu95av2ENTREZID[isNA(hgu95av2ENTREZID)]     # usage like is.na()
Lkeys(y)[1:4]

## ----Question #4, echo=FALSE, results='hide'--------------------------

count.mappedLkeys(hgu95av2GO)
Llength(hgu95av2GO) - count.mappedLkeys(hgu95av2GO)
mappedLkeys(hgu95av2GO)[1]
toTable(hgu95av2GO["1000_at"])

## ----revmapUseCases---------------------------------------------------
x <- hgu95av2CHR
Rkeys(x)
chroms <- Rkeys(x)[23:24]
chroms
Rkeys(x) <- chroms
toTable(x)

## ----easy-------------------------------------------------------------
z <- as.list(revmap(x)[chroms])
names(z)
z[["Y"]]

## ----evilUnlist-------------------------------------------------------
chrs = c("12","6")
mget(chrs, revmap(hgu95av2CHR[1:30]), ifnotfound=NA)

## ----evilUnlist2------------------------------------------------------
unlist(mget(chrs, revmap(hgu95av2CHR[1:30]), ifnotfound=NA))

## ----evilUnlist3------------------------------------------------------
unlist2(mget(chrs, revmap(hgu95av2CHR[1:30]), ifnotfound=NA))

## ----cytogenetic2-----------------------------------------------------
x <- hgu95av2MAP
pbids <- c("38912_at", "41654_at", "907_at", "2053_at", "2054_g_at",
           "40781_at")
x <- subset(x, Lkeys=pbids, Rkeys="18q11.2")
toTable(x)

## ----coerce-----------------------------------------------------------
  pb2cyto <- as.character(x)
  pb2cyto[pbids]

## ----coercWarnings----------------------------------------------------
  cyto2pb <- as.character(revmap(x))

## ----multiProbes------------------------------------------------------
  ## How many probes?
  dim(hgu95av2ENTREZID)
  ## Make a mapping with multiple probes exposed
  multi <- toggleProbes(hgu95av2ENTREZID, "all")
  ## How many probes?
  dim(multi)

## ----multiProbes2-----------------------------------------------------
  ## Make a mapping with ONLY multiple probes exposed
  multiOnly <- toggleProbes(multi, "multiple")
  ## How many probes?
  dim(multiOnly)

  ## Then make a mapping with ONLY single mapping probes
  singleOnly <- toggleProbes(multiOnly, "single")
  ## How many probes?
  dim(singleOnly)

## ----multiProbes3-----------------------------------------------------
  ## Test the multiOnly mapping
  hasMultiProbes(multiOnly)
  hasSingleProbes(multiOnly)

  ## Test the singleOnly mapping
  hasMultiProbes(singleOnly)
  hasSingleProbes(singleOnly)

## ----orgSchema, results='hide'----------------------------------------
org.Hs.eg_dbschema()

## ----connObj, results='hide'------------------------------------------
org.Hs.eg_dbconn()

## ----connObj2, results='hide'-----------------------------------------
query <- "SELECT gene_id FROM genes LIMIT 10;"
result = dbGetQuery(org.Hs.eg_dbconn(), query)
result

## ----Question #5, echo=FALSE, results='hide'--------------------------
sql <- "SELECT gene_id, chromosome FROM genes AS g, chromosomes AS c WHERE g._id=c._id;"
dbGetQuery(org.Hs.eg_dbconn(),sql)[1:10,]

##OR
toTable(org.Hs.egCHR)[1:10,]

## ----complexEnv, eval=FALSE-------------------------------------------
#  ## Obtain SYMBOLS with at least one GO BP
#  ## annotation with evidence IMP, IGI, IPI, or IDA.
#  system.time({
#  bpids <- eapply(hgu95av2GO, function(x) {
#      if (length(x) == 1 && is.na(x))
#        NA
#      else {
#          sapply(x, function(z) {
#              if (z$Ontology == "BP")
#                z$GOID
#              else
#                NA
#              })
#      }
#  })
#  bpids <- unlist(bpids)
#  bpids <- unique(bpids[!is.na(bpids)])
#  g2p <- mget(bpids, hgu95av2GO2PROBE)
#  wantedp <- lapply(g2p, function(x) {
#      x[names(x) %in% c("IMP", "IGI", "IPI", "IDA")]
#  })
#  wantedp <- wantedp[sapply(wantedp, length) > 0]
#  wantedp <- unique(unlist(wantedp))
#  ans <- unlist(mget(wantedp, hgu95av2SYMBOL))
#  })
#  length(ans)
#  ans[1:10]

## ----schema, results='hide'-------------------------------------------
hgu95av2_dbschema()

## ----schema2, results='hide'------------------------------------------
hgu95av2ORGPKG

## ----schema3, results='hide'------------------------------------------
org.Hs.eg_dbschema()

## ----hgu95av2_org_join, tidy=FALSE------------------------------------
orgDBLoc = system.file("extdata", "org.Hs.eg.sqlite", package="org.Hs.eg.db")
attachSQL = paste("ATTACH '", orgDBLoc, "' AS orgDB;", sep = "")
dbGetQuery(hgu95av2_dbconn(), attachSQL)

## ----complexDb--------------------------------------------------------
system.time({
SQL <- "SELECT DISTINCT probe_id,symbol FROM probes, orgDB.gene_info AS gi, orgDB.genes AS g, orgDB.go_bp AS bp WHERE bp._id=g._id AND gi._id=g._id AND probes.gene_id=g.gene_id AND bp.evidence IN ('IPI', 'IDA', 'IMP', 'IGI')"
zz <- dbGetQuery(hgu95av2_dbconn(), SQL)
})
#its a good idea to always DETACH your database when you are finished...
dbGetQuery(hgu95av2_dbconn(), "DETACH orgDB"         )

## ----Question #6, echo=FALSE, results='hide'--------------------------
sql <- "SELECT gene_id, start_location, end_location, cytogenetic_location FROM genes AS g, chromosome_locations AS c, cytogenetic_locations AS cy WHERE g._id=c._id AND g._id=cy._id"
dbGetQuery(org.Hs.eg_dbconn(),sql)[1:10,]

## ----Question #7, echo=FALSE, results='hide'--------------------------
orgDBLoc = system.file("extdata", "org.Hs.eg.sqlite", package="org.Hs.eg.db")
attachSQL = paste("ATTACH '", orgDBLoc, "' AS orgDB;", sep = "")
dbGetQuery(hgu95av2_dbconn(), attachSQL)

goDBLoc = system.file("extdata", "GO.sqlite", package="GO.db")
attachSQL = paste("ATTACH '", goDBLoc, "' AS goDB;", sep = "")
dbGetQuery(hgu95av2_dbconn(), attachSQL)

SQL <- "SELECT DISTINCT p.probe_id, gi.symbol, gt.go_id, gt.definition
    FROM probes 
        AS p, orgDB.gene_info AS gi, orgDB.genes AS g, orgDB.go_bp 
        AS bp, goDB.go_term AS gt  
    WHERE bp._id=g._id AND gi._id=g._id AND p.gene_id=g.gene_id 
        AND bp.evidence IN ('IPI', 'IDA', 'IMP', 'IGI') AND gt.go_id=bp.go_id"
zz <- dbGetQuery(hgu95av2_dbconn(), SQL)

dbGetQuery(hgu95av2_dbconn(), "DETACH orgDB")
dbGetQuery(hgu95av2_dbconn(), "DETACH goDB")

## ----SessionInfo, echo=FALSE------------------------------------------
sessionInfo()

