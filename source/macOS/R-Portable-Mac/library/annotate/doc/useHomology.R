### R code from vignette source 'useHomology.Rnw'

###################################################
### code chunk number 1: loadAndListPackages
###################################################
library(DBI)
library("hom.Hs.inp.db")
ls("package:hom.Hs.inp.db")


###################################################
### code chunk number 2: listMapping
###################################################
as.list(hom.Hs.inpMUSMU[1:4])


###################################################
### code chunk number 3: mapGeneExample
###################################################
# load the organism annotation data for human
library(org.Hs.eg.db)

# get the entrex gene ID and ensembl protein id for gene symbol "MSX2"
select(org.Hs.eg.db, 
       keys="MSX2", 
       columns=c("ENTREZID","ENSEMBLPROT"), 
       keytype="SYMBOL")

# use the inparanoid package to get the mouse gene that is considered 
# equivalent to ensembl protein ID "ENSP00000239243"
select(hom.Hs.inp.db, 
       keys="ENSP00000239243", 
       columns="MUS_MUSCULUS", 
       keytype="HOMO_SAPIENS")

# load the organism annotation data for mouse
library(org.Mm.eg.db)

# get the entrez gene ID and gene Symbol for "ENSMUSP00000021922"
select(org.Mm.eg.db, 
       keys="ENSMUSP00000021922", 
       columns=c("ENTREZID","SYMBOL"), 
       keytype="ENSEMBLPROT")


###################################################
### code chunk number 4: seedPairExample
###################################################
mget("ENSP00000301011", hom.Hs.inpMUSMU)


###################################################
### code chunk number 5: seedPairExample2
###################################################
# make a connection to the human database
mycon <- hom.Hs.inp_dbconn()
# make a list of all the tables that are available in the DB
head(dbListTables(mycon))
# make a list of the columns in the table of interest
dbListFields(mycon, "mus_musculus")


###################################################
### code chunk number 6: seedPairExample3
###################################################
#make a query that will let us see which clust_id we need
sql <- "SELECT * FROM mus_musculus WHERE inp_id = 'ENSP00000301011';"
#retrieve the data
dataOut <- dbGetQuery(mycon, sql)
dataOut


###################################################
### code chunk number 7: seedPairExample4
###################################################
#make a query that will let us see all the data that is affiliated with a clust id
sql <- "SELECT * FROM mus_musculus WHERE clust_id = '1731';"
#retrieve the data
dataOut <- dbGetQuery(mycon, sql)
dataOut


###################################################
### code chunk number 8: useHomology.Rnw:250-251
###################################################
toLatex(sessionInfo())


