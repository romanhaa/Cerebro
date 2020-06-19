## One strategy is to put some of my examples (with expected output into here)
## - but this is fragile unless I also provide a localized example DB (which I
## probably should do anyways for other reasons, but is not necessarily useful
## here).
## Another strategy is that I should probably have tests for all my helper
## functions to make sure that they are returning what is expected.
## Generally, I want to make tests for each thing that can go wrong due to
## changes elsewhere.  The strategy of writing tests for the helpers will
## catch some of this, but also I need to anticipate things that will change
## in the annotations etc.

##  library(AnnotationDbi);AnnotationDbi:::.test()
require(RSQLite)
require(org.Hs.eg.db)
require(org.At.tair.db)
require(org.Sc.sgd.db)
require(GO.db)
require(hgu95av2.db)
require("RUnit")
x <- org.Hs.eg.db
t <- org.At.tair.db
s <- org.Sc.sgd.db
cols <- c("CHR","PFAM","GO")
keys <- c(1,10)
jointype <- "genes.gene_id" ## changed from 'gene_id'
quiet <- suppressWarnings # quieten warnings from 1:many mappings in select()

## resort and friends are really important as they are generic enough to
## be reused elsewhere.
test_generateExtraRows <- function(){
  ttab = data.frame(warpbreaks[1:10,])
  ttab$breaks <- as.character(ttab$breaks)  
  tkeys = ttab$breaks
  tkeys = c(26, tkeys[1:7], tkeys[7], 30, tkeys[8:10], tkeys[10])
  res <- AnnotationDbi:::.generateExtraRows(ttab, tkeys, jointype='breaks')
  expLen <- sum(table(tkeys) * table(ttab$breaks))
  checkTrue(expLen == dim(res)[1])
}

test_dropUnwantedRows <- function() {
    fun <- AnnotationDbi:::.dropUnwantedRows

    ## no duplicates, no changes
    keys <- letters[1:5]
    tab <- data.frame(x=keys, y=LETTERS[1:5], z=LETTERS[5:1],
                      row.names=NULL)
    checkIdentical(tab, fun(tab, keys, "x"))

    ## duplicate jointype element, duplicate dropped
    tab1 <- tab[c(1:5, 3L),]
    rownames(tab1) <- NULL
    checkIdentical(tab, fun(tab1, keys, "x"))

    ## unique all NA (other than jointype column) _retained_
    tab1 <- tab
    tab1[3, 2:3] <- NA
    rownames(tab1) <- NULL
    checkIdentical(tab1, fun(tab1, keys, "x"))

    ## duplicate all NA, made unique
    tab1 <- tab
    tab1[3, 2:3] <- NA
    tab2 <- tab1[c(1:5, 3),]
    checkIdentical(tab1,  fun(tab2, keys, "x"))

    ## duplicate key, dropped
    keys1 <- keys[c(1:5, 3)]
    checkIdentical(tab, fun(tab, keys1, "x"))
}

test_resort <- function() {
    fun <- resort_base ## from AnnotationDbi

    ## repeat keys returned
    keys <- letters[1:5]
    tab <- data.frame(x=keys, y=LETTERS[1:5], z=LETTERS[5:1],
                      row.names=NULL, stringsAsFactors=FALSE)
    keys1 <- keys[c(1:5, 1)]
    tab1 <- tab[c(1:5, 1),]
    rownames(tab1) <- NULL
    checkIdentical(tab1, fun(tab, keys1, "x", names(tab)))

    ## keys with missing values returned
    tab1 <- tab
    tab1[3, 2:3] <- NA
    keys1 <- tab1[["x"]]
    checkIdentical(tab1, fun(tab1, keys, "x", names(tab)))

    ## multiple keys with missing values returned
    tab1 <- tab[c(3,4,3,4),]
    tab1[c(1,3), 2:3] <- NA
    keys1 <- keys[c(3,4,3,4)] 
    rownames(tab1) <- NULL
    checkIdentical(tab1, fun(tab1[1:2,], keys1, "x", names(tab)))

    cols <- c("CHR","SYMBOL", "PFAM")
    keys <- c(1,10)
    res <- AnnotationDbi:::.extractData(x, cols, keytype="ENTREZID", keys)
    ## jumble res to simulate trouble
    resRO = res[order(sort(res$genes.gene_id,decreasing=TRUE)),]
    reqCols <- c("genes.gene_id","chromosomes.chromosome","gene_info.symbol",
                 "pfam.pfam_id")
    Rres <- fun(resRO, keys, jointype, reqCols)
    checkIdentical(Rres$gene_id,Rres$gene_id)
    checkTrue(class(Rres) =="data.frame")

    ## now what if we have MORE keys?
    keys <- c(1, keys, keys)
    cols <- c("CHR","SYMBOL")
    res <- AnnotationDbi:::.extractData(x, cols, keytype="ENTREZID", keys)
    reqCols <- c("genes.gene_id","chromosomes.chromosome","gene_info.symbol")
    res2 <- fun(res, keys, jointype, reqCols)
    checkIdentical(as.numeric(as.character(res2$genes.gene_id)),keys)
    checkTrue(class(res) =="data.frame")
}

test_keytypes <- function(){
  checkTrue("ENTREZID" %in% keytypes(x))
  checkTrue("TAIR" %in% keytypes(t))
  checkTrue("ENTREZID" %in% keytypes(t))
  checkTrue("ORF" %in% keytypes(s))  
  checkTrue("ENTREZID" %in% keytypes(s))  
}

test_keys <- function(){
  checkException(keys(org.Hs.eg.db, keytype="PROBEID"))
  
  egHskeys <- as.numeric(head(keys(x)))
  checkTrue(length(egHskeys[!is.na(egHskeys)])==6)
  rsHskeys <- head(keys(x, "REFSEQ"))
  checkTrue(any(grepl("N", rsHskeys)))
  
  egAtkeys <- as.numeric(head(keys(t,"ENTREZID")))
  checkTrue(length(egAtkeys[!is.na(egAtkeys)])==6)
  rsAtkeys <- head(keys(t, "REFSEQ"))
  checkTrue(any(grepl("N", rsAtkeys)))
  tairAtkeys <- head(keys(t, "TAIR"))
  checkTrue(any(grepl("AT", tairAtkeys)))

  egSckeys <- as.numeric(head(keys(s, "ENTREZID")))
  checkTrue(length(egSckeys[!is.na(egSckeys)])==6)
  rsSckeys <- head(keys(s, "REFSEQ"))
  checkTrue(any(grepl("N", rsSckeys)))
  orfSckeys <- head(keys(s, "ORF"))
  checkTrue(any(grepl("A", orfSckeys)))
}

test_keys_advancedArgs <- function(){
    k1 <- keys(x, keytype="SYMBOL")
    checkTrue("A1BG" %in% k1)
    
    k2 <- keys(x, keytype="SYMBOL", pattern="BRCA")
    checkTrue("BRCA1" %in% k2)
    checkTrue(!("A1BG" %in% k2))
    checkTrue(length(k2) < length(k1))

    l1 <- length(keys(x, keytype="ENTREZID", column="PATH"))
    l2 <- length(keys(x, keytype="ENTREZID"))
    checkTrue(l1 < l2)
    
    k3 <- keys(x,keytype="ENTREZID",pattern="^MSX",column="SYMBOL")
    res <- select(x, k3, c("ENTREZID","SYMBOL"), "ENTREZID")
    checkTrue(any(grep("^MSX",res$SYMBOL)))
}

#########################################################################
## These ones are to test out some real use cases...
test_select1 <- function(){
  keys <- head(keys(hgu95av2.db, "ALIAS"),n=2)
  cols <- c("SYMBOL","ENTREZID","PROBEID")
  res <- quiet(select(hgu95av2.db, keys, cols, keytype="ALIAS"))
  checkIdentical(c(3L, 4L), dim(res))
  checkIdentical(c("ALIAS","SYMBOL","ENTREZID","PROBEID"), colnames(res))
}

test_select2 <- function(){
  keys <- head(keys(org.Hs.eg.db),n=2)
  cols <- c("PFAM","ENTREZID", "GO")
  res <- quiet(select(org.Hs.eg.db, keys, cols, keytype="ENTREZID"))
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==5)
  checkIdentical(c("ENTREZID","PFAM","GO","EVIDENCE","ONTOLOGY"),
                 colnames(res))
}

test_select3 <- function(){
  keys <- head(keys(org.Hs.eg.db,keytype="OMIM"),n=4)
  cols <- c("SYMBOL", "UNIPROT", "PATH")
  res <- quiet(select(hgu95av2.db, keys, cols, keytype="OMIM"))
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==4)
  checkIdentical(c("OMIM","SYMBOL","UNIPROT","PATH"), colnames(res))
}

test_select4 <- function(){
  keys <- head(keys(org.Hs.eg.db),n=2)
  cols <- c("ACCNUM","REFSEQ")
  res <- quiet(select(org.Hs.eg.db, keys, cols))
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("ENTREZID","ACCNUM","REFSEQ"), colnames(res))
}

test_select5 <- function(){
  keys <- head(keys(GO.db), n=4)
  cols <- c("TERM","ONTOLOGY","DEFINITION")
  res <- select(GO.db, keys, cols)
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==4)
  checkIdentical(c("GOID","TERM","ONTOLOGY","DEFINITION"), colnames(res))
}

test_select6 <- function(){
  keys <- head(keys(hgu95av2.db))
  cols <- c("SYMBOL","ENTREZID", "GO")
  ## tests for bad keys:
  checkException(select(hgu95av2.db, keys, cols, keytype="ENTREZID"))
  ## also catch bogus keytype arguments
  checkException(select(hgu95av2.db, keys, cols, keytype="FOO"))
  checkException(keys(hgu95av2.db, keytype="FOO"))
}

test_select7 <- function(){  
  cols <- c("SYMBOL","ENTREZID") ## 1st of all cols should be 1:1 cols
  keys <- head(keys(org.Hs.eg.db),n=3)
  keys <- c(1, keys, keys)
  res <- select(org.Hs.eg.db, keys, cols)
  checkTrue(class(res) =="data.frame")
  checkIdentical(keys, as.character(t(res$ENTREZID)))
}

test_select8 <- function(){
  cols <- c("ENTREZID")
  keys <- head(keys(org.Hs.eg.db),n=3)
  res <- select(org.Hs.eg.db, keys, cols)
  checkTrue(class(res) =="data.frame")
  checkTrue(dim(res)[2] ==1)  
  checkIdentical(as.character(keys), as.character(t(res$ENTREZID)))
}

test_select9 <- function(){  
  ## What about when we need to throw away extra cols?
  uniKeys <- head(keys(org.Hs.eg.db, keytype="UNIPROT"))
  cols <- c("SYMBOL", "PATH")
  res <- quiet(select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="UNIPROT"))
  checkTrue(class(res) =="data.frame")
  checkTrue(dim(res)[2] ==3)  
  checkIdentical(c("UNIPROT","SYMBOL","PATH"), colnames(res))
}

test_select10 <- function(){
  ## What about when we have to get data from Arabidopsis using various
  ## keytypes?
  cols <- c("SYMBOL","CHR")
  keys <- head(keys(t,"TAIR"))
  res <- quiet(select(t, keys, cols, keytype="TAIR"))
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("TAIR","SYMBOL","CHR"), colnames(res))

  keys <- head(keys(t,"ENTREZID"))
  res <- quiet(select(t, keys, cols, keytype="ENTREZID"))
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("ENTREZID","SYMBOL","CHR"), colnames(res))

  keys=head(keys(t,"REFSEQ"))
  res <- quiet(select(t, keys, cols , keytype="REFSEQ"))
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("REFSEQ","SYMBOL","CHR"), colnames(res))
}

test_select11 <- function(){
  ## how about different keytypes for yeast?
  keys <- head(keys(s, "REFSEQ"))
  cols <- c("CHR","PFAM")
  res <- quiet(select(s, keys, cols, keytype="REFSEQ"))
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("REFSEQ","CHR","PFAM"), colnames(res))
  
  keys <- head(keys(s, "ENTREZID"))
  cols <- c("CHR","PATH")
  res <- quiet(select(s, keys, cols, keytype="ENTREZID"))
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("ENTREZID","CHR","PATH"), colnames(res))
  
  keys <- head(keys(s, "ORF"))
  cols <- c("CHR","SGD")
  res <- select(s, keys, cols, keytype="ORF")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("ORF","CHR","SGD"), colnames(res))

  ## And if you flip things the other way
  cols <- c("SGD","CHR")
  res <- select(s, keys, cols, keytype="ORF")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("ORF","SGD","CHR"), colnames(res))

  ## Martins bug discoveries
  keys <- keys(s, keytype="GENENAME")
  checkTrue(length(keys) > 0)
  checkTrue(is.character(keys))
  keys <- keys(s, keytype="CHRLOC")
  checkTrue(length(keys) > 0)
  checkTrue(is.character(keys))

  res <- select(s, "YAL003W", "GENENAME")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("ORF","SGD","GENENAME"), colnames(res))

  ## This works but is slow (therefore it's tested elsewhere)
  ## res <- select(s, keys="YAL003W", columns(s))

  ## Another test to make sure we can join up to ORF properly
  keys <- keys(s,"ENTREZID")
  res <- select(s, columns="ORF", keys=keys, keytype="ENTREZID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("ENTREZID","ORF","SGD"), colnames(res))
}

test_select12 <- function(){
  ## what happens when we use GO as an ID?
  keys <- "1"
  cols <- c("GO","ENTREZID")
  res <- quiet(select(x, keys, cols, keytype="ENTREZID"))
  checkTrue(dim(res)[1]>0)   
  checkTrue(dim(res)[2]==4)
  checkIdentical(c("ENTREZID","GO","EVIDENCE","ONTOLOGY"), colnames(res))

  keys <- "GO:0000018"
  cols <- c("GO","ENTREZID")
  res <- quiet(select(x, keys, cols, keytype="GO"))
  checkTrue(dim(res)[1]>0)   
  checkTrue(dim(res)[2]==4)
  checkIdentical(c("GO","EVIDENCE","ONTOLOGY","ENTREZID"), colnames(res))

  keys <- "GO:0000023"
  cols <- c("GO","ENTREZID")
  res <- quiet(select(t, keys, cols, keytype="GO"))
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==4)
  checkIdentical(c("GO","EVIDENCE","ONTOLOGY","ENTREZID"), colnames(res)) 

  keys <- "GO:0000023"
  cols <- c("ENTREZID","TAIR","GO")
  res <- quiet(select(t, keys, cols, keytype="GO"))
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==5)
  checkIdentical(c("GO","EVIDENCE","ONTOLOGY","ENTREZID","TAIR"), 
  		colnames(res))
}

test_select13 <- function(){
  ## what happens with dropping unwanted rows?
  sym <- "ITGA7"
  res <- quiet(select(x, sym, "PFAM", keytype="ALIAS"))
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==2)
  ## make sure no NAs are in res$PFAM
  checkTrue(length(res$PFAM)== length(res$PFAM[!is.na(res$PFAM)]))
}

test_select14 <- function(){
  ## what happens when there are no results AT ALL? (should be all NAs)
  keys <- c("1001_at","1006_at","1007_s_at")
  res <- select(hgu95av2.db, keys, "PATH", keytype="PROBEID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==2)
  ## make sure all of res$PATH ARE NAs
  ## If this part fails it is a warning that the test is no longer valid,
  ## which would happen if some of these IDs were to be further annotated for
  ## PATH (unlikely since PATH is basically dead for this repos)
  checkTrue(length(res$PATH)== length(res$PATH[is.na(res$PATH)]))
}

test_select15 <- function(){
    ## Another bug that seems to happen in post-processing...
    ## the code that resolves duplicated values is going a bit insane...
    ## (IOW .replaceValues())
    if(!all(.Platform$OS.type == "windows", .Platform$r_arch == "i386")){
        res <- select(x, keys="100008586", columns(x)) 
        checkTrue(dim(res)[1]>0)
        checkTrue(dim(res)[2]==26)
        exp <- c("ENTREZID", "ACCNUM", "ALIAS", "ENSEMBL", "ENSEMBLPROT",
                 "ENSEMBLTRANS", "ENZYME", "EVIDENCE", "EVIDENCEALL",
                 "GENENAME", "GO", "GOALL", "IPI", "MAP", "OMIM", 
                 "ONTOLOGY", "ONTOLOGYALL", "PATH", "PFAM", "PMID", "PROSITE",
                 "REFSEQ", "SYMBOL", "UCSCKG", "UNIGENE", "UNIPROT")
        checkIdentical(exp, colnames(res))
    }
}


test_select16 <- function(){
    ## What happens if we ask for probes back...
    ## (and pass in something else as a key)
    sk = c( 'MAPK3','TIE1' )
    res <- select(hgu95av2.db, keys=sk, columns = c("PROBEID"), keytype="SYMBOL")
    checkTrue(dim(res)[1]>0)
    checkTrue(dim(res)[2]==2)
    checkIdentical(c('SYMBOL','PROBEID'), colnames(res))
}


## NA values are now OK for legacySelect (in order to be consistent
## with other select methods which have greater respect for incoming NA values
test_select_NAsInNAsOut <- function(){
    ## NAs in should result in NAs out.
    ## Not that we like NAs but just because if possible: we want to
    ## preserve the geometry of the keys coming in.
    k=c('1', NA, NA, '10');
    res <- select(x, k, 'SYMBOL', 'ENTREZID')
    checkTrue(dim(res)[1]==4) 
    checkTrue(dim(res)[2]==2)
    checkIdentical(c('ENTREZID','SYMBOL'), colnames(res))
    checkIdentical(k, res$ENTREZID)
    checkTrue(any(is.na(res$SYMBOL)))
}


## Some new messages from AnnotationDbi:::.generateExtraRows()
## where I am calling it internally...
## these tests just make sure the right kinds of messages are being sent...
test_select_XtoX <- function(){
    k=c('1','1', NA, '10');
    res <- tryCatch(select(org.Hs.eg.db, k, 'SYMBOL', 'ENTREZID'),
                    message = function(x){return(x)})
    checkTrue(grepl('many:1', res$message))
    
    res <- tryCatch(select(org.Hs.eg.db, '1', 'ALIAS', 'ENTREZID'),
                    message = function(x){return(x)})
    checkTrue(grepl('1:many', res$message))

    res <- tryCatch(select(org.Hs.eg.db, k, 'ALIAS', 'ENTREZID'),
                    message = function(x){return(x)})
    checkTrue(grepl('many:many', res$message))

    res <- tryCatch(select(org.Hs.eg.db, c('1','10'), 'SYMBOL', 'ENTREZID'),
                    message = function(x){return(x)})
    checkTrue(grepl('1:1', res$message))
}


## TODO: deal with any fallout from having new messages in select()
## This will also cause me to have to silence select() in many places
## to avoid annoying repeats of the



test_dbconn_and_dbfile <- function(){
    resc <- dbconn(org.Hs.eg.db)
    m <- dbGetQuery(resc, "SELECT * FROM metadata")
    checkTrue(dim(m)[2] ==2)
    checkTrue(dim(m)[1] > 10)
              
    resf <- dbfile(org.Hs.eg.db)
    mf <- dbGetQuery(dbConnect(SQLite(), resf), "SELECT * FROM metadata")
    checkTrue(all(mf == m))
}




## Fast checking:
## BiocGenerics:::testPackage(pattern="^test_select.*\\.R$")
