## unit tests for methods associated with reactome.db
require(hom.Hs.inp.db)
i <- hom.Hs.inp.db

## this tests on ones that I think will always be here (will not have their
## own pkg)
test_cols <- function(){
  res <- columns(i)
  checkTrue(length(res)==100)
  checkTrue("APIS_MELLIFERA" %in% res)
  checkTrue("ANOPHELES_GAMBIAE" %in% res)
  checkTrue("GIARDIA_LAMBLIA" %in% res)
  checkTrue("STRONGYLOCENTROTUS_PURPURATUS" %in% res)
}

test_keytypes <- function(){
  res <- keytypes(i)
  checkTrue(length(res)==100)
  checkTrue("APIS_MELLIFERA" %in% res)
  checkTrue("ANOPHELES_GAMBIAE" %in% res)
  checkTrue("GIARDIA_LAMBLIA" %in% res)
  checkTrue("STRONGYLOCENTROTUS_PURPURATUS" %in% res)
  checkTrue("HOMO_SAPIENS" %in% res)
}

test_keys <- function(){
  res <- head(keys(i, keytype="MUS_MUSCULUS"))
  checkTrue(length(res) > 0)
  checkTrue(is.character(res))
  checkTrue(length(res) == length(unique(res)))
  res2 <- head(keys(i, keytype="HOMO_SAPIENS"))
  checkTrue(length(res2) > 0)
  checkTrue(is.character(res2))
  checkTrue(length(res2) == length(unique(res2)))
}




## Tests ability to get one table/query out.
test_extractWithSimpleInpQuery <- function(){
  table <- "Bos_taurus" ## a table (in this case). (Could also be subquery)  
  k <-  head(keys(i, keytype="BOS_TAURUS"))
  keytype <- "BOS_TAURUS"
  baseSpecies <- AnnotationDbi:::.getBaseSpecies(i)
  baseFiveCode <- AnnotationDbi:::.getBaseFiveCode(baseSpecies)
  checkTrue(baseFiveCode=="HOMSA")
  fiveMap <- AnnotationDbi:::.makeFiveLetterMapping()
  res <- AnnotationDbi:::.extractWithSimpleInpQuery(i, table, k, keytype,
                                                    baseFiveCode, fiveMap)
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==2)
  checkIdentical(c("base.inp_id","BOSTA"), colnames(res))

  ## Then check for next pass usages...
  mergeID <- "base.inp_id"
  mkeytype <- baseSpecies
  mergeKeys <- res[[mergeID]]
  table2 <- "Xenopus_tropicalis"
  res2 <- merge(res,
                AnnotationDbi:::.extractWithSimpleInpQuery(i, table2,
                                                           mergeKeys,
                                                           mkeytype,
                                                           baseFiveCode,
                                                           fiveMap),
                by.x=mergeID, by.y=mergeID,
                all.x=TRUE, all.y=TRUE)

  checkTrue(dim(res2)[1]>0)
  checkTrue(dim(res2)[2]==3)
  checkIdentical(c("base.inp_id","BOSTA","XENTR"), colnames(res2))
}




## Test ability to pull data out in vectorized fashion
test_collateInpQueryResults <- function(){
  ## where keytype and cols are baseSpecies
  tables <- AnnotationDbi:::.UCToStandard(c("HOMO_SAPIENS",
                                            "MUS_MUSCULUS",
                                            "APIS_MELLIFERA"))
  checkIdentical(tables, c("Homo_sapiens","Mus_musculus","Apis_mellifera"))
  k <-  head(keys(i, keytype="HOMO_SAPIENS"), n=6)
  keytype <- "HOMO_SAPIENS"
  fiveMap <- AnnotationDbi:::.makeFiveLetterMapping()
  baseSpecies <- AnnotationDbi:::.getBaseSpecies(i)
  baseFiveCode <- AnnotationDbi:::.getBaseFiveCode(baseSpecies)
  res <- AnnotationDbi:::.collateInpQueryResults(i, tables, keys=k,
                                                 keytype,fiveMap,
                                                 baseFiveCode, baseSpecies)
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("HOMSA","MUSMU","APIME"), colnames(res))

  ## where only the keytype is baseSpecies
  tables <- c("Homo_sapiens","Rattus_norvegicus","Apis_mellifera")
  res2 <- AnnotationDbi:::.collateInpQueryResults(i, tables, keys=k,
                                                 keytype,fiveMap,
                                                 baseFiveCode, baseSpecies)
  checkTrue(dim(res2)[1]>0)
  checkTrue(dim(res2)[2]==3)
  checkIdentical(c("HOMSA","RATNO","APIME"), colnames(res2))

  ## where only a col is baseSpecies
  tables <- c("Mus_musculus","Rattus_norvegicus","Homo_sapiens")
  keytype <- "MUS_MUSCULUS"
  k <-  head(keys(i, keytype="MUS_MUSCULUS"), n=6)
  res3 <- AnnotationDbi:::.collateInpQueryResults(i, tables, keys=k,
                                                 keytype,fiveMap,
                                                 baseFiveCode, baseSpecies)
  checkTrue(dim(res3)[1]>0)
  checkTrue(dim(res3)[2]==4)
  checkIdentical(c("HOMSA","MUSMU","RATNO","HOMSA"), colnames(res3))
  
  ## neither keytype or col is the baseSpecies
  tables <- c("Mus_musculus","Rattus_norvegicus","Apis_mellifera")
  res4 <- AnnotationDbi:::.collateInpQueryResults(i, tables, keys=k,
                                                 keytype,fiveMap,
                                                 baseFiveCode, baseSpecies)
  checkTrue(dim(res4)[1]>0)
  checkTrue(dim(res4)[2]==4)
  checkIdentical(c("HOMSA","MUSMU","RATNO","APIME"), colnames(res4))

}





## and tests for select:
test_select_otherKeytype <- function(){
  k <- head(keys(i, "MUS_MUSCULUS"))
  c <- c("APIS_MELLIFERA","AEDES_AEGYPTI")
  res <-  select(i, keys=k, columns=c, keytype="MUS_MUSCULUS")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("MUS_MUSCULUS","AEDES_AEGYPTI", "APIS_MELLIFERA"),
                 colnames(res))
}

test_select_baseSpeciesKeytype <- function(){
  k <- head(head(keys(i, keytype="HOMO_SAPIENS")))
  c <- c("APIS_MELLIFERA","MUS_MUSCULUS")
  res <-  select(i, keys=k, columns=c, keytype="HOMO_SAPIENS")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("HOMO_SAPIENS","APIS_MELLIFERA","MUS_MUSCULUS"),
                 colnames(res))
}

test_select_baseSpeciesKeytype <- function(){
  k <- head(head(keys(i, keytype="HOMO_SAPIENS")))
  c <- c("APIS_MELLIFERA","HOMO_SAPIENS")
  res <-  select(i, keys=k, columns=c, keytype="HOMO_SAPIENS")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==2)
  checkIdentical(c("HOMO_SAPIENS","APIS_MELLIFERA"),
                 colnames(res))
}


test_select_baseSpeciesCols <- function(){
  k <- head(keys(i, "MUS_MUSCULUS"))
  c <- c("APIS_MELLIFERA","HOMO_SAPIENS")
  res <-  select(i, keys=k, columns=c, keytype="MUS_MUSCULUS")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("MUS_MUSCULUS","APIS_MELLIFERA","HOMO_SAPIENS"),
                 colnames(res))
}




