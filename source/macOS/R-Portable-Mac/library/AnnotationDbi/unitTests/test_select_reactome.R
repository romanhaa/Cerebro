## unit tests for methods associated with reactome.db
require(reactome.db)
r <- reactome.db


test_cols <- function(){
  res <- columns(r)
  checkTrue(length(res) >4)
  checkTrue("ENTREZID" %in% res)
  checkTrue("GO" %in% res)
  checkTrue("PATHNAME" %in% res)
  checkTrue("REACTOMEID" %in% res)
}

test_keytypes <- function(){
  res <- keytypes(r)
  checkTrue(length(res) >4)
  checkTrue("ENTREZID" %in% res)
  checkTrue("GO" %in% res)
  checkTrue("PATHNAME" %in% res)
  checkTrue("REACTOMEID" %in% res)
}

test_keys <- function(){
  res <- head(keys(r, keytype="ENTREZID"))
  checkTrue(length(res) > 0)
  checkTrue(length(res) == length(unique(res)))
  res2 <- head(keys(r, keytype="PATHNAME"))
  checkTrue(is.character(res2))
  checkTrue(length(res2) == length(unique(res2)))
}

## test function that gets table names
test_getTables <- function(){
  c <- c("ENTREZID", "PATHNAME")
  res <- AnnotationDbi:::.getTables(c, retVal="table") ##default for retVal
  checkTrue(length(res) ==2)
  checkIdentical(res, c("pathway2gene","pathway2name"))
  
  res2 <- AnnotationDbi:::.getTables(c, retVal="colname")
  checkTrue(length(res2) ==2)
  checkIdentical(res2, c("gene_id","path_name"))
}

## Tests ability to get one table/query out.
test_extractWithSimpleQuery <- function(){
  table <- "pathway2gene" ## a table (in this case). (Could also be subquery)
  colType <- "gene_id" ## column we are interested in.
  k <-  head(keys(r, keytype="ENTREZID"), n=2)
  res <- AnnotationDbi:::.extractWithSimpleQuery(r, table, colType, k)
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==2)
  checkIdentical(c("DB_ID","gene_id"), colnames(res))
}

## Test ability to pull data out in vectorized fashion
test_collateQueryResults <- function(){
  tables <- c("pathway2gene", "reactome2go", "pathway2name")
  colType <- "gene_id"
  k <-  head(keys(r, keytype="ENTREZID"), n=2)
  mergeID = "DB_ID"
  res <- AnnotationDbi:::.collateQueryResults(r, tables, colType, k, mergeID)
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==4)
  checkIdentical(c("DB_ID","gene_id","go_id","path_name"), colnames(res))
}





## and tests for select:
test_select_TYPICAL <- function(){
  k <- head(keys(r, keytype="ENTREZID"))
  c <- c("ENTREZID","PATHNAME","GO","REACTOMEID")
  res <-  select(r, keys=k, columns=c, keytype="ENTREZID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==4)
  checkIdentical(c("ENTREZID","PATHNAME","GO","REACTOMEID"), colnames(res))
}

test_select_MISSING_EG <- function(){
  k <- head(keys(r, keytype="ENTREZID"))
  c <- c("PATHNAME","GO")
  res <-  select(r, keys=k, columns=c, keytype="ENTREZID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("ENTREZID","PATHNAME","GO"), colnames(res))
}


test_select_ONE_COL <- function(){
  k <- head(keys(r, keytype="ENTREZID"))
  c <- c("ENTREZID")
  res <-  select(r, keys=k, columns=c, keytype="ENTREZID")  ## Boom no warning
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==1)
  checkIdentical(c("ENTREZID"), colnames(res))
}


test_select_OTHERKEYTYPE <- function(){
  ## This also checks if that we handle "BS keys" OK
  k <- head(keys(r, keytype="REACTOMEID"))
  k <- c(k[1:4], "109688")
  c <- c("ENTREZID","PATHNAME","GO")
  res <-  select(r, keys=k, columns=c, keytype="REACTOMEID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==4)
  checkIdentical(c("REACTOMEID","ENTREZID","PATHNAME","GO"), colnames(res))  
}

test_select_PATHNAME <- function(){
  k <- head(keys(r, "PATHNAME"))
  suppressWarnings(res <- select(r, k, "ENTREZID", "PATHNAME"))
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==2)
  checkIdentical(c("PATHNAME","ENTREZID"), colnames(res)) 
}




## for convenience...
## debug(AnnotationDbi:::.selectReact)
## debug(AnnotationDbi:::.collateQueryResults)
## debug(AnnotationDbi:::.extractWithSimpleQuery)
