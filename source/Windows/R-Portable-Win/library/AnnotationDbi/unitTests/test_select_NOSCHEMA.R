## this will install a testDb stashed in the

## ## this is the package name
## pkgName <- "org.testing.db"

## ## Get the package path
## pkgPath <- system.file("extdata", pkgName, package="AnnotationDbi")

## ## Then install it
## install.packages(pkgPath, repos=NULL)
## and load it
#####install.packages(system.file('extdata','org.testing.db', package='AnnotationDbi'), repos=NULL)


dir.create(testlib <- tempfile())
old_libPaths <- NULL

.setUp <- function()
{
    installed <- rownames(installed.packages(testlib))
    if ("org.testing.db" %in% installed)
        return()
    pkg <- system.file("extdata", "org.testing.db", package="AnnotationDbi")
    suppressPackageStartupMessages(install.packages(
        pkg, lib = testlib, repos=NULL, type="source",
        INSTALL_opts="--no-test-load", verbose = FALSE, quiet = TRUE
    ))
    old_libPaths <<- .libPaths()
    .libPaths(c(testlib, old_libPaths))
    finchCsomes <<- c(as.character(1:15),as.character(17:28),
                     "MT","Un","W","Z","4A","1A","1B")
    finchCols <<- c("CHROMOSOME","SYMBOL","GENENAME","GID","GO","EVIDENCE",
                   "ONTOLOGY","GOALL","EVIDENCEALL","ONTOLOGYALL")
}

.tearDown <- function()
    .libPaths(old_libPaths)

## lower level tests (more useful)
test_keysLow <- function(){
    x <- org.testing.db::org.testing.db
    res <- unique(AnnotationDbi:::.noSchemaKeys(x, "CHROMOSOME"))
    checkTrue(all(sort(res) == sort(finchCsomes)))
}


test_selectLow <- function(){
    x <- org.testing.db::org.testing.db
    keys <- "100008579"
    cols <- "SYMBOL"
    keytype <- "GID"
    res <- AnnotationDbi:::.noSchemaSelect(x, keys, cols, keytype)
    checkTrue(all(res==c("100008579","EGR1")))
    checkTrue(all(colnames(res)==c("GID","SYMBOL")))

    keys <- "brain-derived neurotrophic factor"
    cols <- c("SYMBOL","GID")
    keytype <- "GENENAME"
    res <- AnnotationDbi:::.noSchemaSelect(x, keys, cols, keytype)
    checkTrue(all(res==c("brain-derived neurotrophic factor","BDNF","751584")))
    checkTrue(all(colnames(res)==c("GENENAME","SYMBOL","GID")))

    keys <- "brain-derived neurotrophic factor"
    cols <- c("GO","GID")
    keytype <- "GENENAME"
    res <- head(AnnotationDbi:::.noSchemaSelect(x, keys, cols, keytype),n=1)
    checkTrue(all(res==c("brain-derived neurotrophic factor","GO:0001657",
                    "751584")))
    checkTrue(all(colnames(res)==c("GENENAME","GO","GID")))    
}


## high level tests (does this dispatch right etc.?)
test_columns <- function(){
    x <- org.testing.db::org.testing.db
    res <- columns(x)
    checkTrue(all(sort(res) == sort(finchCols)))
}

test_keytypes <- function(){
    x <- org.testing.db::org.testing.db
    res <- keytypes(x)
    checkTrue(all(sort(res) == sort(finchCols)))
}

test_keys<- function(){                                          ## BOOM
    x <- org.testing.db::org.testing.db
    ## most basic case
    res <- keys(x, "CHROMOSOME")
    checkTrue(all(sort(res) == sort(finchCsomes)))
    
    res <- head(keys(x, "GID"), n=2)
    checkTrue(all(res==c("751582", "751583")))
    
    res <- head(keys(x, "SYMBOL", pattern="BDNF"))
    checkTrue(res=="BDNF")
    
    res <- head(keys(x, "GID", pattern="BDNF", column="SYMBOL"))
    checkTrue(res=="751584")
    
    res <- head(keys(x, "SYMBOL", column="GID"),n=2)
    checkTrue(all(res==c("ACT5C","AHSA2")))
}


test_select <- function(){
    x <- org.testing.db::org.testing.db
    ## most basic case
    res <- select(x, keys="100008579",
                  columns="SYMBOL", keytype="GID")
    checkTrue(all(res==c("100008579","EGR1")))
    checkTrue(all(colnames(res)==c("GID","SYMBOL")))

    ## return more than one column
    res <- select(x, keys="100008579",
                  columns=c("SYMBOL","CHROMOSOME"), keytype="GID")
    checkTrue(all(res==c("100008579","EGR1","13")))
    checkTrue(all(colnames(res)==c("GID","SYMBOL","CHROMOSOME")))

    ## return GO and evidence codes
    suppressWarnings(res <- head(select(x, keys="100008579",
                       columns=c("GO","EVIDENCE"), keytype="GID"),n=1))
    checkTrue(all(res==c("100008579","GO:0000122","IEA")))
    checkTrue(all(colnames(res)==c("GID","GO","EVIDENCE")))

    ## test lookup from alt-key
    res <- select(x, keys="BDNF",
                  columns="GENENAME", keytype="SYMBOL")
    checkTrue(all(res==c("BDNF","brain-derived neurotrophic factor")))
    checkTrue(all(colnames(res)==c("SYMBOL","GENENAME")))
    
}
