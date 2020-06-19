.broadSets <- function()
    getBroadSets(system.file("extdata", "Broad.xml", package="GSEABase"))

.gsc <- function() {
    gs1 <- GeneSet(setName="set1", setIdentifier="id1",
                   geneIds=LETTERS[1:5])
    gs2 <- GeneSet(setName="set2", setIdentifier="id2",
                   geneIds=letters[1:5])
    GeneSetCollection(list(gs1, gs2))
}

test_GSC_list_constructor <- function() {
    gs1 <- GeneSet(setName="123", setIdentifier="456",
                   geneIds=LETTERS[1:5])
    gs2 <- GeneSet(setName="234", setIdentifier="567",
                   geneIds=letters[1:5])
    gsc <- GeneSetCollection(list(gs1, gs2))
    checkTrue(validObject(gsc))
    checkEquals(2, length(gsc))
    checkIdentical(gs1, gsc[[1]])
    checkIdentical(gs2, gsc[[2]])

    ## duplicate entries
    checkException(GeneSetCollection(gs1, gs1), silent=TRUE)
}

test_GSC_list_constructor_strips_names <- function() {
    gsc <- .gsc()
    lst <- list(gsc1=gsc[[1]], gsc2=gsc[[2]])
    checkTrue(is.null(attr(GeneSetCollection(lst),"names")))
}

test_GSC_docs_constructor <- function() {
    gs1 <- GeneSet(setName="123", setIdentifier="456",
                   geneIds=LETTERS[1:5])
    gs2 <- GeneSet(setName="234", setIdentifier="567",
                   geneIds=letters[1:5])
    gsc <- GeneSetCollection(gs1, gs2)
    checkTrue(validObject(gsc))
    checkEquals(2, length(gsc))
    checkIdentical(gs1, gsc[[1]])
    checkIdentical(gs2, gsc[[2]])
}

test_GSC_idAndSetType_constructor <- function() {
    gsc1 <- GeneSetCollection(idType=AnnotationIdentifier("hgu95av2"),
                              setType=KEGGCollection())
    checkEquals(length(reverseSplit(as.list(hgu95av2PATH))),
                length(gsc1))
    i1 <- incidence(gsc1)
    i1 <- i1[sort(rownames(i1)), sort(colnames(i1))]
    gsc2 <- GeneSetCollection(idType=AnnotationIdentifier("hgu95av2.db"),
                              setType=KEGGCollection())
    i2 <- incidence(gsc2)
    i2 <- i2[sort(rownames(i2)), sort(colnames(i2))]
    checkIdentical(i1, i2)
}

test_GSC_ExpressionSet_constructor <- function() {
    data(sample.ExpressionSet)
    gss <- GeneSetCollection(sample.ExpressionSet[200:220,],
                             setType=KEGGCollection())
    checkTrue(is(collectionType(gss[[1]]), "KEGGCollection"))
    checkTrue(is(geneIdType(gss[[1]]), "AnnotationIdentifier"))

    kids <- mget(featureNames(sample.ExpressionSet[200:220,]),
                 hgu95av2PATH)
    kids <- kids[!is.na(kids)]
    ukids <- unique(unlist(kids))
    checkEquals(length(ukids), length(gss))
    checkEquals(table(sapply(reverseSplit(kids), length)),
                table(sapply(lapply(gss, geneIds), length)))

    gss <- GeneSetCollection(sample.ExpressionSet[200:220,],
                             setType=GOCollection())
    checkTrue(is(collectionType(gss[[1]]), "GOCollection"))
    checkTrue(is(geneIdType(gss[[1]]), "AnnotationIdentifier"))

    kids <- mget(featureNames(sample.ExpressionSet[200:220,]),
                 hgu95av2GO)
    kids <- kids[!is.na(kids)]
    ukids <- unique(unlist(sapply(kids, names)))
    checkTrue(all(sort(ukids)==sort(names(gss))))
    rkids <- lapply(kids, lapply, "[[", "GOID")
    checkEquals(table(sapply(reverseSplit(lapply(rkids, unique)),
                             length)),
                table(sapply(lapply(gss, geneIds), length)))
}

test_GSC_validity <- function() {
    gsc <- .gsc()
    gsc@.Data <- append(gsc@.Data, 1)
    checkException(validObject(gsc), silent=TRUE)
}

test_GSC_length <- function() {
    checkTrue(length(.gsc())==2)
}

test_GSC_names <- function() {
    checkTrue(all(c("set1", "set2")==names(.gsc())))
}

test_GSC_subset_by_name<- function() {
    gsc <- .gsc()

    gsc1 <- gsc["set1"]
    checkTrue(is(gsc1, "GeneSetCollection"))
    checkTrue(validObject(gsc1))
    checkEquals(1, length(gsc1))
    checkEquals("set1", names(gsc1))

    gsc21 <- gsc[c("set2", "set1")]
    checkTrue(validObject(gsc21))
    checkEquals(2, length(gsc21))
    checkTrue(all(c("set2", "set1")==names(gsc21)))

    checkException(gsc["set3"], silent=TRUE) # no element
    checkException(gsc[c("set1", "set1")], silent=TRUE) # duplicate entries
}

test_GSC_subset_by_index <- function() {
    gsc <- .gsc()

    gsc1 <- gsc[1]
    checkTrue(is(gsc1, "GeneSetCollection"))
    checkTrue(validObject(gsc1))
    checkEquals(1, length(gsc1))
    checkEquals("set1", names(gsc1))

    gsc21 <- gsc[2:1]
    checkTrue(validObject(gsc21))
    checkEquals(2, length(gsc21))
    checkTrue(all(c("set2", "set1")==names(gsc21)))

    checkException(gsc[3], silent=TRUE) # no element
    checkException(gsc[c(1,1)], silent=TRUE) # duplicate entries
}

test_GSC_subset_by_logical <- function() {
    gsc <- .gsc()
    checkException(gsc[rep(TRUE, 3)], silent=TRUE) # out-of-bounds
}

test_GSC_subset2 <- function() {
    gsc <- .gsc()
    gsc2 <- gsc[[2]]
    checkTrue(is(gsc2, "GeneSet"))
    checkTrue(validObject(gsc2))
    checkTrue("set2"==setName(gsc2))

    gsc2 <- gsc[["set2"]]
    checkTrue(is(gsc2, "GeneSet"))
    checkTrue(validObject(gsc2))
    checkTrue("set2"==setName(gsc2))

    ## subscript out of bounds
    checkException(gsc[[c(1,2)]], silent=TRUE)
    checkException(gsc[[c("set1", "set2")]], silent=TRUE)
    checkException(gsc[[3]], silent=TRUE)
    checkException(gsc[["set3"]], silent=TRUE)
}

## test_GSC_subset_assign <- function() {
##     checkTrue(FALSE)
## }

## test_GSC_subset2_assign <- function() {
##     checkTrue(FALSE)
## }

test_GSC_incidence <- function() {
    gss <- .broadSets()
    res <- incidence(gss)
    checkTrue(all(dim(res)==c(2, 215)))
    checkTrue(sum(res)== 215)
    res <- incidence(gss, gss)
    checkTrue(all(dim(res)==c(4, 215)))
    checkTrue(sum(res)== 430)
}

test_GSC_logic <- function() {
    gsc <- GeneSetCollection(list(GeneSet(letters[1:3], setName="A"),
                                  GeneSet(letters[3:5], setName="B"),
                                  GeneSet(letters[5:7], setName="C")))
    expected <- list(letters[1:3], "c", character(0))
    checkEquals(expected,
                sapply(gsc & geneIds(gsc[[1]]), geneIds))
    checkEquals(expected,
                sapply(gsc & gsc[[1]], geneIds))
    checkEquals(expected,
                sapply(gsc[[1]] & gsc, geneIds))

    expected <- list(letters[1:3],
                     c(letters[3:5], letters[1:2]),
                     c(letters[5:7], letters[1:3]))
    checkEquals(expected, sapply(gsc | geneIds(gsc[[1]]), geneIds))
    checkEquals(expected, sapply(gsc | gsc[[1]], geneIds))
    checkEquals(expected, sapply(gsc[[1]] | gsc, geneIds))
}

test_GSC_mapIdentifiers <- function() {
    data(sample.ExpressionSet)
    gsc <- GeneSetCollection(sample.ExpressionSet[200:205],
                             setType=GOCollection())
    gsc1 <- mapIdentifiers(gsc, EntrezIdentifier())
    checkTrue(is(gsc1, "GeneSetCollection"))
    checkEquals(length(gsc), length(gsc1))
    checkTrue(all(sapply(gsc1, function(x) {
        is(geneIdType(x), "EntrezIdentifier")
    })))
    checkEquals(length(unlist(geneIds(gsc))), length(unlist(geneIds(gsc1))))
}

test_GSC_GOCollection_ontology <- function() {
    idType <- AnnotationIdentifier("org.Hs.eg.db")
    eids <- as.character(1:2)
    setType <- GOCollection()
    gsc <- GeneSetCollection(eids, idType=idType, setType=setType)
    checkIdentical(length(gsc), 31L)
    tbl <- table(unlist(eapply(GOTERM[names(gsc)], Ontology)))
    checkIdentical(as.integer(tbl), c(10L, 9L, 12L))

    setType <- GOCollection(ontology="BP")
    gsc <- GeneSetCollection(eids, idType=idType, setType=setType)
    checkIdentical(length(gsc), 10L)
    tbl <- table(unlist(eapply(GOTERM[names(gsc)], Ontology)))
    checkIdentical(as.integer(tbl), 10L)
}
