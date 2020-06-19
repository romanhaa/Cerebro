testCombineFeatureData <- function() {
    data(sample.ExpressionSet)          # use as source for exprs data
    suppressWarnings(obj <- updateObject(sample.ExpressionSet)[1:20,1:10])

    obj1 <- new("ExpressionSet", phenoData=phenoData(obj), exprs=exprs(obj))
    obj2 <- obj1

    pData(featureData(obj1))[["x"]] <- FALSE
    pData(featureData(obj1))[["y"]] <- FALSE
    varMetadata(featureData(obj1)) <-
      data.frame(labelDescription=c("the x", "the y"), row.names=c("x", "y"))
    validObject(obj1)

    sampleNames(obj2) <- letters[1:dim(obj1)[[2]]]
    pData(featureData(obj2))[["y"]] <- FALSE
    pData(featureData(obj2))[["z"]] <- TRUE
    varMetadata(featureData(obj2)) <-
        data.frame(labelDescription=c("the y", "the z"), row.names=c("y", "z"))
    validObject(obj2)
    obj <- combine(obj1,obj2)
    checkTrue(all(varLabels(featureData(obj1)) %in% varLabels(featureData(obj))))
    checkTrue(all(varLabels(featureData(obj2)) %in% varLabels(featureData(obj))))

    ## conflicting feature pData
    pData(featureData(obj2))[["y"]] <- TRUE
    validObject(obj2)
    checkException(suppressWarnings(combine(obj1, obj2)), silent=TRUE)
}

testCombineRows <- function() {
    data(sample.ExpressionSet)
    obj <- sample.ExpressionSet

    checkEquals(obj, combine(obj[1:250,], obj[251:500,]))
    checkEquals(obj, combine(obj[,1:13], obj[,14:26]))
    ## overlapping
    checkEquals(obj, combine(obj[1:300,], obj[250:500,]))
    checkEquals(obj, combine(obj[,1:20], obj[,15:26]))
}

testAddTextNotes <- function() {
    eset <- new("ExpressionSet")
    notes(eset) <- "a note"
    checkTrue(identical(notes(eset), list("a note")))
    notes(eset) <- "another"
    checkTrue(identical(notes(eset), list("a note", "another")))
}

testExtraSlotExpressionClassInitialize1 <- function() {
    setClass("ExtraSlotExpressionSet", contains="ExpressionSet",
             representation=representation(
               extraSlot="character"),
             where=.GlobalEnv)
    ## pass if no error
    checkTrue(validObject(new("ExtraSlotExpressionSet")))
    removeClass("ExtraSlotExpressionSet", where=.GlobalEnv)
}

testExtraSlotExpressionClassInitialize2 <- function() {
    setClass("ExtraSlotExpressionSet", contains="ExpressionSet",
             representation=representation(
               extraSlot="character"),
             where=.GlobalEnv)
    e <- new("ExtraSlotExpressionSet",
             exprs=new("matrix"),
             extraSlot="hello",
             storage.mode="environment")
    checkEquals("hello", e@extraSlot)
    checkEquals("exprs", ls(assayData(e)))
    checkEquals("environment", storageMode(e))
    removeClass("ExtraSlotExpressionSet", where=.GlobalEnv)
}

testExtraSlotExpressionClassInitialize3 <- function() {
    setClass("ExtraSlotExpressionSet", contains="ExpressionSet",
             representation=representation(
               extraSlot="character"),
             where=.GlobalEnv)
    e <- new("ExtraSlotExpressionSet",
             assayData=assayDataNew(
               exprs=new("matrix"),
               storage.mode="environment"),
             extraSlot="hello")
    checkEquals("hello", e@extraSlot)
    checkEquals("exprs", ls(assayData(e)))
    checkEquals("environment", storageMode(e))
    removeClass("ExtraSlotExpressionSet", where=.GlobalEnv)
}

testDollar <- function() {
    data(sample.ExpressionSet)
    s1 <- sample.ExpressionSet$sex
    s2 <- sample.ExpressionSet$se       # we expect partial matching to work
    checkTrue(!is.null(s1), msg="$sex broken")
    checkTrue(!is.null(s2), msg="$se broken (pmatch)")
    checkEquals(s1, s2, msg="pmatch equality")
}

testSubset2 <- function() {
    data(sample.ExpressionSet)
    es <- sample.ExpressionSet
    x <- runif(ncol(es))
    ldesc <- "Random variate"
    es[["RVar", labelDescription=ldesc]] <- x
    checkEquals(es[["RVar"]], x)
    checkEquals(varMetadata(es)["RVar", "labelDescription"], ldesc)
}

testHarmonizeAssayDataDimnames <- function() {
    checkHarmonizeOne <- function(exprs) {
        es <- new("ExpressionSet", exprs=exprs)
        checkTrue(validObject(es))
    }
    checkHarmonizeTwo <- function (exprs, se.exprs) {
        es <- new("ExpressionSet", exprs=exprs, se.exprs=se.exprs)
        checkTrue(validObject(es))
        okNames <- list(featureNames(featureData(es)),
                        sampleNames(phenoData(es)))
        dimNames <- Biobase:::.assayDataDimnames(assayData(es))
        checkTrue(all(sapply(dimNames, identical, okNames)))
    }
    checkCreation <- function(exprs, se.exprs) {
        checkHarmonizeOne(exprs)
        checkHarmonizeTwo(exprs, se.exprs)

        ## names on both dimnames
        nexprs <- exprs
        dimnames(nexprs) <-
            lapply(dimnames(nexprs), function(x) {
            names(x) <- as.vector(letters[seq(1, length(x))])
            x
        })        
        checkHarmonizeOne(nexprs)
        checkHarmonizeTwo(nexprs, se.exprs)

        ## names on colnames
        nexprs <- exprs
        cnms <- colnames(nexprs)
        names(cnms) <- letters[seq(1, length(cnms))]
        colnames(nexprs) <- cnms
        checkHarmonizeOne(nexprs)
        checkHarmonizeTwo(nexprs, se.exprs)

        ## names on rownames
        nexprs <- exprs
        rnms <- rownames(nexprs)
        names(rnms) <- letters[seq(1, length(rnms))]
        rownames(nexprs) <- rnms
        checkHarmonizeOne(nexprs)
        checkHarmonizeTwo(nexprs, se.exprs)
    }

    se.exprs <- matrix(0, 5, 2)
    exprs <- matrix(0, 5, 2)
    dimnames(exprs) <- list(LETTERS[1:5], letters[1:2])

    dimnames(se.exprs) <- NULL
    checkCreation(exprs, se.exprs)

    dimnames(se.exprs) <- list(LETTERS[1:5], NULL)
    checkCreation(exprs, se.exprs)

    dimnames(se.exprs) <- list(NULL, letters[1:2])
    checkCreation(exprs, se.exprs)

    ## errors
    dimnames(se.exprs) <- list(letters[1:5], letters[1:2])
    checkException(checkCreation(exprs, se.exprs), silent=TRUE)

    dimnames(se.exprs) <- list(LETTERS[1:5], LETTERS[1:2])
    checkException(checkCreation(exprs, se.exprs), silent=TRUE)

    dimnames(se.exprs) <- list(letters[1:5], LETTERS[1:2])
    checkException(checkCreation(exprs, se.exprs), silent=TRUE)
}

testExprsReplacement <- function() {
    exprs <- se.exprs <- matrix(1:50, 10, 5)
    eset <- ExpressionSet(list2env(list(exprs=exprs, se.exprs=se.exprs)))
    exprs(eset) <- exprs(eset)
    checkTrue(validObject(eset))

    ## shuffled names ok
    exprs(eset) <- exprs(eset)[sample(rownames(eset)), sample(colnames(eset))]
    checkTrue(validObject(eset))

    checkException({ exprs(eset) <- exprs(eset)[, 1:3] }, silent=TRUE)
    checkException({ exprs(eset) <- exprs(eset)[, c(1:4, 1)] }, silent=TRUE)
}
