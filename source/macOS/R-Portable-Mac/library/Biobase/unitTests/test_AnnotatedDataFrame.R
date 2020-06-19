dat <- list(x=factor(1:10), y=I(1:10), z=I(letters[1:10]))
obj1 <- new("AnnotatedDataFrame",
            data=data.frame(x=factor(1:10), y=I(1:10), z=I(letters[1:10]),
              row.names=LETTERS[1:10]),
            varMetadata=data.frame(
              labelDescription=names(dat),
              class=sapply(dat, class), typeof=sapply(dat, typeof),
              mode=sapply(dat, mode),
              row.names=c("x","y","z")))

obj2 <- local({
    obj2 <- obj1
    sampleNames(obj2) <- letters[1:dim(obj1)[[1]]]
    obj2
})

checkAsp <- function(obj1, obj) {
    cidx <- colnames(obj1)
    checkTrue(all(sapply(obj[,cidx, drop=FALSE], typeof)==sapply(obj1, typeof)))
    checkTrue(all(sapply(obj[,cidx, drop=FALSE], class)==sapply(obj1, class)))
    checkTrue(all(sapply(obj[,cidx, drop=FALSE], mode)==sapply(obj1, mode)))
}
checkData <- function(obj1, obj2, obj) {
    checkTrue(all(colnames(obj1) %in% colnames(obj)))
    checkTrue(all(colnames(obj2) %in% colnames(obj)))
    checkAsp(obj1, obj)
    checkAsp(obj2, obj)
}
checkVarMetadata <- function(obj1, obj2, obj) {
    checkTrue(all(colnames(obj1) %in% colnames(obj)))
    checkTrue(all(colnames(obj2) %in% colnames(obj)))
    checkTrue(all(rownames(obj1) %in% rownames(obj)))
    checkTrue(all(rownames(obj2) %in% rownames(obj)))

    checkAsp(obj1, obj)
    checkAsp(obj2, obj)
}
check <- function(obj1, obj2, obj) {
    checkData(pData(obj1), pData(obj2), pData(obj))
    checkVarMetadata(varMetadata(obj1), varMetadata(obj2), varMetadata(obj))
}
checkUnchangedPData <- function(p1, p2, p) {
    checkTrue(all(sapply(colnames(p1), function(nm) {
        identical(p[1:dim(p1)[[1]],nm], p1[,nm])
    })))
    checkTrue(all(sapply(colnames(p2), function(nm) {
        identical(p[dim(p1)[[1]] + 1:dim(p1)[[1]],nm], p2[,nm])
    })))
}

testEmptyCombine <- function() {
    obj <- new("AnnotatedDataFrame")
    checkTrue(identical(obj, combine(obj, obj)))
}

testAnnotatedDataFrameCombine <- function() {
    oldw <- options("warn")
    ## duplicate sampleNames
    checkTrue(identical(obj1, combine(obj1,obj1)))

    ## two distint pData
    obj <- combine(obj1, obj2)
    check(obj1, obj2, obj)
    checkTrue(identical(varMetadata(obj1), varMetadata(obj))) # varMetadata unchanged
    checkTrue(identical(varMetadata(obj2), varMetadata(obj)))

    ## warning about coercing pData factors
    obj2a <- obj2
    pData(obj2a)[["x"]] <- factor(letters[1:10])
    on.exit(options(oldw))
    options(warn=2)
    checkException(combine(obj1,obj2a), silent=TRUE)
    options(oldw)

    ## varMetadata with different numbers of columns
    obj4 <- obj2
    varMetadata(obj4)[,"int"] <- 1:dim(varMetadata(obj1))[[1]]
    varMetadata(obj4)[,"char"] <- I(letters[1:dim(varMetadata(obj1))[[1]]])
    obj <- combine(obj1, obj4)
    checkUnchangedPData(pData(obj1), pData(obj4), pData(obj))
    check(obj1, obj4, obj)
    
    ## varMetadata content mismatch
    obj3 <- obj2
    varMetadata(obj3)[,1] <- varMetadata(obj3)[,4]
    colnames(varMetadata(obj3)) <- colnames(varMetadata(obj2))
    checkException(combine(obj1, obj3), silent=TRUE)

    ## varMetadata multi-column mismatch
    obj3 <- obj2
    varMetadata(obj3)[,1:2] <- varMetadata(obj3)[,4:3]
    colnames(varMetadata(obj3)) <- colnames(varMetadata(obj2))
    checkException(suppressWarnings(combine(obj1, obj3)), silent=TRUE)

    ## varMetadata extra columns
    obj5 <- obj2
    pData(obj5)[,"int"] <- 1:dim(pData(obj1))[[1]]
    pData(obj5)[,"char"] <- I(letters[1:dim(pData(obj1))[[1]]])
    varMetadata(obj5)[c("int","char"),] <- NA
    obj <- combine(obj1, obj5)
    checkUnchangedPData(pData(obj1), pData(obj5), pData(obj))
    check(obj1, obj5, obj)

    ## varMetadata with conflicting information (NAs)
    obj6 <- obj2
    varMetadata(obj6)[2,"typeof"] <- NA
    checkException(suppressWarnings(combine(obj1, obj6)), silent=TRUE)
}

testVarMetadataAssign <- function() {
    ## previously coerced labelData to 'character'
    obj2 <- obj1
    varMetadata(obj2) <- varMetadata(obj1)
    checkTrue(identical(obj1, obj2))

    to <- AnnotatedDataFrame(data.frame(Sample=1:5))
    df <- data.frame(labelDescription="foo", row.names="Sample")
    varMetadata(to) <- df
    checkTrue(validObject(to))

    ## avoid varMetatadata duplication via partial match
    ## https://stat.ethz.ch/pipermail/bioconductor/2014-February/057883.html
    adf <- AnnotatedDataFrame(data.frame(xx=1:5))
    varMetadata(adf)["xx", "labelDescription"] <- "lbl"
    adf$x <- 1:5
    checkIdentical(c("lbl", NA), varMetadata(adf)$labelDescription)
}

testMetadataFactors <- function() {
    pd  = data.frame(covar="Z")
    vmd = data.frame(labelDescription=I("Meta 'covar'"))
    rownames(vmd) = colnames(pd)
    rownames(pd) = "Z"
    a = new("AnnotatedDataFrame", data=pd,  varMetadata=vmd)

    ## factor recode should throw a warning
    oldw=options("warn")
    on.exit(options(oldw))
    options(warn=2)
    pd  = data.frame(covar=LETTERS[1])
    rownames(pd) = LETTERS[1]
    b = new("AnnotatedDataFrame", data=pd,  varMetadata=vmd)
    checkException(combine(a,b), silent=TRUE)
    options(oldw)
}

testNoSharedCols <- function() {
    obj1 <- new("AnnotatedDataFrame",
                data=data.frame(x=factor(1:10)),
                varMetadata=data.frame(labelDescription="x", row.names=c("x")))
    obj2 <- new("AnnotatedDataFrame",
                data=data.frame(y=factor(1:10), row.names=letters[1:10]),
                varMetadata=data.frame(labelDescription="y", row.names=c("y")))
    obj <- combine(obj1,obj2)
    checkTrue(all(pData(obj)[1:10,colnames(pData(obj1)),drop=FALSE]==pData(obj1)))
    checkTrue(all(pData(obj)[11:20,colnames(pData(obj2)),drop=FALSE]==pData(obj2)))
    
}

testPhenoDataFactors <- function() {
    data(sample.ExpressionSet)
    suppressWarnings(obj1 <- updateObject(sample.ExpressionSet))
    obj2 <- obj1
    sampleNames(obj2) <- letters[1:dim(obj1)[[2]]]
    obj1 <- phenoData(obj1)
    obj2 <- phenoData(obj2)
    obj <- combine(obj1, obj2)
    checkTrue(all(pData(obj)[1:nrow(obj1),colnames(pData(obj1)),drop=FALSE]==
                  pData(obj1)))
    checkTrue(all(pData(obj)[nrow(obj1)+1:nrow(obj2),colnames(pData(obj2)),
                             drop=FALSE] == pData(obj2)))
}

testDimLabels <- function() {
    x <- new("AnnotatedDataFrame")
    y <- x
    y@dimLabels <- c("x","y")
    checkException(combine(x,y), silent=TRUE)
}

testNewCovariate <- function() {
    x <- new("AnnotatedDataFrame",data=data.frame(x=1:10))
    x$y <- 1:10
    checkTrue(validObject(x))
    x[["z"]] <- 1:10
    checkTrue(validObject(x))

    x <- new("AnnotatedDataFrame",data=data.frame(x=1:10))
    varMetadata(x)$meta1 <- TRUE
    x[["w"]] <- letters[1:10]
    checkTrue(identical(dim(varMetadata(x)), as.integer(c(2,2))))
    checkTrue(identical(varMetadata(x)["x",,drop=TRUE],
                        list(labelDescription=as.character(NA),meta1=TRUE)))

    x <- new("AnnotatedDataFrame",data=data.frame(x=1:10))
    pData(x) <- pData(x)[1:5,,drop=FALSE]
    checkTrue(validObject(x))
    checkTrue(identical(as.vector(dim(x),"integer"), as.integer(c(5,1))))
}

testReplaceCovariates <- function() {
    ## previously tried to update rather than replace varMetadata
    adf <- new("AnnotatedDataFrame", data=data.frame(x=1:3))
    pData(adf) <- data.frame(y=1:3)
    checkTrue(validObject(adf))
    checkEquals("y", varLabels(adf))

    pData(adf)[,"z"] <- 3:1
    checkTrue(validObject(adf))
    checkEquals(c("y","z"), varLabels(adf))
    checkEquals(data.frame(y=1:3, z=3:1), pData(adf))

    pData(adf)[,"y"] <- NULL
    checkTrue(validObject(adf))
    checkEquals("z", varLabels(adf))
    checkEquals(data.frame(z=3:1), pData(adf))
}

testNewCovariateOnEmptyADF <- function() {
    adf <- new("AnnotatedDataFrame",
               data=data.frame(1:3)[,FALSE,drop=FALSE])
    ## was failing to create varMetadata labelDescription
    pData(adf)$x <- 1:3
    checkTrue(validObject(adf, complete=TRUE))

    obj <- obj1
    pData(obj)[["Z"]] <- NA
    checkTrue(validObject(obj, complete=TRUE))
}

testBadInitializeArugments <- function() {
    checkException(new("AnnotatedDataFrame", data=NULL), silent=TRUE)
    checkException(new("AnnotatedDataFrame", varMetadata=NULL), silent=TRUE)
    checkException(new("AnnotatedDataFrame", data=data.frame(), varMetadata=NULL), silent=TRUE)
}

testNewWithVarMetadata <- function() {
    df <- data.frame(x=1:6,
                     y=rep(c("Low", "High"),3),
                     z=I(LETTERS[1:6]),
                     row.names=paste("Sample", 1:6, sep="_"))
    metaData <- data.frame(labelDescription=c(
                             "Numbers",
                             "Factor levels",
                             "Character"))
    ## standard
    obj <- new("AnnotatedDataFrame",
               data=df, varMetadata=metaData)
    checkTrue(validObject(obj))
    ## varMetadata with inconsistent row names -- silent conversion
    row.names(metaData) <- 1:3
    obj <- new("AnnotatedDataFrame",
               data=df, varMetadata=metaData)
    checkTrue(validObject(obj))
    checkTrue(all(row.names(varMetadata(obj))==names(pData(obj))))
}

testAnnotatedDataFrameFrom <- function() {
    ## empty matrix
    m <- matrix(0,0,0)
    a <- annotatedDataFrameFrom(m, byrow=TRUE)
    checkTrue(validObject(a))
    checkTrue(all.equal(c(0,0), as.vector(dim(a))))
    a <- annotatedDataFrameFrom(m, byrow=FALSE)
    checkTrue(validObject(a))
    checkTrue(all.equal(c(0,0), as.vector(dim(a))))
    ## matrix
    m <- matrix(0,5,10, dimnames=list(letters[1:5], LETTERS[1:10]))
    a <- annotatedDataFrameFrom(m, byrow=TRUE)
    checkIdentical(letters[1:5], sampleNames(a))
    a <- annotatedDataFrameFrom(m, byrow=FALSE)
    checkIdentical(LETTERS[1:10], sampleNames(a))
    ## assayData -- empty env
    ad <- assayDataNew()
    checkTrue(validObject(annotatedDataFrameFrom(ad, byrow=TRUE)))
    checkTrue(validObject(annotatedDataFrameFrom(ad, byrow=FALSE)))
    ## assayData -- empty list
    ad <- assayDataNew(storage.mode="list")
    checkIdentical("list", storageMode(ad))
    checkTrue(validObject(annotatedDataFrameFrom(ad, byrow=TRUE)))
    checkTrue(validObject(annotatedDataFrameFrom(ad, byrow=FALSE)))
    ## assayData -- non-empty env
    ad <- assayDataNew(m=m)
    checkIdentical("lockedEnvironment", storageMode(ad))
    a <- annotatedDataFrameFrom(ad, byrow=TRUE)
    checkIdentical(letters[1:5], sampleNames(a))
    a <- annotatedDataFrameFrom(ad, byrow=FALSE)
    checkIdentical(LETTERS[1:10], sampleNames(a))
    ## assayData -- non-empty list
    ad <- assayDataNew(m=m, storage.mode="list")
    checkIdentical("list", storageMode(ad))
    a <- annotatedDataFrameFrom(ad, byrow=TRUE)
    checkIdentical(letters[1:5], sampleNames(a))
    a <- annotatedDataFrameFrom(ad, byrow=FALSE)
    checkIdentical(LETTERS[1:10], sampleNames(a))
}

testAnnotatedDataFrameDimnames <- function() {
    adf0 <- AnnotatedDataFrame(data.frame(x=1:5, y=5:1, row.names=letters[1:5]),
                              data.frame(foo=1:2, row.names=c("x", "y")))
    adf <- adf0
    dimnames(adf) <- dimnames(adf)
    checkIdentical(adf0, adf)

    adf <- adf0
    exp <- list(LETTERS[seq_len(nrow(adf0))], letters[seq_len(ncol(adf0))])
    dimnames(adf) <- exp
    checkIdentical(exp, dimnames(adf))

    df0 <- varMetadata(adf0)
    rownames(df0) <- exp[[2]]
    checkIdentical(df0, varMetadata(adf))

    adf <- adf0
    dimnames(adf) <- NULL
    checkTrue(validObject(adf))
    checkIdentical(list(as.character(1:5), colnames(adf0)), dimnames(adf))
}
    
