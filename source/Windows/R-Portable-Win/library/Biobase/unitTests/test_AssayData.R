checkAssayDataCombine <- function(nr, nc) {
    obj1 <- assayDataNew(exprs=
                         matrix(runif(nr*nc), nrow=nr, ncol=nc,
                                dimnames=list(
                                  if (nr > 0) letters[1:nr] else NULL,
                                  if (nc > 0) LETTERS[1:nc] else NULL)))
    obj <- combine(obj1,obj1)
    checkTrue(!identical(obj, obj1)) # different environments
    checkTrue(identical(obj1$exprs, obj$exprs))

    storageMode(obj1) <- "list"
    obj <- combine(obj1,obj1)
    checkTrue(identical(obj1, obj)) # same list
    checkTrue(identical(obj1$exprs, obj$exprs))

    if (nc > 2) {
        ## combine distinct cols
        obj1 <- assayDataNew(exprs=
                             matrix(runif(nr*nc), nrow=nr, ncol=nc,
                                    dimnames=list(
                                      if (nr > 0) letters[1:nr] else nr,
                                      LETTERS[1:nc])))
        obj2 <- obj1
        sampleNames(obj2)[3] <- letters[3]
        obj <- combine(obj1, obj2)
        checkTrue(all(dim(obj$exprs)==c(nr,nc+1)))
        checkTrue(identical(obj$exprs[,1:nc],obj1$exprs))
        checkTrue(identical(obj$exprs[,nc+1], obj2$exprs[,3]))
    }

    if (nc > 1) {
        ## inconsistent data -- list, otherwise both copies change!
        storageMode(obj1) <- "list"
        obj2 <- obj1
        obj2$exprs[,1] <- runif(nr)
        checkException(combine(obj1, obj2), silent=TRUE)
    }
}

testAssayDataNew_named_dims <- function()
{
    nms0 <- list(letters[1:5], LETTERS[1:2])
    nms <- Map(setNames, nms0, nms0)
    exprs <- matrix(0, nrow=5, ncol=2, dimnames=nms)
    checkIdentical(nms0[[1]], featureNames(assayDataNew(exprs=exprs)))
    checkIdentical(nms0[[2]], sampleNames(assayDataNew(exprs=exprs)))
}

testAssayDataCombine <- function() {
    checkAssayDataCombine(5,3)
    checkAssayDataCombine(0,0)
    checkAssayDataCombine(1,0)
    checkAssayDataCombine(0,1)
}

testAssayDataCombineRows <- function() {
    m <- matrix(1:20, nrow=5,
                dimnames=list(LETTERS[1:5], letters[1:4]))
    obj <- assayDataNew(exprs=m)
    obj1 <- assayDataNew(exprs=m[1:3,])
    obj2 <- assayDataNew(exprs=m[4:5,])
    checkEquals(obj, combine(obj1, obj2))
    obj3 <- assayDataNew(exprs=m[3:5,])
    checkEquals(obj, combine(obj1, obj3)) # overlapping
}
