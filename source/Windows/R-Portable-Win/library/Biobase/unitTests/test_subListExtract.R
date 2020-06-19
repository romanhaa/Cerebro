library("Biobase")

## simple test list
innerL <- list(i=4L,
               d=1.23,
               l=FALSE,
               n=NULL,
               na=NA,
               s="foo",
               c=3+4i,
               r=charToRaw("g"))
len <- 20L
L <- rep(list(innerL), len)
names(L) <- paste("X", 1:len, sep="")

## unit tests for basic API
test_extract_all_types <- function() {
    for (w in names(innerL)) {
        if (is.null(innerL[[w]]))
          want <- vector(mode="list", length=len)
        else
          want <- as.list(rep(innerL[[w]], len))
        names(want) <- names(L)
        checkEquals(want, subListExtract(L, w))
    }
}

test_simplify_extract_all_types <- function() {
    for (w in names(innerL)) {
        if (is.null(innerL[[w]]))
          next
        else
          want <- rep(innerL[[w]], len)
        names(want) <- names(L)
        checkEquals(want, subListExtract(L, w, simplify=TRUE),
                    msg=paste("extracting", w))
    }
}

test_basic_extract <- function() {
    want <- as.list(rep(4L, len))
    names(want) <- names(L)

    checkEquals(want, subListExtract(L, "i"))
    checkEquals(names(L), names(subListExtract(L, "i")))
    checkEquals(NULL, names(subListExtract(L, "i", keep.names=FALSE)))
}

test_simplify_extract <- function() {
    want <- rep(4L, len)
    names(want) <- names(L)

    checkEquals(want, subListExtract(L, "i", simplify=TRUE))
}

test_simplify_fails_extract <- function() {
    L[[1]] <- list(d=matrix(1:10, ncol=2))

    checkException(subListExtract(L, "i", simplify=TRUE),
                   silent=TRUE)
}

test_extract_bad_inner <- function() {
    L[[3]] <- list(foo=1)

    checkException(subListExtract(L, "i"), silent=TRUE)
}

test_extract_from_empty <- function() {
    checkEquals(list(), subListExtract(list(), "foo"))

    ## but you can't simplify
    checkException(subListExtract(list(), "foo", simplify=TRUE))
}
