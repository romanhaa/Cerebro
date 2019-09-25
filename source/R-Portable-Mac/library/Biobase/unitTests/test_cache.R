test_cache_basic <- function() {
    pre <- "tmp_R_TEST_"
    cfile <- file.path(tempdir(), paste(pre, "theVar.RData", sep=""))
    on.exit(file.remove(cfile), add=TRUE)
    
    expensive <- function(n) n + 1
    cache(theVar <- expensive(5), dir=tempdir(), prefix=pre)
    checkEquals(6, theVar)
    checkTrue(file.exists(cfile))
    ## call again
    remove(theVar)
    cache(theVar <- expensive(5), dir=tempdir(), prefix=pre)
    checkEquals(6, theVar)
}


test_cache_infunc <- function() {
    pre <- "tmp_R_TEST_"
    cfile <- file.path(tempdir(), paste(pre, "theVar.RData", sep=""))
    on.exit(file.remove(cfile), add=TRUE)

    expensive <- function(n) n + 1
    aFunc <- function() {
        m <- 1
        cache(theVar <- expensive(m), dir=tempdir(), prefix=pre)
        theVar
    }
    val <- aFunc()
    checkEquals(2, val)
    checkTrue(file.exists(cfile))
    ## call again
    remove(val)
    val <- aFunc()
    checkEquals(2, val)
}
