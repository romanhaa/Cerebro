testSnpSetAccessors <-
    function()
{
    set.seed(1)
    nms <- list(paste("rs", 1:100, sep="_"), LETTERS[1:10])
    m <- matrix(runif(1000), 100, dimnames=nms)
    n <- matrix(runif(1000), 100, dimnames=nms)
    sset <- new("SnpSet", call=m, callProbability=n)

    checkEquals(m, snpCall(sset))
    checkEquals(n, snpCallProbability(sset))

    m <- matrix(runif(1000), 100, dimnames=nms)
    n <- matrix(runif(1000), 100, dimnames=nms)
    snpCall(sset) <- m
    snpCallProbability(sset) <- n
    checkEquals(m, snpCall(sset))
    checkEquals(n, snpCallProbability(sset))
}
