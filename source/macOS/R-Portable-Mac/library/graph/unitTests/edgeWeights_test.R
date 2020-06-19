egGraphAM <- function() {
    mat <- matrix(c(0, 0, 1, 1,
                    0, 0, 1, 0,
                    0, 0, 0, 0,
                    0, 1, 1, 0),
                  byrow=TRUE, ncol=4)
    rownames(mat) <- colnames(mat) <- letters[1:4]
    gam <- graphAM(adjMat=mat, edgemode="directed")
}

testDefaultIsOnes <- function() {
    gam <- egGraphAM()
    expect <- list(a=c(c=1, d=1), b=c(c=1), c=numeric(0),
                   d=c(b=1, c=1))
    checkEquals(expect, edgeWeights(gam))
    checkEquals(expect[c("a", "d")], edgeWeights(gam, c("a", "d")))
    checkEquals(expect[c("a", "d")], edgeWeights(gam, c(1, 4)))

    ## Also test alternate attr name when undefined
    checkEquals(expect, edgeWeights(gam, attr="foobar"))
}


testSettingDefaultValue <- function() {
    gam <- egGraphAM()
    expect <- list(a=c(c=4, d=4), b=c(c=4), c=numeric(0),
                   d=c(b=4, c=4))
    checkEquals(expect, edgeWeights(gam, default=4))
    checkEquals(expect[c("a", "d")], edgeWeights(gam, c("a", "d"), default=4))
    checkEquals(expect[c("a", "d")], edgeWeights(gam, c(1, 4), default=4))
}


testTypeChecker <- function() {
    gam <- egGraphAM()
    simple <- list(a=c(c=1, d=1), b=c(c=1), c=numeric(0),
                   d=c(b=1, c=1))
    expect <- list(a=c(c=5:5, d=6:6), b=c(c=4:4), c=numeric(0),
                    d=c(b=4:4, c=4:4))
    edgeDataDefaults(gam, attr="foo") <- 4:4
    edgeData(gam, from="a", attr="foo") <- c(5:5, 6:6)

    checkException(edgeWeights(gam, attr="foo", type.checker=is.double),
                   silent=TRUE)
    checkException(edgeWeights(gam, type.checker=is.integer))
    checkEquals(expect, edgeWeights(gam, attr="foo"))
    checkEquals(expect, edgeWeights(gam, attr="foo", type.checker=is.integer))
    checkEquals(simple, edgeWeights(gam))
}
