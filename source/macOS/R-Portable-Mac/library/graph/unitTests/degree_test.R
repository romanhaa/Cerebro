library("graph")

data(graphExamples)
data(apopGraph)

test_degree_undirected <- function() {
    g <- graphExamples[[1]]
    want <- as.integer(c(5, 5, 1, 5, 5, 5, 0, 6, 0, 0))
    names(want) <- nodes(g)
    checkEquals(want, degree(g))

    gam <- as(g, "graphAM")
    checkEquals(want, degree(g))
}

test_degree_directed <- function() {
    want_in <- c(TRF1=1L, "NF-kB"=4L, CASP2=2L, Daxx=0L)
    want_out <- c(TRF1=1L, "NF-kB"=1L, CASP2=0L, Daxx=1L)

    got <- degree(apopGraph, c("TRF1", "NF-kB", "CASP2", "Daxx"))
    checkEquals(want_in, got[["inDegree"]])
    checkEquals(want_out, got[["outDegree"]])
}
