library("graph")

data(apopGraph)
data(graphExamples)

test_leaves_undirected <- function() {
    want <- "c"
    checkEquals(want, leaves(graphExamples[[1]]))
}

test_leaves_directed_in <- function() {
    want <- c("trkA", "CASP2", "CASP6", "DNA fragmentation",
              "Nucleus", "CASP9", "TRF3", "CASP10")
    checkTrue(setequal(want, leaves(apopGraph, degree.dir="out")))
}

test_leaves_directed_in <- function() {
    want <- c("TNFa", "TNFb", "FasL", "CD40L", "FAP-1", "NGF",
              "Daxx", "Bax", "CASP4", "Glucocorticoid", "Mtd", "Bad",
              "CASP11", "Bcl-w")
    checkTrue(setequal(want, leaves(apopGraph, degree.dir="in")))
}
