basicCluserGraph <- function() {
    new("clusterGraph", clusters=list(
                          a=c(1,2,3),
                          b=c(4,5,6)))
}

rename_nodes_test <- function() {
    g <- basicCluserGraph()
    checkEquals(as.character(1:6), nodes(g))
    nodes(g) <- letters[1:6]
    checkEquals(letters[1:6], nodes(g))
    checkEquals(letters[1:6], names(edges(g)))
}
