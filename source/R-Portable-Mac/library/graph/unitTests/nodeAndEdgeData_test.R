#
# Test setup
#
simpleInciMat <- function() {
    ## Here's a simple graph for testing
    ##    a           b
    ##    |\         /|
    ##    | \___c___/ |
    ##    |     |     |
    ##    \     |     /
    ##     \____d____/
    ##
    ##
    mat <- matrix(c(0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0),
                  byrow=TRUE, ncol=4)
    rownames(mat) <- letters[1:4]
    colnames(mat) <- letters[1:4]
    mat
}


simpleDirectedGraph <- function() {
    ## Here's a simple graph for testing
    ##    a           b
    ##    |\         /^
    ##    | \__>c<__/ |
    ##    |     ^     |
    ##    \     |     /
    ##     \___>d____/
    ##
    ##
    mat <- matrix(c(0, 0, 1, 1,
                    0, 0, 1, 0,
                    0, 0, 0, 0,
                    0, 1, 1, 0),
                  byrow=TRUE, ncol=4)
    rownames(mat) <- letters[1:4]
    colnames(mat) <- letters[1:4]
    mat
    graphAM(adjMat=mat, edgemode="directed")
}


testNodeDataDefaults <- function() {
    mat <- simpleInciMat()
    g1 <- graphAM(adjMat=mat)

    ## If no attributes have been defined, empty list.
    checkEquals(list(), nodeDataDefaults(g1))

    ## Can assign a named list
    myEdgeAttributes <- list(foo=1, bar="blue")
    nodeDataDefaults(g1) <- myEdgeAttributes
    checkEquals(myEdgeAttributes, nodeDataDefaults(g1))

    checkEquals(myEdgeAttributes$foo, nodeDataDefaults(g1, attr="foo"))

    nodeDataDefaults(g1, attr="size") <- 400
    checkEquals(400, nodeDataDefaults(g1, attr="size"))

    checkException(nodeDataDefaults(g1, attr="NOSUCHATTRIBUTE"), silent=TRUE)
    checkException(nodeDataDefaults(g1) <- list(1, 3, 4),
                   silent=TRUE) ## must have names
}


testEdgeDataDefaults <- function() {
    mat <- simpleInciMat()
    g1 <- graphAM(adjMat=mat)

    ## If no attributes have been defined, empty list.
    checkEquals(list(), edgeDataDefaults(g1))

    ## Can assign a named list
    myEdgeAttributes <- list(foo=1, bar="blue")
    edgeDataDefaults(g1) <- myEdgeAttributes
    checkEquals(myEdgeAttributes, edgeDataDefaults(g1))

    checkEquals(myEdgeAttributes$foo, edgeDataDefaults(g1, attr="foo"))

    edgeDataDefaults(g1, attr="size") <- 400
    checkEquals(400, edgeDataDefaults(g1, attr="size"))

    checkException(edgeDataDefaults(g1, attr="NOSUCHATTRIBUTE"), silent=TRUE)
    checkException(edgeDataDefaults(g1) <- list(1, 3, 4),
                   silent=TRUE) ## must have names
}


testNodeDataGetting <- function() {
    mat <- simpleInciMat()
    g1 <- graphAM(adjMat=mat)
    myAttributes <- list(size=1, dim=c(3, 3), name="fred")
    nodeDataDefaults(g1) <- myAttributes

    checkEquals("fred", nodeData(g1, "a", attr="name")[[1]])

    someNodes <- c("a", "b")
    expect <- as.list(c(1, 1))
    names(expect) <- someNodes
    checkEquals(expect, nodeData(g1, n=someNodes, attr="size"))

    expect <- as.list(rep("fred", length(nodes(g1))))
    names(expect) <- nodes(g1)
    checkEquals(expect, nodeData(g1, attr="name"))

    checkEquals(myAttributes, nodeData(g1, n="a")[[1]])

    everything <- nodeData(g1)
    for (alist in everything)
      checkEquals(myAttributes, alist)
}


testNodeDataSetting <- function() {
    mat <- simpleInciMat()
    g1 <- graphAM(adjMat=mat)
    myAttributes <- list(size=1, dim=c(3, 3), name="fred")
    nodeDataDefaults(g1) <- myAttributes

    ## unknown node is error
    checkException(nodeData(g1, n="UNKNOWN_NODE", attr="size") <- 5, silent=TRUE)

    ## unknown attr is error
    checkException(nodeData(g1, n="a", attr="UNKNOWN") <- 5, silent=TRUE)

    nodeData(g1, n="a", attr="size") <- 5
    checkEquals(5, nodeData(g1, n="a", attr="size")[[1]])

    nodeData(g1, n=c("a", "b", "c"), attr="size") <- 50
    expect <- myAttributes
    expect[["size"]] <- 50
    checkEquals(list(a=expect, b=expect, c=expect),
                nodeData(g1, n=c("a", "b", "c")))

    nodeData(g1, n=c("a", "b", "c"), attr="size") <- c(1, 2, 3)
    checkEquals(c(1, 2, 3),
                as.numeric(nodeData(g1, n=c("a", "b", "c"), attr="size")))

    nodeData(g1, attr="name") <- "unknown"
    expect <- as.list(rep("unknown", length(nodes(g1))))
    names(expect) <- nodes(g1)
    checkEquals(expect, nodeData(g1, attr="name"))                  
}


testEdgeDataGetting <- function() {
    mat <- simpleInciMat()
    g1 <- graphAM(adjMat=mat)
    myAttributes <- list(size=1, dim=c(3, 3), name="fred")
    edgeDataDefaults(g1) <- myAttributes

    checkEquals("fred", edgeData(g1, from="a", to="d", attr="name")[[1]])

    fr <- c("a", "b")
    to <- c("c", "c")
    expect <- as.list(c(1, 1))
    names(expect) <- c("a|c", "b|c")
    checkEquals(expect, edgeData(g1, fr, to, attr="size"))

    expect <- rep("fred", sum(sapply(edges(g1), length)))
    checkEquals(expect, as.character(edgeData(g1, attr="name")))

    checkEquals(myAttributes, edgeData(g1, from="a", to="c")[[1]])

    everything <- edgeData(g1)
    for (alist in everything)
      checkEquals(myAttributes, alist)

    got <- edgeData(g1, from="d", attr="size")
    checkEquals(3, length(got))
    checkEquals(rep(1, 3), as.numeric(got))

    got <- edgeData(g1, to="d", attr="size")
    checkEquals(3, length(got))
    checkEquals(rep(1, 3), as.numeric(got))

    expect <- c("a|c", "a|d", "d|a", "d|b", "d|c")
    checkEquals(expect, names(edgeData(g1, from=c("a", "d"), attr="name")))
}


testEdgeDataToOnlyUndir <- function() {
    mat <- simpleInciMat()
    mat[1, 3] <- mat[3, 1] <- 100
    mat[1, 4] <- mat[4, 1] <- 200
    g1 <- graphAM(adjMat=mat, values=list(weight=1))
    got <- edgeData(g1, to=c("a", "b"), attr="weight")
    expect <- c("c|a", "d|a", "c|b", "d|b")
    checkEquals(expect, names(got))
}


testEdgeDataToOnlyDir <- function() {
    g1 <- simpleDirectedGraph()
    edgeDataDefaults(g1, attr="weight") <- 1
    edgeData(g1, from=c("a", "b"), to=c("c", "c"), attr="weight") <- c(10, 20)
    got <- edgeData(g1, to=c("a", "b"), attr="weight")
    expect <- c("d|b")
    checkEquals(expect, names(got))
}


testEdgeDataSettingDirected <- function() {
    g1 <- simpleDirectedGraph()
    myAttributes <- list(size=1, dim=c(3, 3), name="fred")
    edgeDataDefaults(g1) <- myAttributes

    edgeData(g1, from="a", to="d", attr="name") <- "Joe"
    expect <- myAttributes
    expect[["name"]] <- "Joe"
    checkEquals(expect, edgeData(g1, from="a", to="d")[[1]])

    fr <- c("a", "b")
    to <- c("c", "c")
    expect <- as.list(c(5, 5))
    names(expect) <- c("a|c", "b|c")
    edgeData(g1, fr, to, attr="size") <- 5
    checkEquals(expect, edgeData(g1, fr, to, attr="size"))

    expect <- as.list(c(10, 20))
    names(expect) <- c("a|c", "b|c")
    edgeData(g1, fr, to, attr="size") <- c(10, 20)
    checkEquals(expect, edgeData(g1, fr, to, attr="size"))

    edgeData(g1, from="a", attr="size") <- 555
    checkEquals(rep(555, 2), as.numeric(edgeData(g1, from="a", attr="size")))

    edgeData(g1, to="b", attr="size") <- 111
    checkEquals(111, as.numeric(edgeData(g1, to="b", attr="size")))
}


testEdgeDataSettingUndirected <- function() {
    mat <- simpleInciMat()
    g1 <- graphAM(adjMat=mat)
    myAttributes <- list(size=1, dim=c(3, 3), name="fred")
    edgeDataDefaults(g1) <- myAttributes

    edgeData(g1, from="a", to="d", attr="name") <- "Joe"
    expect <- myAttributes
    expect[["name"]] <- "Joe"
    checkEquals(expect, edgeData(g1, from="a", to="d")[[1]])
    ## verify reciprocal edge data was set
    checkEquals("Joe", edgeData(g1, from="d", to="a", attr="name")[[1]])

    fr <- c("a", "b")
    to <- c("c", "c")
    expect <- as.list(c(5, 5))
    names(expect) <- c("a|c", "b|c")
    edgeData(g1, fr, to, attr="size") <- 5
    checkEquals(expect, edgeData(g1, fr, to, attr="size"))
    names(expect) <- c("c|a", "c|b")
    checkEquals(expect, edgeData(g1, to, fr, attr="size"))

    expect <- as.list(c(10, 20))
    names(expect) <- c("a|c", "b|c")
    edgeData(g1, fr, to, attr="size") <- c(10, 20)
    checkEquals(expect, edgeData(g1, fr, to, attr="size"))
    names(expect) <- c("c|a", "c|b")
    checkEquals(expect, edgeData(g1, to, fr, attr="size"))
    
    edgeData(g1, from="a", attr="size") <- 555
    checkEquals(rep(555, 2), as.numeric(edgeData(g1, from="a", attr="size")))
    checkEquals(555, edgeData(g1, from="c", to="a", attr="size")[[1]])

    edgeData(g1, to="b", attr="size") <- 111
    checkEquals(rep(111, 2), as.numeric(edgeData(g1, to="b", attr="size")))
    checkEquals(111, edgeData(g1, from="c", to="b", attr="size")[[1]])
}


testEdgeDataSettingFromOnly <- function() {
    mat <- simpleInciMat()
    g1 <- graphAM(adjMat=mat)
    myAttributes <- list(size=1, dim=c(3, 3), name="fred")
    edgeDataDefaults(g1) <- myAttributes

    expect <- rep("fred", 5)
    got <- unlist(edgeData(g1, from=c("a", "d"), attr="name"), use.names=FALSE)
    checkEquals(expect, got, "precondition check")
    
    edgeData(g1, from=c("a", "d"), attr="name") <- "Sam"
    expect <- rep("Sam", 5)
    got <- unlist(edgeData(g1, from=c("a", "d"), attr="name"), use.names=FALSE)
    checkEquals(expect, got, "use from only in assign")
}


testNormalizeEdges <- function() {
    checkException(graph:::.normalizeEdges(c("b", "d"), c("a", "b", "c")), silent=TRUE)
    checkException(graph:::.normalizeEdges(c("a", "b", "c"), c("a", "e")), silent=TRUE)

    f <- letters[1:10]
    t <- letters[11:20]
    checkEquals(list(from=f, to=t), graph:::.normalizeEdges(f, t))

    checkEquals(list(from=c("a", "a", "a"), to=c("a", "b", "c")),
                graph:::.normalizeEdges("a", c("a", "b", "c")))

    checkEquals(list(from=c("a", "b", "c"), to=c("d", "d", "d")),
                graph:::.normalizeEdges(c("a", "b", "c"), "d"))
}

