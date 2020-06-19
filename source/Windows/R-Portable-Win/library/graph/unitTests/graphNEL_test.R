

##.setUp <- function() RNGkind("default", "default")

simpleGraphNEL <- function() {
     V <- letters[1:4]
     edL <- vector("list", length=4)
     names(edL) <- V
     edL[["a"]] <- list(edges=c(3, 4), weights=c(.13, .14))
     edL[["b"]] <- list(edges=c(3, 4), weights=c(.23, .24))
     edL[["c"]] <- list(edges=c(1, 2, 4), weights=c(.13, .23, .34))
     edL[["d"]] <- list(edges=c(1, 2, 3), weights=c(.14, .24, .34))
     gR <- graphNEL(nodes=V, edgeL=edL)
     gR
 }


simpleDirectedGraphNEL <- function() {
    set.seed(123)
     V <- letters[1:4]
     edL <- vector("list", length=4)
     names(edL) <- V
     edL[["a"]] <- list(edges=c(3, 4), weights=c(.13, .14))
     edL[["b"]] <- list(edges=c(3), weights=.23)
     edL[["c"]] <- list(edges=numeric(0), weights=numeric(0))
     edL[["d"]] <- list(edges=c(2, 3), weights=c(.42, .43))
     gR <- graphNEL(nodes=V, edgeL=edL, edgemode="directed")
     gR
}

testConstructorFunction <- function() {
    nodes <- LETTERS[1:4]
    edgeL <- list(A=c("B", "C"), B="C", C="D")

    ## no-argument constructor
    target <- new("graphNEL")
    checkIdentical(target, graphNEL())

    ## node / edgeList constructor
    target <- new("graphNEL", nodes=nodes, edgeL=edgeL, edgemode="directed")
    checkIdentical(target, graphNEL(nodes, edgeL, "directed"))

    ## edgemode default == "undirected"
    edgeL2 <- list(A = c("B", "C"), B = c("A", "C"),
                   C = c("A", "B", "D"), D = "C")
    target <- new("graphNEL", nodes=nodes, edgeL=edgeL2)
    checkIdentical(target, graphNEL(nodes, edgeL2, "undirected"))
}

testCreateBadNodeNames <- function() {
    badNodeName <- paste("foo", graph:::EDGE_KEY_SEP, "bar", sep="")
    checkException(graphNEL(nodes=badNodeName), silent=TRUE)
    checkException(graphNEL(nodes=c(NA, "b")), silent=TRUE)
    checkException(graphNEL(nodes=c("a", "")), silent=TRUE)
}


testIsAdjacent <- function() {
    g1 <- simpleGraphNEL()

    checkEquals(FALSE, isAdjacent(g1, "a", "b"))
    checkEquals(TRUE, isAdjacent(g1, "a", "c"))

    expect <- c(FALSE, TRUE, TRUE)
    got <- isAdjacent(g1, c("a", "a", "a"), c("b", "c", "d"))
    checkEquals(expect, got)
}


testNumEdges <- function() {
    mat <- matrix(c(1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0), ncol=4)
    rownames(mat) <- letters[1:4]
    colnames(mat) <- letters[1:4]
    g <- as(mat, "graphNEL")
    checkEquals(4, numEdges(g))
}


testInEdges <- function() {
    g <- simpleDirectedGraphNEL()
    expectedInEdges <- list(a=character(0), b="d", c=c("a", "b", "d"),
                            d="a")
    checkEquals(expectedInEdges, inEdges(g))
    checkEquals(expectedInEdges, inEdges(object=g))
    n <- c("a", "d")
    checkEquals(expectedInEdges[n], inEdges(n, g))

    ## verify unknown node behavior
    ans <- tryCatch(inEdges("not-a-node", g),
                    error = function(e) e)
    checkEquals("not a node: 'not-a-node'", conditionMessage(ans))
}

testEmptyGraph <- function() {
    g <- graphNEL()
    checkEquals(0, numEdges(g))
    checkEquals(0, numNodes(g))
}


testCreateGraphNoEdges <- function() {
    g <- graphNEL(nodes=c("a", "b"))
    checkEquals(0, numEdges(g))
    checkEquals(2, numNodes(g))

    g <- graphNEL(nodes=c("a", "b"), edgeL=list())
    checkEquals(0, numEdges(g))
    checkEquals(2, numNodes(g))

    checkEquals(2, length(edges(g)))
    checkEquals(nodes(g), names(edges(g)))
    checkEquals(0, sum(sapply(edges(g), length)))
}


testConstructor <- function() {
    g <- simpleGraphNEL()
    g2 <- graphNEL(nodes=nodes(g), edgeL=edges(g))
    checkEquals(nodes(g), nodes(g2))
    checkEquals(edges(g), edges(g2))

    ## We also support the more complicated list structure for describing graph
    ## edges.
    g2 <- graphNEL(nodes=nodes(g), edgeL=g@edgeL)
    checkEquals(nodes(g), nodes(g2))
    checkEquals(edges(g), edges(g2))
}


testNullHandlingInEdgeL <- function() {
    g <- simpleDirectedGraphNEL()
    eL <- g@edgeL
    eL <- c(eL[c("a", "b", "c")], list(d=NULL))
    g2 <- graphNEL(nodes(g), eL, "directed")
    checkTrue(all(sapply(g2@edgeL, function(x) !is.null(x))))
}


testCaptureWeightsWithEdgeLUndirected <- function() {
    g <- simpleGraphNEL()
    expect <- as.list(c(.13, .14))
    names(expect) <- c("a|c", "a|d")
    checkEquals(expect, edgeData(g, from="a", attr="weight"))
}


testCaptureWeightsWithEdgeLDirected <- function() {
    g <- simpleDirectedGraphNEL()
    expect <- as.list(c(.13, .14))
    names(expect) <- c("a|c", "a|d")
    checkEquals(expect, edgeData(g, from="a", attr="weight"))
}


testAddNode <- function() {
    g1 <- simpleGraphNEL()
    newNodes <- c("r", "s", "a", "b")
    checkException(addNode(newNodes, g1), silent=TRUE)

    newNodes <- c("r", "s")
    expect <- c(nodes(g1), newNodes)
    g1 <- addNode(newNodes, g1)
    checkEquals(expect, nodes(g1))
}


testAddNodeWithEdges <- function() {
    g1 <- simpleGraphNEL()
    newNodes <- c("r", "s", "t")
    newEdges <- list(r=c("a", "s"), s="b", t=character(0))
    g2 <- addNode(newNodes, g1, newEdges)

    checkEquals(c(nodes(g1), newNodes), nodes(g2))
    expect <- list(r=c("a", "s"))
    checkEquals(expect, edges(g2)["r"])
    expectEdges <- edges(g1)
    expectEdges[["a"]] <- c(expectEdges[["a"]], "r")
    expectEdges[["b"]] <- c(expectEdges[["b"]], "s")
    expectEdges[["r"]] <- c("a", "s")
    expectEdges[["s"]] <- c("r", "b")
    expectEdges[["t"]] <- character(0)
    checkEquals(expectEdges, edges(g2))
}


testAddNodeWithEdgesAndWeights <- function() {
    g1 <- simpleGraphNEL()
    newNodes <- c("r", "s", "t")
    newEdges <- list(r=c(a=11, s=22), s=c(b=33), t=numeric(0))
    g2 <- addNode(newNodes, g1, newEdges)

    checkEquals(c(nodes(g1), newNodes), nodes(g2))
    expect <- list(r=c("a", "s"))
    checkEquals(expect, edges(g2)["r"])
    expectEdges <- edges(g1)
    expectEdges[["a"]] <- c(expectEdges[["a"]], "r")
    expectEdges[["b"]] <- c(expectEdges[["b"]], "s")
    expectEdges[["r"]] <- c("a", "s")
    expectEdges[["s"]] <- c("r", "b")
    expectEdges[["t"]] <- character(0)
    checkEquals(expectEdges, edges(g2))
}


testAddNodeBadNodeName <- function() {
    g1 <- simpleGraphNEL()
    badNodeName <- paste("foo", graph:::EDGE_KEY_SEP, "bar", sep="")
    checkException(addNode(badNodeName, g1), silent=TRUE)
}


testSubGraphNoEdges <- function() {
    g1 <- simpleGraphNEL()
    g1 <- removeEdge("a", c("c", "d"), g1)
    g2 <- subGraph("a", g1) ## g2 has no edges
    checkEquals(0, numEdges(g2))
    checkEquals(1, numNodes(g2))
}


testSubGraphNoEdgesDirected <- function() {
    g1 <- simpleDirectedGraphNEL()
    g1 <- removeEdge("a", c("c", "d"), g1)
    g2 <- subGraph("a", g1) ## g2 has no edges
    checkEquals(0, numEdges(g2))
    checkEquals(1, numNodes(g2))
 }


testSubGraphAttributes <- function() {
    g1 <- simpleDirectedGraphNEL()
    nodeDataDefaults(g1) <- list(w=NA, n="")
    nodeData(g1, n=c("a", "b"), attr="w") <- c(1, 2)
    nodeData(g1, n=c("a", "b"), attr="n") <- c("A", "B")
    edgeDataDefaults(g1) <- list(x=NA)
    edgeData(g1, from="a", to="d", attr="x") <- 6
    edgeData(g1, from="a", to="c", attr="x") <- 7

    g2 <- subGraph(c("a", "d"), g1)
    checkEquals("a", names(g2@nodeData))

    g3 <- subGraph(c("a", "b", "c"), g1)
    checkEquals(c("a|c", "b|c"), names(g3@edgeData))
}


testRemoveEdgeUndirected <- function() {
    g <- simpleGraphNEL()
    g1 <- removeEdge("a", c("c", "d"), g)
    checkEquals(3, numEdges(g1))
    eD <- edges(g1)
    checkEquals(character(0), eD$a)
    checkEquals(c("c", "d"), eD$b)

    g2 <- removeEdge(c("c", "d"), "a", g)
    checkEquals(3, numEdges(g2))
    eD <- edges(g2)
    checkEquals(character(0), eD$a)
    checkEquals(c("c", "d"), eD$b)
}


testRemoveEdgeDirected <- function() {
    g1 <- simpleDirectedGraphNEL()
    f <- c("a", "a")
    t <- c("c", "d")
    g2 <- removeEdge(from=f, to=t, g1)
    checkEquals(3, numEdges(g2))
    checkTrue(!length(edges(g2)[["a"]]))
}


testRemoveEdgeLarge <- function() {
    ## This test is from Denise Scholtens
    set.seed(678)
    N <- 500
    numEdges <- 2500
    nodes <- paste("n", 1:500, sep="")
    g <- randomEGraph(nodes, edges=numEdges)
    edgemode(g) <- "directed"
    checkEquals(numEdges*2, numEdges(g))
     from <- c("n1","n2","n2","n3","n5","n7","n7","n8","n8","n8","n9","n9",
               "n9")
    to <- c("n255","n383","n261","n381","n234","n225","n315","n38","n296",
            "n78","n310","n19","n422")

    g1 <- removeEdge(from, to, g)
    checkEquals(numEdges*2 - length(from), numEdges(g1))
}

testRemoveEdgeLarge2 <- function() {
    ## This test is from a bug discovered by Dan Bebber
    From <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
              2, 5, 5, 8, 8, 8, 8, 8, 11, 11, 12, 12, 12, 14, 14, 16, 16, 20,
              20, 23, 23, 24, 24, 25, 25, 29, 29, 32, 32, 32, 32, 38, 38, 41,
              41, 41, 43, 43, 54, 54, 59, 59, 60, 68, 68, 69, 72, 82, 83, 88,
              88, 88, 88, 89, 89, 90, 90, 96, 97, 98, 98, 98)

    To <- c(2, 5, 8, 11, 12, 14, 16, 20, 23, 37, 38, 54, 57, 68, 72, 81,
            86, 87, 97, 88, 100, 32, 102, 38, 41, 49, 51, 53, 58, 59, 60,
            63, 67, 71, 72, 75, 76, 84, 85, 24, 29, 25, 28, 26, 27, 30, 31,
            33, 34, 35, 36, 40, 41, 42, 43, 47, 44, 45, 55, 56, 60, 62, 61,
            70, 71, 71, 73, 84, 84, 89, 90, 93, 94, 97, 141, 95, 141, 141,
            98, 99, 141, 158)

    FT <- matrix(c(From, To), ncol=2) #create a 'from-to' matrix

    g <- ftM2graphNEL(FT, edgemode="undirected")

    gr <- removeEdge(from=as.character(From[1:2]),
                     to=as.character(To[1:2]), g)

    checkEquals(numEdges(g) - 2, numEdges(gr))

    gr <- removeEdge(from=as.character(From[1:20]),
                     to=as.character(To[1:20]), g)
    checkEquals(numEdges(g) - 20, numEdges(gr))
}



test_eWV <- function() {
    V <- LETTERS[1:4]
    gR <- graphNEL(nodes=V)
    gX <- addEdge("A", "C", gR, 0.2)

    ans <- eWV(gX, edgeMatrix(gX), useNNames = TRUE)
    checkEquals(c("A--C"=0.2), ans)
}


testEdgeWeightsNoEdges <- function() {
    g <- graphNEL(nodes=letters[1:6])
    expect <- lapply(edges(g), as.numeric)
    checkEquals(expect, edgeWeights(g))
}


testRemoveNode1 <- function() {
    ## using the example from the removeNode help page
    V <- LETTERS[1:4]
    edL2 <- vector("list", length=4)
    names(edL2) <- V
    for(i in 1:4)
      edL2[[i]] <- list(edges=c(2,1,2,1)[i],
                        weights=sqrt(i))
    gR2 <- graphNEL(nodes=V, edgeL=edL2, edgemode="directed")

    gX <- removeNode("C", gR2)
    checkEquals(c("A", "B", "D"), nodes(gX))

    gY <- removeNode(c("A","D"), gX)
    checkEquals("B", nodes(gY))

    gZ <- removeNode(c("A","C","D"), gR2)
    checkEquals("B", nodes(gZ))

    ## XXX: using direct slot access to verify that edge attributes
    ##      have been completely removed.
    checkTrue(length(gZ@edgeData@data) == 0)
}


testRemoveNode2 <- function() {
    g <- simpleDirectedGraphNEL()
    nds <- nodes(g)
    for (n in nds) {
        g2 <- removeNode(n, g)
        checkEquals(nds[nds != n], nodes(g2))
    }
}


test_ugraph <- function() {
    g <- simpleDirectedGraphNEL()
    ug <- ugraph(g)
    eg <- simpleGraphNEL()
    checkTrue(isDirected(g))
    checkTrue(!isDirected(ug))
    checkEquals(nodes(g), nodes(ug))
    checkEquals(nodes(eg), nodes(ug))
    ## verify edges
    eGot <- edges(ug)[nodes(g)]
    eExp <- edges(eg)[nodes(g)]
    for (n in nodes(g)) {
        checkTrue(setequal(eExp[[n]], eGot[[n]]))
    }
}


test_rename_nodes_edgeWeights <- function() {
    g <- simpleGraphNEL()
    ew <- edgeWeights(g)
    ew <- lapply(ew, function(x) {
        names(x) <- toupper(names(x))
        x
    })
    names(ew) <- toupper(names(ew))
    nodes(g) <- LETTERS[1:4]
    checkEquals(LETTERS[1:4], nodes(g))
    checkEquals(ew, edgeWeights(g))
}


test_rename_nodes_nodeData <- function() {
    g <- simpleGraphNEL()
    nodeDataDefaults(g) <- list(type=NA)
    nodeData(g, n="a", attr="type") <- "the first one"
    nodeData(g, n="d", attr="type") <- "the last one"
    ndDef <- nodeDataDefaults(g)
    nd <- nodeData(g, attr="type")
    names(nd) <- toupper(names(nd))

    nodes(g) <- toupper(nodes(g))

    checkEquals(nd, nodeData(g, attr="type"))
}

test_subgraph_attrs <- function() {
    x <- graphNEL(nodes=c("a", "b"),
             edgeL=list(a="b", b="b"),
             edgemode="directed")
    defs <- list(tag="NONE")
    nodeDataDefaults(x) <- defs
    edgeDataDefaults(x) <- defs
    nodeData(x, n="a", attr="tag") <- "zoo"
    edgeData(x, "a", "b", attr="tag") <- "yes"

    gg <- subGraph(c("a", "b"), x)
    checkEquals(defs, nodeDataDefaults(gg))
    checkEquals(defs, edgeDataDefaults(gg))
    checkEquals("zoo", nodeData(gg, "a", attr="tag")[[1]])
    checkEquals("yes", edgeData(gg, "a", "b", attr="tag")[[1]])
}

test_ftM2_with_self_edges <- function() {

    ft <- cbind(c(1:5,1,5),c(1:5,3,2))
    W <- c(1:5,7,9)
    ## this failed till 2008-06-26:
    gr <- ftM2graphNEL(ft, W, edgemode="undirected")
    m <- as(gr, "matrix")
    g2 <- as(m, "graphNEL")
    m2 <- as(g2, "matrix")
    checkEquals(m2, m)
    checkEquals(which(m2 != 0), c(1,3,7,10,11,13,19,22,25))
}

test_coerce_matrix_round_trip <- function()
{
    V <- LETTERS[1:4]
    g <- graphNEL(nodes=V, edgemode="directed")
    g <- addEdge(V[1+0],V[1+1],g, 3)
    g <- addEdge(V[1+0],V[2+1],g, 1.5)
    g <- addEdge(V[1+0],V[3+1],g, 1.8)
    g <- addEdge(V[1+1],V[2+1],g, 4.3)
    g <- addEdge(V[1+2],V[3+1],g, 2.2)

    mat0 <- matrix(c(0, 0, 0, 0, 3, 0, 0, 0, 1.5, 4.3, 0, 0, 1.8, 0, 2.2, 0),
                   ncol=4, dimnames = list(LETTERS[1:4], LETTERS[1:4]))

    checkEquals(mat0, as(g, "matrix"))
    checkEquals(mat0, as(as(mat0, "graphNEL"), "matrix"))
}
