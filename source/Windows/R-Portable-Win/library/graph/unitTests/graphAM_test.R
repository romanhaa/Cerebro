simpleAdjMat <- function() {
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
    graphAM (adjMat=mat, edgemode="directed")
}

testConstructorFunction <- function() {
    ## no-argument constructor
    adjMat <- matrix(integer(), 0, 0)
    target <- new("graphAM", adjMat=adjMat)
    checkIdentical(target, graphAM())

    ## adjMat constructor
    adjMat <- simpleAdjMat()
    target <- new("graphAM", adjMat=adjMat, edgemode="directed")
    checkIdentical(target, graphAM(adjMat, "directed"))
    target <- new("graphAM", adjMat=adjMat, edgemode="undirected")
    checkIdentical(target, graphAM(adjMat, "undirected"))
    checkIdentical(target, graphAM(adjMat))

    ## values:  adjacency matrix non-zero (and necessarily positive) values
    ## will be used as edge weights, an edge attribute. the -1 is ignored
    ## in all cases, but retrievable via edgeDataDefaults
    ## values is to contain exactly one name, and one value
    ## the name indicates the edge attribute, whose values come from the
    ## elements of the adjacency matrix:  often, but not necessarily, set to 1

    values <- list(weight=-1)  
    target <- new("graphAM", adjMat=adjMat, edgemode="directed",
                  values=values)
    checkEquals (edgeData(target, 'a', 'c', attr='weight')[[1]], 1)
    checkEquals (edgeDataDefaults(target, 'weight'), -1)
    checkIdentical(target, graphAM(adjMat, "directed", values))
}

testInvalidNonSquare <- function() {
    mat <- cbind(c(0, 0, 1), c(1, 1, 1))
    checkException(graphAM (adjMat=mat), silent=TRUE)
}


testInvalidNegativeValues <- function() {
    mat <- matrix(c(0, 1, -4, -1), ncol=2)
    checkException(graphAM (adjMat=mat), silent=TRUE)
}


testInvalidNonSymmetric <- function() {
    mat <- matrix(c(0, 1, 1,
                    0, 0, 1,
                    0, 0, 0), ncol=3, byrow=TRUE)
    colnames(mat) <- letters[1:3]
    checkException(graphAM (adjMat=mat), silent=TRUE)
    checkException(graphAM (adjMat=mat, edgemode="undirected"), silent=TRUE)
    g1 <- graphAM (adjMat=mat, edgemode="directed")
}


testInvalidBadNodeNames <- function() {
    mat <- simpleAdjMat()
    n <- paste(letters[1:4], 1:4, sep=graph:::EDGE_KEY_SEP)
    colnames(mat) <- rownames(mat) <- n
    checkException(graphAM (adjMat=mat), silent=TRUE)

    colnames(mat) <- rownames(mat) <- c("a", "b", NA, "c")
    checkException(graphAM (adjMat=mat), silent=TRUE)

    colnames(mat) <- rownames(mat) <- c("a", "f", "", "d")
    checkException(graphAM (adjMat=mat), silent=TRUE)
}

test_empty_graph <- function() {
    mat <- matrix(integer(0), nrow=0, ncol=0)
    g <- graphAM (adjMat = mat)
    checkEquals(0L, numNodes(g))
    checkEquals(0L, numEdges(g))
    checkEquals(character(0), nodes(g))
    checkEquals(list(), edges(g))

    m <- as(g, "matrix")
    checkEquals(c(0L, 0L), dim(m))
    checkEquals(0L, length(m))

    g <- graphAM (adjMat = mat, values = list(weight = 1L))
    checkEquals(0L, numNodes(g))
    checkEquals(0L, numEdges(g))
    checkEquals(character(0), nodes(g))
    checkEquals(list(), edges(g))

    m <- as(g, "matrix")
    checkEquals(c(0L, 0L), dim(m))
    checkEquals(0L, length(m))
}

test_no_edge_graph <- function() {
    mat <- matrix(0L, nrow=3, ncol=3,
                  dimnames=list(letters[1:3], letters[1:3]))
    g <- graphAM (adjMat = mat)
    checkEquals(letters[1:3], nodes(g))
    checkEquals(0L, numEdges(g))
    want <- list(a = character(0), b = character(0), c = character(0))
    checkEquals(want, edges(g))
    m <- as(g, "matrix")
    checkEquals(c(3L, 3L), dim(m))
    checkTrue(all(m == 0L))

    g <- graphAM (adjMat = mat, values = list(weight = 1L))
    checkEquals(letters[1:3], nodes(g))
    checkEquals(0L, numEdges(g))
    checkEquals(want, edges(g))

    m <- as(g, "matrix")
    checkEquals(c(3L, 3L), dim(m))
    checkTrue(all(m == 0L))
}

testValuesToAttr <- function() {
    mat <- matrix(c(0, 0, 1, 2,
                    0, 0, 3, 0,
                    0, 0, 0, 0,
                    0, 4, 5, 0),
                  byrow=TRUE, ncol=4)
    rownames(mat) <- letters[1:4]
    colnames(mat) <- letters[1:4]
    g1 <- graphAM (adjMat=mat, edgemode="directed",
              values=list(weight=1))
    checkEquals(4, edgeData(g1, "d", "b", attr="weight")[[1]])
    checkEquals(3, edgeData(g1, "b", "c", attr="weight")[[1]])
    checkEquals(2, edgeData(g1, "a", "d", attr="weight")[[1]])
    checkEquals(1, edgeData(g1, "a", "c", attr="weight")[[1]])

    checkException(graphAM (adjMat=mat, edgemode="directed",
                         values=list(weight=1, not=2)), silent=TRUE)
    checkException(graphAM (adjMat=mat, edgemode="directed",
                         values=list("must", "name")), silent=TRUE)
    checkException(graphAM (adjMat=mat, edgemode="directed",
                         values="weight"), silent=TRUE)

    g1 <- graphAM (adjMat=mat, edgemode="directed",
              values=list(type=4))
    checkEquals(4, edgeData(g1, "d", "b", attr="type")[[1]])
    checkEquals(3, edgeData(g1, "b", "c", attr="type")[[1]])
    checkEquals(2, edgeData(g1, "a", "d", attr="type")[[1]])
    checkEquals(1, edgeData(g1, "a", "c", attr="type")[[1]])
}


testEdges <- function() {
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)
    got <- edges(g1)
    expect <- list(a=c("c", "d"), b=c("c", "d"), c=c("a", "b", "d"),
                   d=c("a", "b", "c"))
    checkEquals(expect, got)

    got <- edges(g1, c("a", "d"))
    expect <- expect[c("a", "d")]
    checkEquals(expect, got)
}


testEdgesDirected <- function() {
    g1 <- simpleDirectedGraph()
    expect <- list(a=c("c", "d"), b="c", c=character(0),
                   d=c("b", "c"))
    checkEquals(expect, edges(g1))
}


testEdgesSubset <- function() {
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)
    got <- edges(g1)
    expect <- list(a=c("c", "d"), d=c("a", "b", "c"))
    got <- edges(g1, c("a", "d"))
    checkEquals(expect, got)
}


testNodeNames <- function() {
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)
    got <- nodes(g1)
    expect <- letters[1:4]
    checkEquals(expect, got)
}


testNodeNamesReplace <- function() {
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)
    nodes(g1) <- LETTERS[1:4]
    expect <- LETTERS[1:4]
    checkEquals(expect, nodes(g1))
}


testNumNodes <- function() {
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)
    checkEquals(nrow(mat), numNodes(g1))
}


testNumEdges <- function() {
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)
    checkEquals(5, numEdges(g1))

    edgemode(g1) <- "directed"
    checkEquals(10, numEdges(g1))
}


testNumEdgesWithSelfLoop <- function() {
    mat <- matrix(c(1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0), ncol=4)
    g1 <- graphAM (adjMat=mat)
    checkEquals(4, numEdges(g1))
}

testIsAdjacent <- function() {
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)

    checkEquals(TRUE, isAdjacent(g1, "a", "c"))
    checkEquals(TRUE, isAdjacent(g1, "c", "a"))
    checkEquals(FALSE, isAdjacent(g1, "a", "b"))
    checkEquals(FALSE, isAdjacent(g1, "b", "a"))
    checkException(isAdjacent(g1, "z", "a"), silent=TRUE)
    checkException(isAdjacent(g1, "a", "z"), silent=TRUE)
}


testIsAdjacentVectorized <- function() {
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)

    fr <- c("a", "c", "a", "b")
    to <- c("c", "a", "b", "a")
    expect <- c(TRUE, TRUE, FALSE, FALSE)
    checkEquals(expect, isAdjacent(g1, fr, to))
    checkEquals(expect, isAdjacent(g1, to, fr))
}


## testSubgraph <- function() {
##     mat <- simpleAdjMat()
##     g1 <- graphAM (adjMat=mat)
##     g2 <- subgraph(c("a", "b", "c"), ffff)

##                }


testSimpleEdgeWeights <- function() {
    mat <- simpleAdjMat()
    g <- graphAM (mat)
    checkEquals(nodes(g), names(edgeWeights(g)))
    expect <- c(c=1:1, d=1:1)
    checkEquals(expect, edgeWeights(g)$a)
}


testAddNode <- function() {
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)

    newNodes <- c("r", "s", "a", "b")
    checkException(addNode(newNodes, g1), silent=TRUE)

    newNodes <- c("r", "s")
    expect <- c(nodes(g1), newNodes)
    g1 <- addNode(newNodes, g1)
    checkEquals(expect, nodes(g1))

    badNodeName <- paste("foo", graph:::EDGE_KEY_SEP, "bar", sep="")
    checkException(addNode(badNodeName, g1), silent=TRUE)
}


testAddEdge <- function() {
    ## I would like different order of the args in the generic, but not sure it is
    ## worth it to change... but would seem more consistent.
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)
    g1 <- addNode("e", g1)
    checkEquals(FALSE, isAdjacent(g1, "b", "e"))
    g1 <- addEdge(graph=g1, from="b", to="e")
    checkEquals(TRUE, isAdjacent(g1, "b", "e"))
}

testAddEdgeMultiple <- function()
{
    a <- matrix(0L, nrow=8, ncol=8)
    dimnames(a) <- list(letters[1:8], letters[1:8])
    G <- graphAM (adjMat=a, edgemode = "directed")
    GU <- graphAM (adjMat=a)
    ## make sure we don't warn for this call
    tryCatch({
        H <- addEdge(from=c("a", "b", "c"), to=c("d", "e", "f"), G)
        HU <- addEdge(from=c("a", "b", "c"), to=c("d", "e", "f"), GU)
    }, warning = function(w)
             stop("unwanted warning message: ", conditionMessage(w)))
    expect <- a
    fr <- c("a", "b", "c")
    to <- c("d", "e", "f")
    wh <- cbind(match(fr, letters[1:8]), match(to, letters[1:8]))
    expect[wh] <- 1L
    checkEquals(expect, as(H, "matrix"))
    expectU <- expect
    expectU[wh[ , c(2L, 1L)]] <- 1L
    checkEquals(expectU, as(HU, "matrix"))
}

testClearNode <- function() {
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)
    edgeDataDefaults(g1, attr="weight") <- 1
    edgeData(g1, "a", "c", attr="weight") <- 400

    checkEquals(TRUE, isAdjacent(g1, "a", "c"))
    checkEquals(TRUE, isAdjacent(g1, "a", "d"))
    checkEquals(400, edgeData(g1, "a", "c", attr="weight")[[1]])
    g1 <- clearNode("a", g1)
    checkEquals(FALSE, isAdjacent(g1, "a", "c"))
    checkEquals(FALSE, isAdjacent(g1, "a", "d"))
    checkException(edgeData(g1, "a", "c", attr="weight"), silent=TRUE)
}


testRemoveNode <- function() {
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)
    origNodes <- nodes(g1)

    checkEquals(TRUE, "b" %in% origNodes)
    g1 <- removeNode("b", g1)
    checkEquals(FALSE, "b" %in% nodes(g1))
}


testRemoveEdge <- function() {
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)

    checkEquals(TRUE, isAdjacent(g1, "a", "c"))
    g1 <- removeEdge("a", "c", g1)
    checkEquals(FALSE, isAdjacent(g1, "a", "c"))
}

testRemoveEdgeWithWeights <- function() {
    mat <- simpleAdjMat()
    mat[mat != 0] <- runif(length(mat[mat != 0]))
    g <- graphAM (adjMat = mat, edgemode = "directed",
             values = list(weight = 1.0))
    weights <- unlist(edgeData(g, attr = "weight"))
    toRemove <- names(weights[weights < 0.5])
    expect <- numEdges(g) - length(toRemove)
    fromTo <- do.call(rbind, strsplit(toRemove, "|", fixed = TRUE))
    g2 <- removeEdge(fromTo[, 1], fromTo[, 2], g)
    checkEquals(expect, numEdges(g2))
    apply(fromTo, 1, function(row) {
        checkEquals(FALSE, isAdjacent(g2, row[1], row[2]))
    })
}


testGraphAMCloning <- function() {
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)
    origNodes <- nodes(g1)

    g2 <- g1

    ## modify g1
    g1 <- addNode("NEW", g1)
    edgeDataDefaults(g1, "weight") <- 2
    edgeDataDefaults(g1, "color") <- "green"
    ## g2 should not have changed
    checkEquals(list(), edgeDataDefaults(g2))
    checkEquals(origNodes, nodes(g2))
}


testUndirectedAsGraphNEL <- function() {
    mat <- simpleAdjMat()
    g1 <- graphAM (adjMat=mat)
    gNel <- as(g1, "graphNEL")
    checkEquals(edges(g1), edges(gNel))
    checkEquals(nodes(g1), nodes(gNel))
    checkEquals(edgemode(g1), edgemode(gNel))
    checkEquals(edgeDataDefaults(g1), edgeDataDefaults(gNel))
    checkEquals(nodeDataDefaults(g1), nodeDataDefaults(gNel))
}


testDirectedAsGraphNEL <- function() {
    g1 <- simpleDirectedGraph()
    gNel <- as(g1, "graphNEL")
    checkEquals(edges(g1), edges(gNel))
    checkEquals(nodes(g1), nodes(gNel))
    checkEquals(edgemode(g1), edgemode(gNel))
    checkEquals(edgeDataDefaults(g1), edgeDataDefaults(gNel))
    checkEquals(nodeDataDefaults(g1), nodeDataDefaults(gNel))
}


testDirectedAsGraphAM <- function() {
    g1 <- simpleDirectedGraph()
    gNel <- as(g1, "graphNEL")
    g2 <- as(gNel, "graphAM")
    checkEquals(edges(g1), edges(g2))
    checkEquals(nodes(g1), nodes(g2))
    checkEquals(edgemode(g1), edgemode(g2))
    checkEquals(edgeDataDefaults(g1), edgeDataDefaults(g2))
    checkEquals(nodeDataDefaults(g1), nodeDataDefaults(g2))
}


testInEdges <- function() {
    g1 <- simpleDirectedGraph()
    expected <- list(a=character(0), b="d", c=c("a", "b", "d"), d="a")
    checkEquals(expected, inEdges(g1), msg="gramAM")
    checkEquals(expected, inEdges(object=g1), msg="gramAM")
    checkEquals(expected, inEdges(node=g1), msg="gramAM")
}


testNoEdges <- function() {
    m <- matrix(0, nrow=3, ncol=3)
    g <- graphAM (m)
    checkEquals(0, numEdges(g))
    checkEquals(3, length(edges(g)))
    checkEquals(nodes(g), names(edges(g)))
    checkEquals(0, sum(sapply(edges(g), length)))
}


testAsMatrix <- function() {
    mat <- rbind(c(0, 0, 12, 1),
                 c(0, 0, 1, 1),
                 c(12, 1, 0, 1),
                 c(1, 1, 1, 0))
    rownames(mat) <- colnames(mat) <- letters[1:4]
    ## If no values arg, then matrix just converted to 0/1
    g1 <- graphAM (adjMat=mat, edgemode="undirected")
    mat1 <- mat
    mat1[mat1 != 0] <- 1:1
    checkEquals(mat1, as(g1, "matrix"))

    ## With values arg, matrix values stored as edge attribute
    ## which gets restored for as(<.>, "matrix")
    g2 <- graphAM (adjMat=mat, edgemode="undirected",
              values=list(weight=1))
    checkEquals(mat, as(g2, "matrix"))
}

test_coerce_matrix_to_graphAM <- function()
{
    mat <- matrix(c(0, 0, 1, 2,
                    0, 0, 3, 0,
                    0, 0, 0, 0,
                    0, 4, 5, 0),
                  byrow=TRUE, ncol=4,
                  dimnames=list(letters[1:4], letters[1:4]))

    g <- as(mat, "graphAM")
    checkEquals(mat, as(g, "matrix"))

    g2 <- graphAM (adjMat=mat, edgemode="directed",
              values=list("weight"=1))
    checkEquals(as(g, "matrix"), as(g2, "matrix"))
}


test_edgeMatrix <- function() {
    ugam <- graphAM (adjMat=simpleAdjMat(), edgemode="undirected")
    gam <- simpleDirectedGraph()

    expect <- c("1+3", "1+4", "2+3", "2+4", "3+4")
    got <- edgeMatrix(ugam)
    checkTrue(setequal(expect, paste(got[1, ], got[2, ], sep="+")))
    checkEquals(list(c("from", "to"), NULL), dimnames(got))

    expect <- c("1+3", "1+4", "2+3", "4+2", "4+3")
    got <- edgeMatrix(gam)
    checkTrue(setequal(expect, paste(got[1, ], got[2, ], sep="+")))
    ## duplicates should have no effect on directed graph
    got <- edgeMatrix(gam, duplicates=TRUE)
    checkTrue(setequal(expect, paste(got[1, ], got[2, ], sep="+")))

    expect <- c("1+3", "1+4", "2+3", "2+4", "3+4",
                "3+1", "4+1", "3+2", "4+2", "4+3")
    got <- edgeMatrix(ugam, duplicates=TRUE)
    checkTrue(setequal(expect, paste(got[1, ], got[2, ], sep="+")))
}


test_rename_nodes_edgeWeights <- function() {
    mat <- matrix(c(0, 0, 1, 2,
                    0, 0, 3, 0,
                    0, 0, 0, 0,
                    0, 4, 5, 0),
                  byrow=TRUE, ncol=4)
    rownames(mat) <- letters[1:4]
    colnames(mat) <- letters[1:4]
    g <- graphAM (adjMat=mat, edgemode="directed",
              values=list(weight=1))
    ew <- edgeWeights(g)
    ew <- lapply(ew, function(x) {
        if (length(x))
          names(x) <- toupper(names(x))
        x
    })
    names(ew) <- toupper(names(ew))
    nodes(g) <- LETTERS[1:4]
    checkEquals(LETTERS[1:4], nodes(g))
    checkEquals(ew, edgeWeights(g))
}


test_rename_nodes_nodeData <- function() {
    g <- simpleDirectedGraph()
    nodeDataDefaults(g) <- list(type=NA)
    nodeData(g, n="a", attr="type") <- "the first one"
    nodeData(g, n="d", attr="type") <- "the last one"
    ndDef <- nodeDataDefaults(g)
    nd <- nodeData(g, attr="type")
    names(nd) <- toupper(names(nd))

    nodes(g) <- toupper(nodes(g))
    checkEquals(nd, nodeData(g, attr="type"))
}

