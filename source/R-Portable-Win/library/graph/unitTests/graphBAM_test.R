## library("graph")
set.seed(0x12a9b)
   library(graph)
   library(RUnit)

randBAMGraph <- function(numNodes = 10 , numEdges = 10)
{
    df <-  graph:::randFromTo(numNodes, numEdges)
    df$ft$weight = seq_len(numNodes)
    g <- graphBAM(df$ft, nodes = df$nodes, edgemode = "directed")
    g
}

make_smallBAM <- function() {
    from = c("a", "a", "a", "x", "x", "c")
    to   = c("b", "c", "x", "y", "c", "a")
    weight=c(3.4, 2.6, 1.7, 5.3, 1.6, 7.9)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "directed")
    g1
}

make_unDirectedBAM <- function() {

    from = c("a", "a", "a", "x", "x", "c")
    to   = c("b", "c", "x", "y", "c", "d")
    weight=c(3.4, 2.6, 1.7, 5.3, 1.6, 7.9)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "undirected")
    g1

}

create_bigBAM <- function()
{
    r1 <- randFromTo(100, 100)
    r1$ft$weight <- seq_len(100)
    g1 <- graphBAM(r1$ft, r1$nodes, edgemode="directed")
    g1
}

test_create_graphBAMSmall <- function() {

    from = c("a", "d", "d", "b")
    to = c("b", "a", "d", "c")
    weight= c(1.5, 3.1, 5.4, 1)
    nodes = c("a","b","c","d")
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, nodes, edgemode = "directed")
    g2 <- graphBAM(df, nodes, edgemode = "undirected")

    checkEquals(4L, numEdges(g1))
    checkEquals(isDirected(g1), TRUE)
    checkEquals(isAdjacent(g1, c("a", "d", "b"), c("b", "d", "c") ), c(TRUE,TRUE,TRUE))
    checkEquals(names(edges(g1)), c("a", "b", "c", "d"))
    k <- edges(g1)
    checkEquals(list(k$a, k$b, k$c, k$d), list("b", "c", character(0), c("a", "d")))
    w <- edgeWeights(g1)
    checkEquals(names(w), c("a", "b", "c", "d"))
    checkEquals(list(w$a, w$b, w$c, w$d), list(structure(1.5, names="b"),
            structure(1, names="c"), numeric(0), structure(c(3.1, 5.4),
            names= c("a", "d"))))
    checkEquals(4L, numNodes(g2))
    checkEquals(4L, numEdges(g2))
    checkEquals(isDirected(g2), FALSE)
    checkEquals(isAdjacent(g1, c("a","d","b"), c("b","d","c") ), c(TRUE,TRUE,TRUE))

}

test_BAMNodes <- function() {
    from = c("a", "a", "a", "x", "x", "c")
    to   = c("b", "c", "x", "y", "c", "a")
    weight=c(3.4, 2.6, 1.7, 5.3, 1.6, 7.9)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "directed")
    nds <- nodes(g1)
    checkIdentical(all(nds %in% unique(c(from,to))),TRUE)
    checkIdentical(isDirected(g1),TRUE)

    ## node names 
    from = paste0("X", 8:11)
    to   = paste0("X", 8:11) 
    df <- data.frame(from, to, weight=rep(1, 4))
    g2 <- graphBAM(df) ## no 'nodes'
    checkIdentical(nodes(g2), c("X10", "X11", "X8", "X9"))
    g2 <- graphBAM(df, nodes="X7") ## degree-zero node
    checkIdentical(nodes(g2), c("X10", "X11", "X7", "X8", "X9"))
    g2 <- graphBAM(df, nodes=paste0("X", 8:11)) ## forced ordering
    checkIdentical(nodes(g2), c("X8", "X9", "X10", "X11"))
}


checkBAMSubGraph <- function(g, subG) {
    nds <- nodes(g)
    subNodes <- nodes(subG)
    w1 <- g@edgeSet@weights
    ft1 <- .Call(graph:::graph_bitarray_rowColPos, g@edgeSet@bit_vector)
    origFromTo <- data.frame(from=nds[ft1[,"from"]], to = nds[ft1[,"to"]], weights = w1)

    w2 <- subG@edgeSet@weights
    ft2 <- .Call(graph:::graph_bitarray_rowColPos, subG@edgeSet@bit_vector)
    subFromTo <- data.frame(from = subNodes[ft2[,"from"]], to = subNodes[ft2[,"to"]], weights = w2)

    indx <- (origFromTo$from %in% subNodes) &
    (origFromTo$to %in% subNodes)
    want <- origFromTo[(origFromTo$from %in% subNodes) & (origFromTo$to %in% subNodes),]

    checkEquals(as.character(want$from), as.character(subFromTo$from))
    checkIdentical(as.character(want$to), as.character(subFromTo$to))
    checkEquals(g@edgeSet@weights[indx], subG@edgeSet@weights)
}

test_BAMSubGraph_Small <- function() {
    g1 <- make_smallBAM()
    sg <- subGraph(c("a","x", "y"), g1)
    checkIdentical(isDirected(sg), TRUE)
    checkIdentical(nodes(sg), c("a", "x", "y"))
    checkBAMSubGraph(g1,sg)
}


test_BAMSubGraph_Large  <- function() {
    g1 <- randBAMGraph(100,100)
    sn <- sample(nodes(g1), 55)
    sg <- subGraph( sn, g1)
    checkIdentical(isDirected(sg), TRUE)
    checkBAMSubGraph(g1,sg)
}


test_BAM_edgeWeights <- function() {
    g1 <- make_smallBAM()
    ew1 <- edgeWeights(g1)
    checkEquals(names(ew1), c("a", "b", "c", "x", "y"))
    checkEquals(list(ew1$a, ew1$b, ew1$c, ew1$x, ew1$y),
            list(structure( c(3.4, 2.6, 1.7), names = c("b","c","x")),
            numeric(0), structure(c(7.9), names = "a"),
            structure(c(1.6, 5.3), names= c("c", "y")), numeric(0)))

    ew2 <- edgeWeights(g1,c("a","b")) ##index = char
    checkEquals(names(ew2), c("a","b"))
    checkEquals(list(ew2$a, ew2$b), list(structure( c(3.4, 2.6, 1.7),
                            names = c("b","c","x")), numeric(0)))

    ew2 <- edgeWeights(g1, 1:2) ##index = numeric
    checkEquals(names(ew2), c("a","b"))
    checkEquals(list(ew2$a, ew2$b), list(structure( c(3.4, 2.6, 1.7),
                            names = c("b","c","x")), numeric(0)))
}

test_BAM_edgeWeights_undirected <- function()
{
      from = c("a", "d", "d", "b", "a")
        to = c("b", "a", "d", "c", "c")
    weight = c(1.5, 2.1, 3.4, 4.1, 5.6)
    df <- data.frame(from, to, weight)
    gu <- graphBAM(df, nodes="e", edgemode = "undirected")
    want <- list(a=c(b=1.5, c=5.6, d=2.1),
                 b=c(a=1.5, c=4.1),
                 c=c(a=5.6, b=4.1),
                 d=c(a=2.1, d=3.4),
                 e=numeric(0))
   checkEquals(want, edgeWeights(gu))

   checkEquals(want[c("c", "a")], edgeWeights(gu, c("c", "a")))
}


test_BAM_edges <- function() {
    g1 <- make_smallBAM()
    ew1 <- edges(g1)
    checkEquals(names(ew1), c("a", "b", "c", "x", "y"))
    checkEquals(list(ew1$a, ew1$b, ew1$c, ew1$x, ew1$y),
            list( c("b","c","x"), character(0), "a", c("c", "y"), character(0)))

    ew2 <- edges(g1, c("c", "b"))
    checkEquals(names(ew2), c("c","b"))
    checkEquals(list(ew2$c, ew2$b), list("a", character(0)))
}

test_BAM_adj <- function() {
    g1 <- make_smallBAM()
    ew <- adj(g1, c("c", "b"))
    checkEquals(names(ew), c("c","b"))
    checkEquals(list(ew$c, ew$b), list("a", character(0)))
}

test_BAM_edgeMatrix <- function() {
      g1 <- make_smallBAM()
      em <- edgeMatrix(g1)
      checkEquals(em[1,], c(3, 1, 1, 4, 1, 4))
      checkEquals(em[2,], c(1, 2, 3, 3, 4, 5))
}

test_BAM_adjacencyMatrix <- function() {
      g1 <- make_smallBAM()
      checkEquals(edgemode(g1), "directed")
      checkEquals(nodes(g1),  c("a","b","c","x","y"))
      am <- adjacencyMatrix(g1)
      checkEquals(rownames(am), nodes(g1))
      checkEquals(colnames(am), nodes(g1))
      checkEquals(as.integer(am["a",]), c(0, 1, 1, 1, 0))
      checkEquals(as.integer(am["b",]), c(0, 0, 0, 0, 0))
      checkEquals(as.integer(am["c",]), c(1, 0, 0, 0, 0))
      checkEquals(as.integer(am["x",]), c(0, 0, 1, 0, 1))
      checkEquals(as.integer(am["y",]), c(0, 0, 0, 0, 0))

}

test_BAM_removeEdge_unknown_nodes <- function()
{
    g1 <- make_smallBAM()
    checkException(removeEdge("a", "q", g1), silent=TRUE)
    checkException(removeEdge("q", "a", g1), silent=TRUE)
    checkException(removeEdge("a", c("q", "aa", "tt"), g1), silent=TRUE)
    checkException(removeEdge(c("a", "q", "tt", "aa"),
                              c("a", "q", "aa", "tt"), g1), silent=TRUE)
}

test_BAM_removeEdge <- function()
{
    g1 <- make_smallBAM()
    ## removing nothing does nothing
    c0 <- character(0)
    checkEquals(edges(g1), edges(removeEdge(c0, c0, g1)))
    ## there is no y => a edge, throw error
    checkException(removeEdge("y", "a", g1), silent=TRUE)

    g2 <- removeEdge("c", "a", g1)
    checkEquals(list(c=character(0)), edges(g2, "c"))
    em <- edgeMatrix(g2)
    checkEquals(em[1,], c(1, 1, 4, 1, 4))
    checkEquals(em[2,], c(2, 3, 3, 4, 5))

    g3 <- removeEdge("a", c("b", "x"), g1)
    checkEquals(list(a="c"), edges(g3, "a"))
    checkEquals(edges(g1)[-1], edges(g3)[-1])

    g4 <- removeEdge(c("a", "x"), "c", g1)
    checkEquals(list(a=c("b", "x")), edges(g4, "a"))
    checkEquals(list(x="y"), edges(g4, "x"))
}

test_BAMSmall_edgeData <- function(){
      g1 <- make_smallBAM()
      eg <- edgeData(g1)
      tmp <- paste(c("c", "a", "a", "x", "a", "x"), c("a","b","c","c","x","y"),sep="|")
      checkEquals(names(eg), tmp)
      vals <- sapply( names(eg),function(k){
               eg[[k]]$weight
              })
      checkEquals(names(vals), tmp)
      checkEquals( as.numeric(vals),c(7.9, 3.4, 2.6, 1.6, 1.7, 5.3))

      eg <- edgeData(g1, "a", attr="weight")
      tmp <- paste( c("a", "a", "a"), c("b", "c", "x"), sep = "|")
      checkEquals(names(eg), tmp)
      vals <- sapply( names(eg),function(k){
               eg[[k]] 
              })
      checkEquals(names(vals), tmp)
      checkEquals( as.numeric(vals), c(3.4, 2.6, 1.7))

      checkException(eg <- edgeData(g1, "a", attr="weightsss"), silent=TRUE)

      eg <- edgeData(g1, "a", "b", attr="weight")
      tmp <- paste("a", "b", sep = "|")
      checkEquals(names(eg), tmp)
      vals <- sapply( names(eg),function(k){
               eg[[k]]
              })
      checkEquals(names(vals), tmp)
      checkEquals( as.numeric(vals),3.4)
 }

test_BAM_extractFromToUndirected <- function() {
    g1 <- make_unDirectedBAM()
    ft <- extractFromTo(g1)
    checkEquals(as.character(ft$from), c("a", "a", "c", "a", "c", "x"))
    checkEquals(as.character(ft$to), c("b", "c", "d", "x", "x", "y"))
    checkEquals(ft$weight, c(3.4, 2.6, 7.9, 1.7, 1.6, 5.3))
}

test_BAM_extractFromToDirected <- function() {
    g1 <- make_smallBAM()
    ft <- extractFromTo(g1)
    checkEquals(as.character(ft$from), c("c", "a", "a", "x", "a", "x"))
    checkEquals(as.character(ft$to), c("a", "b", "c", "c", "x", "y"))
    checkEquals(ft$weight, c(7.9, 3.4, 2.6, 1.6, 1.7, 5.3))
}

test_BAM_bamToMatrix_UnDirected <- function() {
    g1 <- make_unDirectedBAM()
    mat <- as(g1, "matrix")
    checkEquals(isSymmetric(mat), TRUE)
    checkEquals(mat[upper.tri(mat)],
          c(3.4, 2.6, 0.0, 0.0, 0.0, 7.9, 1.7, 0.0,
                  1.6, 0.0, 0.0, 0.0, 0.0, 0.0, 5.3))
    checkEquals(rownames(mat),colnames(mat))
    checkEquals(rownames(mat), c("a", "b", "c", "d", "x", "y"))
}

test_BAM_bamToMatrix_Directed <- function() {
    g1 <- make_smallBAM()
    mat <- as(g1, "matrix")
    checkEquals(as.numeric(mat), c(0.0, 0.0, 7.9, 0.0,
                    0.0, 3.4, 0.0, 0.0, 0.0, 0.0, 2.6, 0.0,
                    0.0, 1.6, 0.0, 1.7, 0.0, 0.0, 0.0,0.0,
                    0.0, 0.0, 0.0, 5.3, 0.0))
    checkEquals(rownames(mat),colnames(mat))
    checkEquals(rownames(mat), c("a","b", "c", "x","y"))
}

test_BAM_bamTographAM_unDirected <- function() {
    g1 <- make_unDirectedBAM()
    am <- as(g1,"graphAM")
    checkEquals(nodes(g1), nodes(am))
    checkEquals(edgemode(g1), edgemode(am))
    checkEquals(edges(g1), edges(am))
    w1 <- edgeWeights(g1)
    w2 <- edgeWeights(am)
    checkEquals(names(w1), names(w2))
    checkEquals( w1$a, w2$a)
    checkEquals( w1$b, w2$b)
    checkEquals( sort(w1$c), sort(w2$c))
    checkEquals( w1$d, w2$d)
    checkEquals( sort(w1$x), sort(w2$x))
    checkEquals( w1$y, w2$y)
}

test_BAM_bamTographAM_Directed <- function() {
    g1 <- make_smallBAM()
    am <- as(g1,"graphAM")
    checkEquals(nodes(g1), nodes(am))
    checkEquals(edgemode(g1), edgemode(am))
    checkEquals(edges(g1), edges(am))
    w1 <- edgeWeights(g1)
    w2 <- edgeWeights(am)
    checkEquals(names(w1), names(w2))
    checkEquals( w1$a, w2$a)
    checkEquals( w1$b, w2$b)
    checkEquals( sort(w1$c), sort(w2$c))
    checkEquals( w1$d, w2$d)
    checkEquals( sort(w1$x), sort(w2$x))
    checkEquals( w1$y, w2$y)
}

test_BAM_bamTographNEL_UnDirected <- function() {
    g1 <- make_unDirectedBAM()
    nel <- as(g1,"graphNEL")
    checkEquals(nodes(g1), nodes(nel))
    checkEquals(edgemode(g1), edgemode(nel))
    checkEquals(edges(g1), edges(nel))
    w1 <- edgeWeights(g1)
    w2 <- edgeWeights(nel)
    checkEquals(names(w1), names(w2))
    checkEquals( w1$a, w2$a)
    checkEquals( w1$b, w2$b)
    checkEquals( sort(w1$c), sort(w2$c))
    checkEquals( w1$d, w2$d)
    checkEquals( sort(w1$x), sort(w2$x))
    checkEquals( w1$y, w2$y)
}


test_BAM_bamTographNEL_Directed <- function() {
    g1 <- make_smallBAM()
    nel <- as(g1,"graphNEL")
    checkEquals(nodes(g1), nodes(nel))
    checkEquals(edgemode(g1), edgemode(nel))
    checkEquals(edges(g1), edges(nel))
    w1 <- edgeWeights(g1)
    w2 <- edgeWeights(nel)
    checkEquals(names(w1), names(w2))
    checkEquals( w1$a, w2$a)
    checkEquals( w1$b, w2$b)
    checkEquals( sort(w1$c), sort(w2$c))
    checkEquals( w1$d, w2$d)
    checkEquals( sort(w1$x), sort(w2$x))
    checkEquals( w1$y, w2$y)
}

create_GraphNEL_Directed <- function() {
     set.seed(123)
     V <- letters[1:4]
     edL <- vector("list", length=4)
     names(edL) <- V
     edL[["a"]] <- list(edges=c(3, 4), weights=c(.13, .14))
     edL[["b"]] <- list(edges=c(3), weights=.23)
     edL[["c"]] <- list(edges=numeric(0), weights=numeric(0))
     edL[["d"]] <- list(edges=c(2, 3), weights=c(.42, .43))
     gR <- graphNEL(nodes = V, edgeL = edL, edgemode = "directed" )
     gR
}

create_GraphNEL_UnDirected <- function() {
     set.seed(123)
     V <- letters[1:4]
     edL <- vector("list", length=4)
     names(edL) <- V
     edL[["a"]] <- list(edges=c(2, 3), weights=c(.13, .14))
     edL[["b"]] <- list(edges=c(1), weights=.13)
     edL[["c"]] <- list(edges=c(1), weights=0.14)
     edL[["d"]] <- list(edges= numeric(0), weights=numeric(0))
     gR <- graphNEL(nodes = V, edgeL = edL, edgemode = "undirected" )
     gR
}

test_graphNEL_Directed_To_graphBAM <-function() {
    nel <- create_GraphNEL_Directed()
    bam <- as(nel, "graphBAM")
    checkEquals(nodes(nel), nodes(bam))
    checkEquals(edgemode(nel), edgemode(bam))
    checkEquals(edges(nel), edges(bam))
    w1 <- edgeWeights(nel)
    w2 <- edgeWeights(bam)
    checkEquals(w1,w2)
}

test_graphNEL_Directed_To_graphBAM <- function() {
    nel <- create_GraphNEL_Directed()
    bam <- as(nel, "graphBAM")
    checkEquals(nodes(nel), nodes(bam))
    checkEquals(edgemode(nel), edgemode(bam))
    checkEquals(edges(nel), edges(bam))
    w1 <- edgeWeights(nel)
    w2 <- edgeWeights(bam)
    checkEquals(w1,w2)
}

test_graphNEL_UnDirected_To_graphBAM <- function()  {
   nel <- create_GraphNEL_UnDirected()
   bam <- as(nel, "graphBAM")
   checkEquals(nodes(nel), nodes(bam))
   checkEquals(edgemode(nel), edgemode(bam))
   checkEquals(edges(nel), edges(bam))
   w1 <- edgeWeights(nel)
   w2 <- edgeWeights(bam)
   checkEquals(w1,w2)
}

test_graphAM_Directed_To_graphBAM <- function() {
    nel <- create_GraphNEL_Directed()
    am <- as(nel, "graphAM")
    bam <- as(am, "graphBAM")
    checkEquals(nodes(am), nodes(bam))
    checkEquals(edgemode(am), edgemode(bam))
    checkEquals(edges(am), edges(bam))
    w1 <- edgeWeights(am)
    w2 <- edgeWeights(bam)
    checkEquals(w1,w2)
}

test_graphAM_UnDirected_To_graphBAM<- function() {
   nel <- create_GraphNEL_UnDirected()
   am <- as(nel, "graphAM")
   bam <- as(am, "graphBAM")
   checkEquals(nodes(am), nodes(bam))
   checkEquals(edgemode(am), edgemode(bam))
   checkEquals(edges(am), edges(bam))
   w1 <- edgeWeights(am)
   w2 <- edgeWeights(bam)
   checkEquals(w1, w2)
}

test_BAM_set_edge_weights <- function()
{
    getw <- function(x) unlist(edgeWeights(x))

    g <- make_smallBAM()
    weight0 <- unlist(edgeWeights(g))
    edgeData(g, "c", "a", attr="weight") <- 123.0
    want <- weight0
    want["c.a"] <- 123.0
    checkEquals(want, getw(g))

    g <- make_smallBAM()
    edgeData(g, "a", c("b", "c", "x"), attr="weight") <- c(10, 11, 12)
    want <- weight0
    want[c("a.b", "a.c", "a.x")] <- c(10, 11, 12)
    checkEquals(want, getw(g))
}

test_BAM_Intersect_UnDirected <- function() {
    ## nodes a b c d x y
    from = c("a", "b", "d", "d")
    to   = c("b", "c", "x", "y")
    weight=c(1.2, 2.4, 3.2, 5.4)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "undirected")

    ## nodes a b c d x y z
    from = c("a", "b", "b", "d", "d")
    to   = c("b", "c", "d", "c", "x")
    weight=c(3.2, 1.2, 2.1, 3.2, 3.5)
    df <- data.frame(from, to, weight)
    g2 <- graphBAM(df, nodes = c("a","b","c", "d", "x", "y", "z"),
            edgemode = "undirected")

    g <- graphIntersect(g1,g2)
    checkEquals(intersect(nodes(g1), nodes(g2)), nodes(g))
    checkEquals(FALSE, isDirected(g))
    eg <- edgeData(g)
    vals <- sapply( names(eg),function(k){
               eg[[k]]$weight
              })
    tmp <- paste(c("a", "b", "d", "b", "c", "x"), c("b", "c", "x", "a", "b", "d"), sep= "|")
    checkEquals(tmp, names(vals))
    checkEquals(as.numeric(rep(NA, 6)), as.numeric(vals))
}


test_BAM_Intersect_Directed <- function() {
    ## nodes a b c d x y
    from = c("a", "b", "d", "d")
    to   = c("b", "c", "x", "y")
    weight=c(1.2, 2.4, 3.2, 5.4)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "directed")

    ## nodes a b c d x y z
    from = c("a", "b", "b", "d", "d")
    to   = c("b", "c", "d", "c", "x")
    weight=c(1.2, 1.2, 2.1, 3.2, 3.5)
    df <- data.frame(from, to, weight)
    g2 <- graphBAM(df, nodes = c("a","b","c", "d", "x", "y", "z"),
            edgemode = "directed")

    g <- graphIntersect(g1,g2)
    checkEquals(intersect(nodes(g1), nodes(g2)), nodes(g))
    checkEquals(TRUE, isDirected(g))
    eg <- edgeData(g)
    vals <- sapply( names(eg),function(k){
               eg[[k]]$weight
              })

    tmp <- paste(c("a", "b", "d"), c("b", "c", "x"), sep= "|")
    checkEquals(tmp, names(vals))
    checkEquals(c(1.2, NA, NA), as.numeric(vals))

}

test_BAM_Intersect_UnDirected2 <- function() {
    ## nodes a b d x y
    from = c("a", "d", "d")
    to   = c("b", "x", "y")
    weight=c(1.2, 3.2, 5.4)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "undirected")

    ## nodes a b c d x y z
    from = c("a", "b", "b", "d", "d")
    to   = c("b", "c", "d", "c", "x")
    weight=c(3.2, 1.2, 2.1, 5.2, 3.2)
    df <- data.frame(from, to, weight)
    g2 <- graphBAM(df, nodes = c("a","b","c", "d", "x", "y", "z"),
            edgemode = "undirected")

    g <- graphIntersect(g1,g2)
    checkEquals(intersect(nodes(g1), nodes(g2)), nodes(g))
    checkEquals(FALSE, isDirected(g))
    eg <- edgeData(g)
    vals <- sapply( names(eg),function(k){
               eg[[k]]$weight
              })
    tmp <- paste(c("a", "d", "b", "x"), c("b", "x", "a", "d"), sep= "|")
    checkEquals(tmp, names(vals))
    checkEquals(rep(c(NA,3.2),2), as.numeric(vals))
}

test_BAM_Intersect_EmptyEdges <- function() {

    from = c("a", "d", "d")
    to   = c("b", "x", "y")
    weight=c(1.2, 3.2, 5.4)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "directed")

    from = c("h", "i", "j")
    to   = c("b", "x", "y")
    weight=c(1.2, 3.2, 5.4)
    df <- data.frame(from, to, weight)
    g2 <- graphBAM(df, edgemode = "directed")

    g <- graphIntersect(g1,g2)
    checkEquals(nodes(g), intersect(nodes(g1), nodes(g2)))
    checkEquals(isDirected(g), TRUE)
    eg <- edgeWeights(g)
    checkEquals(c("b", "x", "y"), names(eg))
    checkEquals(list(numeric(0), numeric(0), numeric(0)),list(eg$b, eg$x, eg$y))
}

test_BAM_Intersect_EmptyNodes <- function() {

    from = c("a", "d", "d")
    to   = c("b", "x", "y")
    weight=c(1.2, 3.2, 5.4)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "unirected")

    from = c("h", "i", "j")
    to   = c("s", "h", "l")
    weight=c(1.2, 3.2, 5.4)
    df <- data.frame(from, to, weight)
    g2 <- graphBAM(df, edgemode = "undirected")

    g <- graphIntersect(g1,g2)
    checkEquals(intersect(nodes(g1), nodes(g2)), nodes(g))
    checkEquals(FALSE, isDirected(g))
    eg <- edgeWeights(g)
    checkEquals(list(), eg)
}

test_BAM_isAdjacent <- function()
{
     from = c("a", "d", "d", "b", "a")
     to   = c("b", "a", "d", "c", "c")
     weight= c(1.5, 2.1, 3.4, 4.1, 5.6)
     df <- data.frame(from, to, weight)
     gd <- graphBAM(df, nodes="e", edgemode = "directed")

     ## single edges
     for (i in seq_len(nrow(df))) {
         checkEquals(TRUE, isAdjacent(gd, from[i], to[i]))
     }

     ## vectorized
     checkEquals(c(FALSE, TRUE, TRUE, FALSE, FALSE),
                 isAdjacent(gd, "a", letters[1:5]))

     checkEquals(c(FALSE, FALSE, FALSE, TRUE, FALSE),
                 isAdjacent(gd, letters[1:5], "a"))
}

test_BAM_Union_UnDirected <- function() {
    ## nodes a b c d x y
    from = c("a", "b", "d", "d")
    to   = c("b", "c", "x", "y")
    weight=c(1.2, 2.4, 3.5, 5.4)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "undirected")

    ## nodes a b c d x y z 
    from = c("a", "b", "b", "d", "d")
    to   = c("b", "c", "d", "c", "x")
    weight=c(3.2, 1.2, 2.1, 3.2, 3.5)
    df <- data.frame(from, to, weight)
    g2 <- graphBAM(df, nodes = c("a","b","c", "d", "x", "y", "z"), 
        edgemode = "undirected")
    g <- graphUnion(g1,g2)
    checkEquals(union(nodes(g1), nodes(g2)), nodes(g))
    checkEquals(FALSE, isDirected(g))
    df <- extractFromTo(g)
    tmp <- data.frame(from = c("a", "b", "b", "c", "d", "d"),
        to = c("b", "c", "d", "d", "x", "y"),
        weight = c( NA, NA, 2.1, 3.2, 3.5, 5.4))
    checkEquals(tmp, df)
}


test_BAM_Union_Directed <- function() {
    ## nodes a b c d x y
    from = c("a", "b", "d", "d")
    to   = c("b", "c", "x", "y")
    weight=c(1.2, 2.4, 3.5, 5.4)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "directed")

    ## nodes a b c d x y z 
    from = c("a", "b", "b", "d", "d")
    to   = c("b", "c", "d", "c", "x")
    weight=c(1.2, 1.2, 2.1, 3.2, 3.5)
    df <- data.frame(from, to, weight)
    g2 <- graphBAM(df, nodes = c("a","b","c", "d", "x", "y", "z"), 
            edgemode = "directed")

    g <- graphUnion(g1,g2)
    checkEquals(union(nodes(g1), nodes(g2)), nodes(g))
    checkEquals(TRUE, isDirected(g))

    df <- extractFromTo(g)
    tmp <- data.frame(from = c("a", "b", "d", "b", "d", "d"),
                        to = c("b", "c", "c", "d", "x", "y"),
                    weight = c( 1.2, NA, 3.2, 2.1, 3.5, 5.4))
    checkEquals(tmp, df)

}

test_BAM_Union_Mixed <- function() {
    ## nodes a b d x y
    from = c("a", "d", "d")
    to   = c("b", "x", "y")
    weight=c(1.2, 3.2, 5.4)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "undirected")

    ## nodes a b c d x y z 
    from = c("a", "b", "b", "d", "d")
    to   = c("b", "c", "d", "c", "x")
    weight=c(3.2, 1.2, 2.1, 3.2, 3.5)
    df <- data.frame(from, to, weight)
    g2 <- graphBAM(df, nodes = c("a","b","c", "d", "x", "y", "z"), 
            edgemode = "directed")

    checkException(g <- graphUnion(g1,g2), silent=TRUE)
}

test_BAM_inEdges <- function()
{
      from = c("a", "d", "d", "b", "a")
        to = c("b", "a", "d", "c", "c")
    weight = c(1.5, 2.1, 3.4, 4.1, 5.6)
      df <- data.frame(from, to, weight)
      ## directed
      gd <- graphBAM(df, nodes="e", edgemode = "directed")

      want <- list(a="d",
                   b="a",
                   c=c("a", "b"),
                   d="d",
                   e=character(0))
      checkEquals(want, inEdges(nodes(gd), gd))

      ## undirected
      gu <- graphBAM(df, nodes="e", edgemode = "undirected")
      checkEquals(edges(gu), inEdges(nodes(gu), gu))
}


test_BAM_directed_attrs <- function() {

    from = c("a", "a", "a", "x", "x", "c")
    to   = c("b", "c", "x", "y", "c", "a")
    weight = c(2, 1, 3, 4, 5, 6)
    df <- data.frame(from, to, weight)
    bam <- graphBAM(df, edgemode = "directed")
    
    checkException(edgeData(bam,from="a", attr="code"), silent=TRUE)
    edgeDataDefaults(bam, attr ="weight") <- 1
    edgeDataDefaults(bam, attr = "code") <- "plain"

    res <- unlist(edgeData(bam,from="a", attr="code"))
    nmres <- paste(c("a","a","a"), c ("b", "c", "x"), sep="|")  
    checkEquals(names(res), nmres)
    checkEquals(as.character(res), c("plain", "plain", "plain"))

    edgeData(bam,from = "a", to = "x", attr= "code") <- "red"
    res <- unlist(edgeData(bam, from = "a", attr = "code"))
    checkEquals(names(res), nmres)
    checkEquals(as.character(res), c("plain", "plain", "red"))

    edgeData(bam,to = "c", attr= "code") <- "yellow"
    res <- unlist(edgeData(bam, to= "c", attr = "code"))
    nmres <- paste(c("a", "x"), c("c", "c"), sep = "|")
    checkEquals(names(res), nmres)
    checkEquals(as.character(res), c("yellow", "yellow"))
}

test_BAM_undirected_attrs <- function() {

    from = c("a", "a", "a", "x", "x")
    to   = c("b", "c", "x", "y", "c")
    weight = c(2, 1, 3, 4, 5)
    df <- data.frame(from, to, weight)
    bam <- graphBAM(df, edgemode = "undirected")
    checkException(edgeData(bam,from="a", attr="code"), silent=TRUE)

    edgeDataDefaults(bam, attr = "weight") <- 1
    edgeDataDefaults(bam, attr = "code") <- "plain"

    res <- unlist(edgeData(bam,from="a", attr="code"))
    nmres <- paste(c("a","a","a"), c ("b", "c", "x"), sep="|")  
    checkEquals(names(res), nmres)
    checkEquals(as.character(res), c("plain", "plain", "plain"))

    edgeData(bam,from = "a", to = "x", attr= "code") <- "red"
    res <- unlist(edgeData(bam, from = "a", attr = "code"))
    checkEquals(names(res), nmres)
    checkEquals(as.character(res), c("plain", "plain", "red"))

    edgeData(bam,to = "c", attr= "code") <- "yellow"
    res <- unlist(edgeData(bam, to= "c", attr = "code"))
    nmres <- paste(c("a", "x"), c("c", "c"), sep = "|")
    checkEquals(names(res), nmres)
    checkEquals(as.character(res), c("yellow", "yellow"))
}


test_graphBAM_detailed_Attribute_Intersection <- function() {

    ## nodes a b c d x y
    from = c("a", "b", "d", "d")
    to   = c("b", "c", "y", "x")
    weight=c(1.2, 2.4, 5.4, 3.2)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "directed")
    edgeData(g1, from = from, to = to ,attr = "weight")  <- c(1.2, 2.4, 5.4, 3.2)


    edgeDataDefaults(g1, attr = "color") <- "unknown"
    edgeDataDefaults(g1, attr ="type") <- "unknown"
    edgeData(g1, from = from, to = to ,attr = "color") <-  c("red", "blue", NA, "green")
    edgeData(g1, from = from, to = to , attr = "type") <-  c("high", "low", "high", NA)
    ## nodes a b c d x y z
    from = c("a", "b", "b", "d", "d")
    to   = c("b", "c", "d", "c", "x")
    weight=c(1.2, 4.2, 5.6, 2.1, 3.2)
    df <- data.frame(from, to, weight)
    g2 <- graphBAM(df, nodes = c("a","b","c", "d", "x", "y", "z"),
        edgemode = "directed")
    edgeDataDefaults(g2, attr = "color") <- "unknown"
    edgeData(g2, from = from, to = to,  attr = "color") <- c("red", "blue", NA, "red",
        "yellow")
    g <- graphIntersect(g1, g2)
    df <- extractFromTo(g)
    tmp <- data.frame( from = c("a", "b", "d"), 
        to = c("b", "c", "x"), 
        weight = c(1.2, NA, 3.2))

    checkEquals(tmp, df)

    attColor <- edgeData(g, attr = "color")
    nms <- paste(c("a", "b", "d"),  c("b", "c", "x"), sep = "|")
    target <- structure( c("red", "blue", NA), names = nms)
    checkEquals(target, unlist(attColor))

    checkException(edgeData(g, attr = "type"), silent=TRUE)

    weightFun <- function(x, y) {
        return(x +y )
    }

    colorFun <- function(x,y) {
        if(x=="red" && y == "red")
            return("white")
        else
            return("black")
        }

    setClass("myType", representation = representation(typ ="character")) 
    myType <- function(typ){ new("myType", typ = typ)}
    typeFun <- function(x,y) {
            if(is(x, "myType")  && is(y, "myType")){
                if(x@typ =="low" || y@typ == "med")
                    return("low")
                else
                    return("high")
                }
            else {return (NA)}
    }
    nodeDataDefaults(g1, attr ="color") <- "unknown"
    nodeDataDefaults(g1, attr ="type") <- "unknown"
    nodeDataDefaults(g2, attr ="color") <- "unknown"
    nodeDataDefaults(g2, attr ="type") <- "unknown"

    nodeData(g1,n = c("a", "b", "c"), attr ="color") <- c("red", "green", "blue")
    nodeData(g1,n = c("b", "c"), attr ="type") <- c(myType("low"), myType("high"))
    nodeData(g2,n = c("a", "b", "c"), attr ="color") <- c("red", "green", "red")
    nodeData(g2,n = c("b", "c"), attr ="type") <- c(myType("med"), myType("low"))
    g <- graphIntersect(g1, g2, nodeFun = list(type = typeFun),
        edgeFun = list(weight = weightFun, color = colorFun))

    attWeight <- edgeData(g, attr = "weight")
    nms <- paste(c("a", "b", "d"),  c("b", "c", "x"), sep = "|")
    target <- structure( c( 2.4, 6.6, 6.4), names = nms)
    checkEquals(target, unlist(attWeight))

    attColor <- edgeData(g, attr = "color")
    nms <- paste(c("a", "b", "d"),  c("b", "c", "x"), sep = "|")
    target <- structure( c( 2.4, 6.6, 6.4), names = nms)
    checkEquals(target, unlist(attWeight))

    nodeColor <- nodeData(g, attr = "color")
    target <-  as.list(structure(c("red", "green", NA, "unknown", "unknown",
                "unknown"), names = c("a", "b", "c", "d", "x", "y")))
    checkEquals(target, nodeColor)

    nodeType <- nodeData(g, attr = "type")
    target <-  as.list(structure(c("unknown", "low", "high", "unknown",
                "unknown", "unknown"), names = c("a", "b", "c", "d", "x", "y")))
    checkEquals(target, nodeType)
}

test_graphBAM_detailed_Attribute_Union <- function() {

    ## nodes a b c d x y
    from = c("a", "b", "d", "d")
    to   = c("b", "c", "y", "x")
    weight=c(1.2, 2.4, 5.4, 3.2)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "directed")
    edgeData(g1, from = from, to = to ,attr = "weight")  <- c(1.2, 2.4, 5.4, 3.2)
    
    edgeDataDefaults(g1, attr = "color")  <- "cyan"
    edgeDataDefaults(g1, attr = "type")  <- "unknown"
    edgeData(g1, from = from, to = to ,attr = "color") <-  c("red", "blue", NA, "green")
    edgeData(g1, from = from, to = to , attr = "type") <-  c("high", "low", "high", NA)
  
    ## nodes a b c d x y z
    from = c("a", "b", "b", "d", "d")
    to   = c("b", "c", "d", "c", "x")
    weight=c(1.2, 4.2, 5.6, 2.1, 3.2)
    df <- data.frame(from, to, weight)
    g2 <- graphBAM(df, nodes = c("a","b","c", "d", "x", "y", "z"),
        edgemode = "directed")
    edgeDataDefaults(g2, attr = "color")  <- "cyan"

    edgeData(g2, from = from, to = to,  attr = "color") <- c("red", "blue", NA, "red",
        "yellow")
    g <- graphUnion(g1, g2)
    df <- extractFromTo(g)
    tmp <- data.frame( from = c("a", "b", "d", "b", "d", "d"), 
        to = c("b", "c", "c", "d", "x", "y"), 
        weight = c(1.2, NA, 2.1, 5.6, 3.2, 5.4))
    checkEquals(tmp, df)

    attColor <- edgeData(g, attr = "color")
    nms <- paste(c("a", "b", "d", "b", "d", "d"),  c("b", "c", "c", "d", "x", "y"), sep = "|")
    target <- structure( c("red", "blue", "red", NA, NA, NA), names = nms)
    checkEquals(target, unlist(attColor))

    attType <- edgeData(g, attr = "type")
    nms <- paste(c("a", "b", "d", "b", "d", "d"),  c("b", "c", "c", "d", "x", "y"), sep = "|")
    target <- structure( c("high", "low", NA, NA, NA, "high"), names = nms)
    checkEquals(target, unlist(attType))


    weightFun <- function(x, y) {
        return(x + y )
    }

    colorFun <- function(x,y) {
        if(x=="red" || y == "red")
            return("white")
        else
            return("black")
    }

    setClass("myType", representation = representation(typ ="character")) 
    myType <- function(typ){ new("myType", typ = typ)}
    typeFun <- function(x,y) {
        if(is(x, "myType")  && is(y, "myType")){
            if(x@typ =="low" || y@typ == "med")
                return("low")
            else
                return("high")
            }
        else {return (NA)}
    }
   
    nodeDataDefaults(g1, attr ="color") <- "cyan"
    nodeDataDefaults(g1, attr="type") <- "unknown"
    nodeData(g1,n = c("a", "b", "c"), attr ="color") <- c("red", "green", "blue")
    nodeData(g1,n = c("b", "c"), attr ="type") <- c(myType("low"), myType("high"))
    
    nodeDataDefaults(g2, attr ="color") <- "cyan"
    nodeDataDefaults(g2, attr="type") <- "unknown"
    nodeDataDefaults(g2, attr="test") <- "missing"

    nodeData(g2,n = c("a", "b", "c", "z"), attr ="color") <- c("red", "green", "red","pink")
    nodeData(g2,n = c("b", "c"), attr ="type") <- c(myType("med"), myType("low"))
    nodeData(g2,n = c("a", "b", "c"), attr = "test") <- c("pass", "fail", "pass")


    g <- graphUnion(g1, g2, edgeFun = list(weight = weightFun, color = colorFun))

    attWeight <- edgeData(g, attr = "weight")
    nms <- paste(c("a", "b", "d", "b", "d", "d"),  c("b", "c", "c", "d", "x", "y"), sep = "|")
    target <- structure( c( 2.4, 6.6, 2.1, 5.6, 6.4, 5.4), names = nms)
    checkEquals(target, unlist(attWeight))

    attColor <- edgeData(g, attr = "color")
    nms <- paste(c("a", "b", "d", "b", "d", "d"),  c("b", "c", "c", "d", "x", "y"), sep = "|")
    target <- structure(c( "white", "black", "red", NA, "black", NA), names = nms)
    checkEquals( target, unlist(attColor))

    attType <- edgeData(g, attr = "type")
    nms <- paste(c("a", "b", "d", "b", "d", "d"),  c("b", "c", "c", "d", "x", "y"), sep = "|")
    target <- structure( c("high", "low", NA, NA, NA, "high"), names = nms)
    checkEquals(target, unlist(attType))


    attType <- edgeData(g, attr = "type")
    nms <- paste(c("a", "b", "d", "b", "d", "d"),  c("b", "c", "c", "d", "x", "y"), sep = "|")
    target <- structure( c("high", "low", NA, NA, NA, "high"), names = nms)
    checkEquals(target, unlist(attType))

}


test_graphBAM_removeEdgesByWeight <- function() {
    from = c("a", "b", "d", "d")
    to   = c("b", "c", "y", "x")
    weight=c(2.2, 2.0, 0.4, 0.2)
    df <- data.frame(from, to, weight)
    g <- graphBAM(df, edgemode = "directed")

    edgeDataDefaults(g, attr="color") <- "pink"
    edgeData(g, from = from, to = to ,attr = "color") <-  c("red", "blue", NA, "green")

    res <- removeEdgesByWeight(g, lessThan = 2.0)
    checkEquals(attr(res@edgeSet@bit_vector, "nbitset"), 2)
    checkEquals(res@edgeSet@weights, c(2.2, 2.0))
    current <- unlist( edgeData(res, attr = "color"))
    target <- structure(c("red", "blue"), 
        names = paste(c("a", "b"), c("b", "c"), sep = "|"))
    checkEquals(target, current)

    res <- removeEdgesByWeight(g, greaterThan = 1.9)
    checkEquals(attr(res@edgeSet@bit_vector, "nbitset"), 2)
    checkEquals(res@edgeSet@weights, c(0.2, 0.4))
    current <- unlist( edgeData(res, attr = "color"))
    target <- structure(c("green", NA), 
        names = paste(c("d", "d"), c("x", "y"), sep = "|"))
    checkEquals(target, current)

    res <- removeEdgesByWeight(g, lessThan =1.0, greaterThan = 2)
    checkEquals(res@edgeSet@weights, c(2.0))
    current <- unlist( edgeData(res, attr = "color"))
    target <- structure(c("blue"), 
        names = paste(  "b", "c", sep = "|"))
    checkEquals(target, current)

    res <- removeEdgesByWeight(g, greaterThan = 0.1)
    checkEquals(res@edgeSet@weights, numeric(0))
    checkEquals(res@edgeSet@edge_attrs$color, character(0))
}

test_graphBAM_nodeAttributes <- function(){

    from = c("a", "b", "d", "d")
    to   = c("b", "c", "y", "x")
    weight=c(2.2, 2.0, 0.4, 0.2)
    df <- data.frame(from, to, weight)
    g <- graphBAM(df, edgemode = "directed")
    nodeDataDefaults(g, attr ="color") <- "blue"
    sg <- subGraph(c("a", "c"), g)
    checkIdentical(unname(unlist(nodeData(sg))), c("blue", "blue"))

    nodeData(g, n = c("d","a"), attr = "color") <- c("red", "green")
    current <- nodeData(g, attr = "color")
    target <- as.list(structure(  c("green", "blue", "blue", "red", "blue", "blue"), 
                        names = c("a", "b", "c", "d", "x", "y")))
    checkEquals(target, current)
 
    nodeDataDefaults(g, attr="mat") <- NA 
    nodeData(g, n= c("x", "y"), attr = "mat") <- df
    current <-  nodeData(g, n= c("x", "y"), attr = "mat")
    target <- list(x = df, y = df)
    checkEquals(target, current)

    sg <- subGraph(c("d","b"), g)
    current <- nodeData(sg, attr = "color")
    target  <- as.list(structure(c("blue", "red"), names = c("b", "d")))
    checkEquals(target, current)
}


test_BAM_directed_attrs_s4 <- function() {

    from = c("a", "a", "a", "x", "x", "c")
    to   = c("b", "c", "x", "y", "c", "a")
    weight = c(2, 1, 3, 4, 5, 6)
    df <- data.frame(from, to, weight)
    bam <- graphBAM(df, edgemode = "directed")
    
    edgeDataDefaults(bam, attr = "weight") <- 1.3
    edgeDataDefaults (bam, attr = "vals") <- df
    edgeData(bam, from = "a", attr= "vals") <- "unknown"

    res <- edgeData(bam, attr="vals")
    nmres <- c("c|a", "a|b", "a|c", "x|c", "a|x", "x|y")
    target <- structure(list(df, "unknown", "unknown", df, "unknown",df), names = nmres)
    checkEquals(res, target)

    edgeDataDefaults(bam, attr = "mat") <- NA
    edgeData(bam,from = "a", to = "x", attr= "mat") <- matrix(1)
    res <- edgeData(bam, from = "a", attr = "mat")
    nmres <- paste(c("a", "a", "a"), c("b", "c", "x"), sep = "|")
    target <- structure( list(NA, NA, matrix(1)), names = nmres)
    checkEquals(res, target)

    edgeDataDefaults(bam, attr = "mk") <- NA
    edgeData(bam,to = "c", attr= "mk") <- matrix(1)
    res <- edgeData(bam, attr = "mk")
    nmres <- paste(c("c", "a", "a", "x", "a", "x"), c("a", "b", "c", "c", "x", "y"), sep ="|")
    target <- structure( list(NA, NA, matrix(1), matrix(1), NA ,NA), names = nmres)
    checkEquals(res, target)
}

test_BAM_undirected_attrs_s4 <- function() {

    from = c("a", "a", "a", "x")
    to   = c("b", "c", "x", "y")
    weight = c(2, 1, 3, 4)
    df <- data.frame(from, to, weight)
    bam <- graphBAM(df, edgemode = "undirected")
   
    edgeDataDefaults(bam, attr = "weight") <- 1.3
    edgeDataDefaults(bam, attr = "vals") <- df

    #  edgeData(bam, attr = "weight") <- 1.3
    # edgeData(bam, attr = "vals") <- df
    edgeData(bam, from = "x", attr = "vals") <- "unknown"

    res <- edgeData(bam, attr="vals")
    nmres <- c("a|b", "a|c", "a|x", "x|y", "b|a", "c|a", "x|a", "y|x")
    target <- structure(list(df, df, "unknown", "unknown", df, df, "unknown", 
                    "unknown"), names = nmres)
    checkEquals(res, target)
  
    edgeDataDefaults(bam, attr ="mat") <- NA 
    edgeData(bam,from = "a", to = "x", attr= "mat") <-  matrix(1)
    res <- edgeData(bam, attr = "mat")
    target <- structure(list(NA, NA, matrix(1), NA, NA, NA, matrix(1), NA), 
            names = nmres)
    checkEquals(res, target)

    edgeDataDefaults(bam, attr = "mk") <- NA
    edgeData(bam,to = "c", attr= "mk") <- matrix(1)
    res <- edgeData(bam, attr = "mk")
    target <- structure( list(NA, matrix(1), NA, NA, NA, matrix(1), NA ,NA), 
            names = nmres)
    checkEquals(res, target)
}


test_graphBAM_S4_Attribute_Intersection <- function() {

    setClass("myColor", representation = representation(col ="character")) 
    setClass("myType", representation = representation(typ ="character")) 
    myColor <- function(col){ new("myColor", col = col)}
    myType <- function(typ){ new("myType", typ = typ)}

    ## nodes a b c d x y
    from = c("a", "b", "d", "d")
    to   = c("b", "c", "y", "x")
    weight=c(1.2, 2.4, 5.4, 3.2)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "directed")
    edgeData(g1, from = from, to = to ,attr = "weight")  <- c(1.2, 2.4, 5.4, 3.2)

    edgeDataDefaults(g1, attr="color") <- "cyan"
    edgeDataDefaults(g1, attr="type") <- "unknown"
    edgeData(g1, from = from, to = to ,attr = "color") <-  c(myColor("red"),
        myColor("blue"), NA, myColor("green"))

    edgeData(g1, from = from, to = to , attr = "type") <- c(myType("high"),
        myType("low"), myType("high"), NA)
    ## nodes a b c d x y z
    from = c("a", "b", "b", "d", "d")
    to   = c("b", "c", "d", "c", "x")
    weight=c(1.2, 4.2, 5.6, 2.1, 3.2)
    df <- data.frame(from, to, weight)
    g2 <- graphBAM(df, nodes = c("a","b","c", "d", "x", "y", "z"),
        edgemode = "directed")
    edgeDataDefaults(g2, attr ="color") <- "cyan"
    edgeData(g2, from = from, to = to,  attr = "color") <- c(myColor("red"),
        myColor("blue"), NA, myColor("red"), myColor("yellow"))

    g <- graphIntersect(g1, g2)
    df <- extractFromTo(g)
    tmp <- data.frame( from = c("a", "b", "d"), to = c("b", "c", "x"), 
        weight = c(1.2, NA, 3.2))
    checkEquals(tmp, df)

    attColor <- edgeData(g, attr = "color")
    nms <- paste(c("a", "b", "d"),  c("b", "c", "x"), sep = "|")
    target <- structure( c(myColor("red"), myColor("blue"), NA), names = nms)
    checkEquals(target, unlist(attColor))

    checkException(edgeData(g, attr = "type"), silent=TRUE)

    weightFun <- function(x, y) {
        return(x + y )
    }
    colorFun <- function(x,y) {
        if(x@col=="red" && y@col == "red")
            return("white")
        else
            return("black")
    }
    g <- graphIntersect(g1, g2, edgeFun =list(weight = weightFun, color = colorFun))

    df <- extractFromTo(g)
    tmp <- data.frame( from = c("a", "b", "d"), 
            to = c("b", "c", "x"), 
            weight = c(2.4, 6.6 , 6.4))
    checkEquals(tmp, df)
    attColor <- edgeData(g, attr = "color")
    nms <- paste(c("a", "b", "d"),  c("b", "c", "x"), sep = "|")
    target <- structure( c("white", "black", "black"), names = nms)
    checkEquals(target, unlist(attColor))
    checkException(edgeData(g, attr = "type"), silent=TRUE)

}

test_graphBAM_S4_Attribute_Union <- function() {

    setClass("myColor", representation = representation(col ="character")) 
    setClass("myType", representation = representation(typ ="character")) 
    myColor <- function(col){ new("myColor", col = col)}
    myType <- function(typ){ new("myType", typ = typ)}

    ## nodes a b c d x y
    from = c("a", "b", "d", "d")
    to   = c("b", "c", "y", "x")
    weight=c(1.2, 2.4, 5.4, 3.2)
    df <- data.frame(from, to, weight)
    g1 <- graphBAM(df, edgemode = "directed")
    edgeData(g1, from = from, to = to ,attr = "weight")  <- c(1.2, 2.4, 5.4, 3.2)

    edgeDataDefaults(g1, attr = "color") <- "cyan"
    edgeDataDefaults(g1, attr = "type") <- "missing"
    edgeData(g1, from = from, to = to ,attr = "color") <- c(myColor("red"), 
        myColor("blue"), NA, myColor("green"))
    edgeData(g1, from = from, to = to , attr = "type") <- c(myType("high"), 
        myType("low"), myType("high"), NA)
    ## nodes a b c d x y z
    from = c("a", "b", "b", "d", "d")
    to   = c("b", "c", "d", "c", "x")
    weight=c(1.2, 4.2, 5.6, 2.1, 3.2)
    df <- data.frame(from, to, weight)
    g2 <- graphBAM(df, nodes = c("a","b","c", "d", "x", "y", "z"),
            edgemode = "directed")
    edgeDataDefaults(g2, attr = "color") <- "cyan"
    edgeData(g2, from = from, to = to,  attr = "color") <- c(myColor("red"), 
        myColor("blue"), NA, myColor("red"), myColor("yellow"))

    g <- graphUnion(g1, g2)
    df <- extractFromTo(g)
    tmp <- data.frame( from = c("a", "b", "d", "b", "d", "d"), 
            to = c("b", "c", "c", "d", "x", "y"), 
            weight = c(1.2, NA, 2.1, 5.6, 3.2, 5.4))
    checkEquals(tmp, df)

    attColor <- edgeData(g, attr = "color")
    nms <- paste(c("a", "b", "d", "b", "d", "d"),  c("b", "c", "c", "d", "x", "y"), sep = "|")
    target <- structure( c(myColor("red"), myColor("blue"), myColor("red"), NA, NA, NA), names = nms)
    checkEquals(target, unlist(attColor))

    attType <- edgeData(g, attr = "type")
    nms <- paste(c("a", "b", "d", "b", "d", "d"),  c("b", "c", "c", "d", "x", "y"), sep = "|")
    target <- structure( c(myType("high"), myType("low"), NA, NA, NA, myType("high")), names = nms)
    checkEquals(target, unlist(attType))

    weightFun <- function(x, y) {
        return(x + y )
    }

    colorFun <- function(x,y) {
        if(x@col =="red" || y@col == "red")
            return("white")
        else
            return("black")
    }

    g <- graphUnion(g1, g2, edgeFun = list(weight = weightFun, color = colorFun))

    attWeight <- edgeData(g, attr = "weight")
    nms <- paste(c("a", "b", "d", "b", "d", "d"),  c("b", "c", "c", "d", "x", "y"), sep = "|")
    target <- structure( c( 2.4, 6.6, 2.1, 5.6, 6.4, 5.4), names = nms)
    checkEquals(target, unlist(attWeight))

    attColor <- edgeData(g, attr = "color")
    nms <- paste(c("a", "b", "d", "b", "d", "d"),  c("b", "c", "c", "d", "x", "y"), sep = "|")
    target <- structure(c( "white", "black", myColor("red"), NA, "black", NA), names = nms)
    checkEquals( target, unlist(attColor))

    attType <- edgeData(g, attr = "type")
    nms <- paste(c("a", "b", "d", "b", "d", "d"),  c("b", "c", "c", "d", "x", "y"), sep = "|")
    target <- structure( c(myType("high"), myType("low"), NA, NA, NA, myType("high")), names = nms)
    checkEquals(target, unlist(attType))

    attType <- edgeData(g, attr = "type")
    nms <- paste(c("a", "b", "d", "b", "d", "d"),  c("b", "c", "c", "d", "x", "y"), sep = "|")
    target <- structure(c( myType("high"), myType("low"), NA, NA, NA, myType("high")), names = nms)
    checkEquals(target, unlist(attType))

}

test_graphBAM_addNode1 <- function(){

    from = c("a", "b", "d", "d")
    to   = c("b", "c", "y", "x")
    weight=c(2.2, 2.0, 0.4, 0.2)
    df <- data.frame(from, to, weight)
    g <- graphBAM(df, edgemode = "directed")
    nodeDataDefaults(g, attr="color") <- "pink"
    nodeData(g, n = c("d","a"), attr = "color") <- c("red", "green")
    nodeDataDefaults(g, attr="type") <- "unknown"
    nodeData(g, n = c("a", "b", "y", "d"), attr = "type") <- c("high", "med", "high", "low") 

    gr <- addNode(c("q", "ss"), g)
    current <- nodeData(gr, attr = "color")

    target <- c("green", "pink", "pink", "red", "pink", "pink", "pink", "pink")
    names(target) <- c("a", "b", "c", "d", "q", "ss", "x", "y")
    checkTrue(all(current[sort(names(current))] == target[sort(names(target))]))

    current <- nodeData(gr, attr = "type")
    target <- c("high", "med", "unknown", "low", "unknown", "unknown", "unknown", "high")
    names(target) <- c("a", "b", "c", "d", "q", "ss", "x", "y")
    checkTrue(all(current[sort(names(current))] == target[sort(names(target))]))

}

test_graphBAM_addNode1 <- function(){

    from = c("a", "b", "d", "d")
    to   = c("b", "c", "y", "x")
    weight=c(2.2, 2.0, 0.4, 0.2)
    df <- data.frame(from, to, weight)
    g <- graphBAM(df, edgemode = "directed")
    nodeDataDefaults(g, attr="color") <- "pink"
    nodeData(g, n = c("d","a"), attr = "color") <- c("red", "green")
    nodeDataDefaults(g, attr="type") <- "unknown"
    nodeData(g, n = c("a", "b", "y", "d"), attr = "type") <- c("high", "med", "high", "low") 

    gr <- addNode(c("q", "ss"), g)
    current <- nodeData(gr, attr = "color")

    target <- c("green", "pink", "pink", "red", "pink", "pink", "pink", "pink")
    names(target) <- c("a", "b", "c", "d", "q", "ss", "x", "y")
    checkTrue(all(current[sort(names(current))] == target[sort(names(target))]))

    current <- nodeData(gr, attr = "type")
    target <- c("high", "med", "unknown", "low", "unknown", "unknown", "unknown", "high")
    names(target) <- c("a", "b", "c", "d", "q", "ss", "x", "y")
    checkTrue(all(current[sort(names(current))] == target[sort(names(target))]))

}

test_graphBAM_addNode2 <- function(){

    from = c("a", "b", "d", "d")
    to   = c("b", "c", "y", "x")
    weight=c(2.2, 2.0, 0.4, 0.2)
    df <- data.frame(from, to, weight)
    g <- graphBAM(df, edgemode = "directed")
    
    edgeDataDefaults(g, attr="color") <- "blue"
    edgeDataDefaults(g, attr="type") <- "unknown"
    edgeData(g, from = c("d","a"), to = c("y", "b"), attr = "color") <- c("red", "green")
    edgeData(g, from = c("a", "b"), to = c("b", "c") , attr = "type") <- c("low", "high")
    g1 <- addEdge(from = c("d", "b"), to = c("c", "x"), g, weights = c(4.0, 10.0))

    current <- edgeData(g1, attr ="weight")
    lbl <- paste(c("a", "b", "d", "b", "d", "d"), c( "b", "c", "c", "x", "x", "y") , sep ="|")
    target <- as.list( structure(c(2.2, 2, 4, 10, 0.2, 0.4), names = lbl))
    checkEquals(target, current)
      
    current <- edgeData(g1, attr ="color")
    lbl <- paste(c("a", "b", "d", "b", "d", "d"), 
                c( "b", "c", "c", "x", "x", "y"), sep ="|")
    target <- as.list( structure(c("green", "blue", "blue", "blue", "blue", "red"),
                 names = lbl))
    checkEquals(target, current)

    current <- edgeData(g1, attr ="type")
    lbl <- paste(c("a", "b", "d", "b", "d", "d"), 
                c( "b", "c", "c", "x", "x", "y") , sep ="|")
    target <- as.list( structure(c("low", "high", "unknown", "unknown", "unknown", "unknown"),
             names = lbl))
    checkEquals(target, current)
}

# in version prior to 1.41.1, the propagation of existing "user" edge attributes,
# that is, anything other than weight, failed if a new node is added which is
# lexically less than the any nodes already in an edge.
# the fix was simple:  see arguments to setBitCell in methods-graphBAM.R,
# the addNode method.  test that fix here

test_graphBAM_addNode_outOfAlphabeticalOrder_copyUserEdgeAttributes <- function(){

    from = c("a", "b", "m", "m")
    to   = c("b", "c", "y", "x")
    weight=c(2.2, 2.0, 0.4, 0.2)
    df <- data.frame(from, to, weight)
    g <- graphBAM(df, edgemode = "directed")
    
    edgeDataDefaults(g, attr="color") <- "blue"
    edgeDataDefaults(g, attr="type") <- "unknown"
    edgeData(g, from = c("m","a"), to = c("y", "b"), attr = "color") <- c("red", "green")
    edgeData(g, from = c("a", "b"), to = c("b", "c") , attr = "type") <- c("low", "high")

    expected.edge.names <-  c("a|b", "b|c", "m|x", "m|y")
    checkEquals(sort(names(edgeData(g, attr="color"))),
                expected.edge.names)

    checkEquals(unlist(edgeData(g, attr="color")[expected.edge.names],
                       use.names=FALSE),
                c("green", "blue", "blue", "red"))
                      
    g2 <- addNode("f", g)
   
       # make sure that the addition of node f does not disrupt
       # edgeData retrieval
    checkEquals(sort(names(edgeData(g2, attr="color"))),
                expected.edge.names)
    checkEquals(unlist(edgeData(g, attr="color")[expected.edge.names],
                       use.names=FALSE),
                c("green", "blue", "blue", "red"))

}



test_graphBAM_nodeUnion_Attributes <- function(use.factors=TRUE){
   
    setClass("myType", representation = representation(typ ="character")) 
    myType <- function(typ){ new("myType", typ = typ)}
    testFun <- function(x,y) {
        if(is(x, "myType")  && is(y, "myType")){
    
            if(x@typ =="aa" || y@typ == "ac")
                return("ax")
            else
                return("ab")
        } else return(as.character(NA))
    }
    funList <- structure(list(testFun), names ="gene")
    ft1 <- data.frame(from=c("a", "a", "a", "b", "b"),
            to  =c("b", "c", "d", "a", "d"),
            weight=c(1, 3.1, 5.4, 1, 2.2),
            stringsAsFactors = use.factors)

    g1 <- graphBAM(ft1, edgemode="directed")
    nodeDataDefaults(g1, attr="color") <- "cyan"
    nodeDataDefaults(g1, attr="type") <- "missing"
    nodeDataDefaults(g1, attr="kp") <- "missing"
    nodeDataDefaults(g1, attr="gene") <- "unknown"
    nodeData(g1, n = c("a", "b", "c") , attr = "color") <- c("red", "green", "blue")
    nodeData(g1, n = c("a", "b"), attr = "type") <- c("low", "high")
    nodeData(g1, n = c("a", "b"), attr = "kp") <- c("kplow", "kphigh")
    nodeData(g1, n = c("a", "b"), attr = "gene") <- c(myType("aa"), myType("bt"))

    ft1 <- data.frame(from=c("a", "a", "b"),
            to=c("b", "x", "z"),
            weight=c(6, 5, 2),
            stringsAsFactors = use.factors)
    g2 <- graphBAM(ft1,nodes = c("a","b", "c", "d", "x", "y", "z"), edgemode = "directed")
    
    nodeDataDefaults(g2, attr ="color") <- "cyan"
    nodeDataDefaults(g2, attr="type") <- "missing"
    nodeDataDefaults(g2, attr="gene") <- "unknown"
    nodeData(g2, n = c("a", "b", "x", "y", "z") , attr = "color") <- c("red", "red", "green", "pink", "yellow")
    nodeData(g2, n = c("a", "b"), attr = "type") <- c("low", "high")
    nodeData(g2, n = c("a", "b"), attr = "gene") <- c(myType("at"), myType("kt"))

    res <- graphUnion(g1, g2, nodeFun = funList)

    current <- nodeData(res, attr = "color")
    cn <- as.character(NA)
    target <- as.list( structure(c("red", cn, cn, "cyan", "green", "pink", "yellow"), 
                    names = c("a", "b", "c", "d", "x", "y", "z")))
    checkEquals(target, current)

    current <- nodeData(res, attr = "type")
    target <- as.list( structure(c("low", "high", "missing", "missing", "missing", "missing", "missing"), 
                    names = c("a", "b", "c", "d", "x", "y", "z")))
    checkEquals(target, current)

    current <- nodeData(res, attr = "kp")
    target <- as.list( structure(c("kplow", "kphigh", "missing", "missing", "missing",
                "missing", "missing"), 
                    names = c("a", "b", "c", "d", "x", "y", "z")))
    checkEquals(target, current)
  
    current <- nodeData(res, n = c("a", "b", "c", "d"), attr ="gene")
    target <- as.list( structure(c("ax", "ab", cn ,cn), names = c("a", "b", "c", "d")))
    checkEquals(target, current)

    current <- nodeData(res, n= c( "x", "y", "z"), attr ="gene")
    target <- as.list( structure(c("unknown","unknown", "unknown"), 
                    names = c("x", "y", "z")))
    checkEquals(target, current)
}


test_graphBAM_removeNode <- function(){

    from = c("a", "b", "d", "d")
    to   = c("b", "c", "y", "x")
    weight=c(2.2, 2.0, 0.4, 0.2)
    df <- data.frame(from, to, weight)
    g <- graphBAM(df, edgemode = "directed")
    nodeDataDefaults(g, attr="name") <- "NN"
    nodeData(g, n = c("a","b", "c", "d", "x", "y"), attr = "name") <-  
             c("a", "b", "c", "d", "x", "y")
    edgeDataDefaults(g, attr="name") <- "EE"
    edgeData(g, from = from, to = to , attr = "name") <-  paste(from, to , sep= "")

    res <- removeNode(c("x","b"), g)
    current <- nodeData(res, attr = "name")
    target <- as.list(structure( c("a", "c", "d", "y"), names =  c("a", "c",
                            "d", "y")))
    checkEquals(target, current)
     
    current <- edgeData(res, attr = "name")
    target <-  as.list(structure( "dy", names =  paste("d", "y", sep = "|")))
    checkEquals(current, target)                       

    res <- removeNode(c("x", "a"), g)
    current <- edgeData(res, attr = "name")
    target <-  as.list(structure( c("bc", "dy"), names =  paste(c("b", "d"),
                            c("c","y"), sep = "|")))
    checkEquals(target, current)
 }

test_edgeDataUndirectedGraph <- function() {

     df <- data.frame(from=c("a", "a", "c"),
                     to=c("b", "c", "d"),
                     weight=rep(1, 3), stringsAsFactors=FALSE)
     g <- graphBAM(df, edgemode="undirected")
     edgeDataDefaults(g, attr="EDA") <- 0
     edgeData(g, from="a", to="b", attr="EDA") <- 1
     edgeData(g, from="a", to="c", attr="EDA") <- 2
     edgeData(g, attr="EDA", from="a")

         # for edges where "a" is the source node, and to unspecified
     checkEquals(edgeData(g, attr="EDA", from="a")[["a|b"]], 1)
     checkEquals(edgeData(g, attr="EDA", from="a")[["a|c"]], 2)

         # specify single values for from and to
     checkEquals(edgeData(g, attr="EDA", from="a", to="b")[[1]], 1)
     checkEquals(edgeData(g, attr="EDA", from="a", to="c")[[1]], 2)

         # multiple target nodes
    x <- edgeData(g, from="a", to=c("b","c"), attr="EDA")
    checkEquals(length(x), 2)
    checkEquals(sort(names(x)), c("a|b", "a|c"))
    checkEquals(as.numeric(edgeData(g, from="a", to=c("b","c"), attr="EDA")),
                c(1,2))

    checkException(edgeData(g, from="a", to="bogus", attr="EDA"), silent=TRUE)
    checkException(edgeData(g, from=c("a", "c"),
                            to=c("bogus", "bagus"), attr="EDA"), silent=TRUE)
}
    
test_edgeMatrix <- function() {

    g <- graphBAM(data.frame(from="1", to="2", weight=1))
    mtx <- edgeMatrix(g, duplicates=FALSE)
    checkEquals(dim(mtx), c(2,1))
    checkEquals(rownames(mtx), c("from", "to"))
    checkEquals(as.numeric(mtx), c(1, 2))

    mtx.dup <- edgeMatrix(g, duplicates=TRUE)
    checkEquals(dim(mtx.dup), c(2,2))
    checkEquals(rownames(mtx.dup), c("from", "to"))
    checkEquals(as.numeric(mtx.dup), c(1, 2, 2, 1))
}

  
test_removeEdge_from_undirectedGraph <- function() {

  g <- graphBAM(data.frame(from="A", to="B", weight=1))
  g <- removeEdge(from="A", to="B", g=g)
  checkEquals(numEdges(g), 0)
  
  g <- graphBAM(data.frame(from="A", to="B", weight=1))
  g <- removeEdge(from="B", to="A", g=g)
  checkEquals(numEdges(g), 0)
}


test_AM2BAM <- function(){

     # test the fix for a nov 2014 bug, in which MultiGraph::.makeMDEdgeSet
     # fails to make edge_sets in the graphBAM constructor when only 1 edge
     # exists.   fix is at line 66 in MultiGraph.R: "drop=FALSE" added to the subset operation
   mtx <- matrix(c(0,1,0,0), ncol=2, byrow=TRUE, dimnames=list(c("A", "B"), c("A", "B")))

      # first create and check a simple (non-binary) adjacency matrix graph
   g.am <- graphAM(mtx, edgemode="directed")
   checkEquals(nodes(g.am),  c("A", "B"))
   checkEquals(edgemode(g.am), "directed")
   checkEquals(edgeNames(g.am), "A~B")

     # now convert to BAM
   g.bam <- as(g.am, "graphBAM")
   checkEquals(nodes(g.bam),  c("A", "B"))
   checkEquals(edgemode(g.bam), "directed")
   checkEquals(edgeNames(g.bam), "A~B")

}

test_isAdjacent <- function()
{
  am <- adjacencyMatrix   # for shorthand

  g <- graphBAM(data.frame(from="B", to="C", weight=1), edgemode="undirected")
  checkEquals(rownames(am(g)), c("B", "C"))
  checkEquals(colnames(am(g)), c("B", "C"))
  checkEquals(am(g)["B","C"], 1)
  checkEquals(am(g)["C","B"], 1)
  checkTrue(isAdjacent(g, "B", "C"))
  checkTrue(isAdjacent(g, "C", "B"))

  checkEquals(as.numeric(edgeMatrix(g)), c(1,2))   # reciprocal edges not stored

    # add a node, then an edge to the undirected graph g
  g <- addNode("A", g)
  checkEquals(nodes(g), c("A", "B", "C"))
    # just one edge
  checkEquals(sum(am(g)), 2)
  checkEquals(am(g)["B", "C"], 1)
  checkEquals(am(g)["C", "B"], 1)

  checkTrue(isAdjacent(g, "B", "C"))
  checkTrue(isAdjacent(g, "C", "B"))

  g <- addEdge(from="C", to="A", graph=g)
  checkEquals(sum(am(g)), 4)
  checkEquals(am(g)["B", "C"], 1)
  checkEquals(am(g)["C", "B"], 1)
  checkEquals(am(g)["A", "C"], 1)
  checkEquals(am(g)["C", "A"], 1)

    # robert's bug:  both of these fail though direct inspection
    # of either edgeMatrix or adjacencyMatrix show correct edges
  checkTrue(isAdjacent(g, "A", "C"))
  checkTrue(isAdjacent(g, "C", "A"))

     # now verify non-reciprocity of B-C edge in a directed graph
  gd <- graphBAM(data.frame(from="B", to="C", weight=1), edgemode="directed")
  checkEquals(rownames(am(gd)), c("B", "C"))
  checkEquals(colnames(am(gd)), c("B", "C"))
  checkEquals(am(gd)["B","C"], 1)
  checkTrue(isAdjacent(gd, "B", "C"))
  checkTrue(!isAdjacent(gd, "C", "B"))

    # add a node, then an edge to the directed graph gd
  gd <- addNode("A", gd)
  checkEquals(nodes(gd), c("A", "B", "C"))
    # just one edge
  checkEquals(sum(am(gd)), 1)
  checkEquals(am(gd)["B", "C"], 1)
  checkTrue(isAdjacent(gd, "B", "C"))

  gd <- addEdge(from="C", to="A", graph=gd)
  checkTrue(isAdjacent(gd, "C", "A"))

} # test_isAdjacent

# incomplete draft test supplied by Robert Castello (November 2014)
test_robertCastelos_addEdge_edgeData_bug <- function() {

  am <- adjacencyMatrix   # for shorthand
  #Sys.setlocale("LC_ALL", "C")
  #checkEquals(Sys.getlocale(), "C")

  g <- graphBAM(data.frame(from="B", to="C", weight=1))
  checkEquals(rownames(am(g)), c("B", "C"))
  checkEquals(colnames(am(g)), c("B", "C"))
  checkEquals(am(g)["B","C"], 1)

  edgeDataDefaults(g, "x") <- NA_real_
  g <- addNode("A", g)

  checkEquals(rownames(am(g)), c("A", "B", "C"))
  checkEquals(colnames(am(g)), c("A", "B", "C"))

  g <- addEdge(from="C", to="A", graph=g)
  checkEquals(am(g)["C", "A"], 1)

  ## this one works fine
  edgeData(g, from="A", to="C", "x") <- 10

  ## however, this one breaks the code:  no longer!
  edgeData(g, from="C", to="A", "x") <- 10
    # ensures that no error was found in the above operations
  checkTrue(TRUE)
}



#  Sys.setlocale("LC_ALL", "C")
#  checkEquals(Sys.getlocale(), "C")
#  g <- graphBAM(data.frame(from="B", to="C", weight=1))
#  edgeDataDefaults(g, "x") <- NA_real_
#  g <- addNode("A", g)
#  g <- addEdge(from="C", to="A", graph=g)
#
#    ## this one works fine
#  edgeData(g, from="A", to="C", "x") <- 10
#
#    ## however, this one breaks the code
#  edgeData(g, from="C", to="A", "x") <- 10


