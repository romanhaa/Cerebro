## Basic Graph Tests

## Each test has:
##     - an input graph specified as an edge list
##     - an operation name (or should it be the function?
##     - an expected result.  If the result is a graph object, it will
##       be expressed as an edge list.
##
## Do we need a way to indicate nodes with no edges?
##

makeFT <- function(from, to) {
    ord <- order(from)
    from <- from[ord]
    to <- to[ord]
    cbind(from, to)
}

basicDirected <- function() {
    from <- c("b", "b", "b", "i", "o", "c", "c", "e")
    to   <- c("e", "i", "s", "o", "c", "i", "o", "c")
    w    <- seq_len(length(to))
    list(nodes = letters[1:20],
         edges = makeFT(from, to),
         weights = w,
         edgemode = "directed")
}

basicUndirected <- function() {
    from <- c("b", "b", "b", "i", "c", "e")
    to   <- c("e", "i", "s", "o", "o", "c")
    w    <- seq_len(length(to))
    list(nodes = letters[1:20],
         edges = makeFT(from, to),
         weights = w,
         edgemode = "undirected")
}

basic_to_ft <- function(g) {
    from <- match(g$edges[, 1L], g$nodes)
    to <- match(g$edges[, 2L], g$nodes)
    ft <- cbind(from, to)
    ft
}

create <-
    list(
         "graphAM" = function(gDesc) {
             ft <- basic_to_ft(gDesc)
             numNodes <- length(gDesc$nodes)
             mat <- matrix(0L, nrow = numNodes, ncol = numNodes,
                           dimnames = list(NULL, gDesc$nodes))
             coord <- graph:::coordToIndex(ft[, 1L], ft[, 2L], numNodes)
             w <- gDesc$weights
             if (gDesc$edgemode == "undirected") {
                 coord <- c(coord,
                            graph:::coordToIndex(ft[, 2L], ft[, 1L], numNodes))
                 w <- c(w, w)
             }
             mat[coord] <- w
             graphAM(adjMat = mat, edgemode = gDesc$edgemode,
                 values = list(weight = 1L))
         },
         "graphNEL" = function(gDesc) {
             edgeL <- split(gDesc$edges[ , 2L], gDesc$edges[ , 1L])
             if (gDesc$edgemode == "undirected") {
                 f <- gDesc$edges[, 1L]
                 t <- gDesc$edges[, 2L]
                 ft <- c(f, t)
                 tf <- c(t, f)
                 edgeL <- split(tf, ft)
             }
             g <- graphNEL(nodes = gDesc$nodes,
                      edgeL = edgeL, edgemode = gDesc$edgemode)
             edgeDataDefaults(g, attr = "weight") <- 1L
             edgeData(g, from = gDesc$edges[, 1L],
                      to = gDesc$edges[, 2L],
                      attr = "weight") <- gDesc$weights
             g

         })

graph2desc <- function(g) {
    nms <- nodes(g)
    ft <- t(edgeMatrix(g))
    from <- nms[ft[, 1L]]
    to <- nms[ft[, 2L]]
    list(nodes = nms,
         edges = makeFT(from, to),
         weights = unlist(edgeWeights(g), use.names = FALSE),
         edgemode = edgemode(g))
}

toGraphDesc <- list("graphAM" = graph2desc,
                    "graphNEL" = graph2desc)

gam <- create$graphAM(basicDirected())
gnel <- create$graphNEL(basicDirected())

gam0 <- toGraphDesc$graphAM(gam)
gnel0 <- toGraphDesc$graphNEL(gnel)

ugam <- create$graphAM(basicUndirected())
ugnel <- create$graphNEL(basicUndirected())

