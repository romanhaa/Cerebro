## ------------------------------------------------------------
## (wh) 05 Feb 2005
## Test and benchmark the three different implementations of
## graph intersection
## Results: for sparse graphs (e.g. nodes=edges=2000), intersection3
##   is fastest; for dense graphs (e.g. nodes=200, edges=10000),
##   intersection2 is faster. With the given parameters, I obtained:
##
## Sparse
## t1  27.74 0.23 58.75    0    0
## t2  27.35 0.11 61.68    0    0
## t3   5.03 0.02 10.98    0    0
##
## Dense:
## t1  2.61 0.00 2.77    0    0
## t2  1.13 0.01 1.57    0    0
## t3  6.28 0.01 7.15    0    0

library("graph")
options(error=recover)

nodes = 2000; edges = 2000  ## sparse
## nodes = 200; edges = 10000  ## dense

V = paste(formatC(1:nodes, width=5, flag="0"))
B = 5

set.seed(123)
g1 <-lapply(1:B, function(i) randomEGraph(V=V, edges=edges))
g2 <-lapply(1:B, function(i) randomEGraph(V=V, edges=edges))

t3 <- system.time(
  i3 <- mapply(intersection3, g1, g2)
)

t1 <- system.time(
  i1 <- mapply(intersection, g1, g2)
)

t2 <- system.time(
  i2 <- mapply(intersection2, g1, g2)
)

identical.graphs = function(g1, g2) {
  if(!identical(nodes(g1), nodes(g2)))
    stop("Baeh 1")
  e1 <- edges(g1)
  e2 <- edges(g2)
  s = mapply(function(x,y) all(sort(x)==sort(y)), e1, e2)
  if(!all(s))
    stop("Baeh 2")
  return(TRUE)       
}

cat("system.time:\n")
print(rbind(t1,t2,t3))

## Check whether all are identical
cat("Now checking:\n")
for(i in seq(along=i1)) {
  stopifnot(identical.graphs(i1[[i]], i2[[i]]))
  stopifnot(identical.graphs(i1[[i]], i3[[i]]))
}

