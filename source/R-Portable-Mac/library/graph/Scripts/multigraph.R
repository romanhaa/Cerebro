
  
  set.seed(123)
  V <- LETTERS[1:4]
  edL <- vector("list", length=4)
  names(edL) <- V
  for(i in 1:4)
    edL[[i]] <- list(edges=5-i, weights=runif(1))

  e1 = new("edgeSetNEL", edgemode="undirected",
    edgeL = edL)

  x = matrix(rnorm(12), nrow=4, dimnames=list(V, NULL))
  d1 = as.matrix(dist(x))
 
  e2 = new("edgeSetAM", edgemode="undirected", 
           adjMat = ifelse(d1>1, 1, 0))


  mg1 = new("multiGraph", nodes = V, edgeL = list(e1, e2))

  mg1

  edges(mg1)

  isDirected(mg1)

  numEdges(mg1)

