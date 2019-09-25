##some new things

library(methods)
library(graph)

library(GO)
library(hgu95a)

 xx <- ls(env = hgu95aGO)
set.seed(1234)
 myGenes <- sample(xx, 100)
 mG <- mget(myGenes, env=hgu95aGO)

makeGoGraph <- function(x) {
   library(GO)
   newNodes <- get(x, env=hgu95aGO)
   if( is.na(x) )
     return(NULL)
   oldEdges <- vector("list", length=0)
   oldNodes <- vector("character", length=0)
   done <- FALSE
   while( !done ) {
     newNodes <- newNodes[!(newNodes %in% oldNodes)]
     if( length(newNodes) == 0 )
       done <- TRUE
     else {
         oldNodes <- c(oldNodes, newNodes)
         numE <- length(newNodes)
         nedges <- vector("list", length=numE)
         names(nedges) <- newNodes
         nedges <- mget(newNodes, env=GOmolecularfunction)
         nedges <- nedges[!is.na(nedges)]
         oldEdges <- c(oldEdges, nedges)
         newNodes <- sort(unique(unlist(nedges)))
     }
   }
   rE <- vector("list", length=length(oldNodes))
   names(rE) <- oldNodes
   rE[names(oldEdges)] <- oldEdges
   return(list(nodes=oldNodes, edges=rE))
}

 Gmf1 <- makeGoGraph(myGenes[1])

##old examples
 library(graph)
 data(pmedu95aAffy)
 pmG <- pmedu95aAffy

 edgeL <- lapply(pmG@edges, function(x) list(edges=x))
 pmG <- graphNEL(nodes=pmG@nodes, edgeL=edgeL)

 xx1<-acc(pmG, "1025_g_at")

 xx2<-acc2(pmG, "1025_g_at")

 set.seed(12345)

 myNodes <- sample(nodes(pmG), 500)

 pmS <- makeSubGraph(pmG, myNodes)

 xx<- acc(pmS, "35990_at")

 xy <- dfs(pmS)

 zz <- acc(pmS, "35990_at")

 zz2 <- acc(pmS, "35990_at")

 pm2 <- makeSubGraph(pmG, myNodes[1:100])

  pm3<-isect(pm2, pmS)


 ##Sept 27 -- trying to test some graph code

 #library(methods)

 #setwd("c:\\cygwin/home/rgentlem\\Software\\graph\\R")

 #source("graph.R")

 #.initGraph(globalenv())

 x<-1:100

 rw <- rep("a", 100)
 for(i in 1:100) rw[i] <- paste(sample(letters, 10, replace=TRUE),
             sep="", collapse="")


 set.seed(121)

 ##this is a directed graph -- the nodes are one way nodes
 y<- vector("list", length=100)
 for(i in 1:100) {
  nnodes<- floor(runif(1)*20)
  y[[i]]  <- list(edges=sample(x, nnodes), weights=runif(nnodes))
 }

 sapply(y, function(x) length(x$weights))

 names(y) <- rw
 g1 <- graphNEL(nodes=rw, edgeL=y)

 z<-lapply(y, function(x) {x$weights<-NULL; x})
 g2 <- graphNEL(nodes=rw, edgeL=z)

 #get the adjcency list for node number 10
 vv<-adj(g1, 10)

 #get the accessibility list for node number 10
 vw <- acc(g1, 2)

 ##for an undirected graph we generate a node list for each node

 set.seed(333)

 y2<- vector("list", length=100)
 for(i in 1:100) y2[[i]] <- list(edges=numeric(0), weights=numeric(0))
 for(i in 1:100) {
   nnodes<- floor(runif(1)*3)
   jj<-sample(x, nnodes)
   for (j in jj) {
     wt <- 18*runif(1)
     y2[[i]]$edges <- c(y2[[i]]$edges, j)
      y2[[i]]$weights <- c(y2[[i]]$weights, wt)
     y2[[j]]$edges <- c(y2[[j]]$edges, i)
     y2[[j]]$weights <- c(y2[[j]]$weights, wt)
   }
 }

 g3 <- graphNEL(nodes=rw, edgeL=y2)

 b1 <- isect(g1, g3)

 sN1 <- sample(1:100, 20)

 g4 <- subGraph(g1, sN1)

 E1 <- edgeL(g3)
 E2 <- edgeL(g3, sN1)
