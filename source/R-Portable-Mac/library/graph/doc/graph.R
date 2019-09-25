### R code from vignette source 'graph.Rnw'

###################################################
### code chunk number 1: g1
###################################################
library(graph)
set.seed(123)
g1 = randomEGraph(LETTERS[1:15], edges=100)
g1


###################################################
### code chunk number 2: simplefuns
###################################################
nodes(g1)
degree(g1)
adj(g1, "A")
acc(g1, c("E", "G"))


###################################################
### code chunk number 3: subG
###################################################

 sg1 = subGraph(c("A", "E", "F","L"), g1)
 boundary(sg1, g1)

 edges(sg1)
 edgeWeights(sg1)



###################################################
### code chunk number 4: example1
###################################################
V <- LETTERS[1:4]
edL1 <- vector("list", length=4)
names(edL1) <- V
for(i in 1:4)
  edL1[[i]] <- list(edges=c(2,1,4,3)[i], weights=sqrt(i))
gR <- graphNEL(nodes=V, edgeL=edL1)

edL2 <- vector("list", length=4)
names(edL2) <- V
for(i in 1:4)
  edL2[[i]] <- list(edges=c(2,1,2,1)[i], weights=sqrt(i))
gR2 <- graphNEL(nodes=V, edgeL=edL2, edgemode="directed")



###################################################
### code chunk number 5: addNodes
###################################################

gX = addNode(c("E", "F"), gR)
gX
gX2 = addEdge(c("E", "F", "F"), c("A", "D", "E"), gX, c(1,2,3))
gX2

gR3 = combineNodes(c("A","B"), gR, "W")
gR3

clearNode("A", gX)



###################################################
### code chunk number 6: combine
###################################################

##find the underlying graph
ugraph(gR2)



###################################################
### code chunk number 7: unions
###################################################
set.seed(123)
gR3 <- randomGraph(LETTERS[1:4], M<-1:2, p=.5)

 x1 <-  intersection(gR,gR3)
 x1

 x2 <-  union(gR,gR3)
 x2

 x3 <- complement(gR)

 x3



###################################################
### code chunk number 8: randomEGraph
###################################################
set.seed(333)
V = letters[1:12]
g1 = randomEGraph(V, .1)
g1
g2 = randomEGraph(V, edges=20)
g2



###################################################
### code chunk number 9: randomGraph
###################################################
set.seed(23)
V <- LETTERS[1:20]
M <- 1:4
g1 <- randomGraph(V, M, .2)


###################################################
### code chunk number 10: randomNodeGraph
###################################################
    set.seed(123)
     c1 <- c(1,1,2,4)
     names(c1) <- letters[1:4]

     g1 <- randomNodeGraph(c1)





###################################################
### code chunk number 11: rGraph
###################################################
 g1

 g1cc <- connComp(g1)
 g1cc

 g1.sub <- subGraph(g1cc[[2]], g1)
 g1.sub



###################################################
### code chunk number 12: dfs
###################################################
DFS(gX2, "E")



###################################################
### code chunk number 13: clusterGraph
###################################################

cG1 <- new("clusterGraph", clusters=list(a=c(1,2,3), b=c(4,5,6)))
cG1
acc(cG1, c("1", "2"))



###################################################
### code chunk number 14: distanceGraph
###################################################
set.seed(123)
x <- rnorm(26)
names(x) <- letters
library(stats)
d1 <- dist(x)
g1 <- new("distGraph", Dist=d1)
g1



