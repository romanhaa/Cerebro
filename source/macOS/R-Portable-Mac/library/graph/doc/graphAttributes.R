### R code from vignette source 'graphAttributes.Rnw'

###################################################
### code chunk number 1: exampleGraph1
###################################################
library("graph")
mat <- matrix(c(0, 0, 1, 1, 
                0, 0, 1, 1, 
                1, 1, 0, 1, 
                1, 1, 1, 0),
              byrow=TRUE, ncol=4)
rownames(mat) <- letters[1:4]
colnames(mat) <- letters[1:4]


###################################################
### code chunk number 2: exampleGraph2
###################################################
g1 <- graphAM(adjMat=mat)


###################################################
### code chunk number 3: foo
###################################################
if (require("Rgraphviz")) {
    gn = as(g1, "graphNEL")
    plot(gn, nodeAttrs=makeNodeAttrs(gn, shape="circle", fillcolor="orange"))
} else {
  plot(1, 1, main="Rgraphviz required for this plot")
}


###################################################
### code chunk number 4: edgeDataDefaults1
###################################################
edgeDataDefaults(g1)


###################################################
### code chunk number 5: edgeDataDefaults2
###################################################
edgeDataDefaults(g1, "weight") <- 1
edgeDataDefaults(g1, "code") <- "plain"
edgeDataDefaults(g1)


###################################################
### code chunk number 6: edgeDataDefaults3
###################################################
edgeDataDefaults(g1, "weight")


###################################################
### code chunk number 7: edgeData1
###################################################
edgeData(g1, from="a", to="d", attr="weight")
edgeData(g1, from="a", attr="weight")
edgeData(g1, to="a", attr="weight")
allAttrsAllEdges <- edgeData(g1)
weightAttrAllEdges <- edgeData(g1, attr="weight")


###################################################
### code chunk number 8: edgeData2
###################################################
edgeData(g1, from="a", to="d", attr="weight") <- 2
edgeData(g1, from="a", attr="code") <- "fancy"
edgeData(g1, from="a", attr="weight")
edgeData(g1, from="a", attr="code")


###################################################
### code chunk number 9: edgeData3
###################################################
f <- c("a", "b")
t <- c("c", "c")
edgeData(g1, from=f, to=t, attr="weight") <- 10
edgeData(g1, from=f, to=t, attr="weight")


###################################################
### code chunk number 10: edgeData4
###################################################
edgeData(g1, from=f, to=t, attr="weight") <- c(11, 22)
edgeData(g1, from=f, to=t, attr="weight")


###################################################
### code chunk number 11: edgeData5
###################################################
edgeData(g1, from="a", to="d", attr="code") <- list(1:10)
edgeData(g1, from=f, to=t, attr="weight") <- mapply(c, f, t, "e", SIMPLIFY=FALSE) 
edgeData(g1, from="a", to="d", attr="code")
edgeData(g1, from=f, to=t, attr="weight")


###################################################
### code chunk number 12: defaultNodeData1
###################################################
nodeDataDefaults(g1)
nodeDataDefaults(g1, attr="weight") <- 1
nodeDataDefaults(g1, attr="type") <- "vital"
nodeDataDefaults(g1)
nodeDataDefaults(g1, "weight")


###################################################
### code chunk number 13: nodeData1
###################################################
nodeData(g1, n="a")
nodeData(g1, n="a", attr="weight") <- 100
nodeData(g1, n=c("a", "b"), attr="weight")
nodeData(g1, n=c("a", "b"), attr="weight") <- 500
nodeData(g1, n=c("a", "b"), attr="weight")
nodeData(g1, n=c("a", "b"), attr="weight") <- c(11, 22)
nodeData(g1, n=c("a", "b"), attr="weight")


###################################################
### code chunk number 14: other
###################################################
## We need to reconcile this
#g2 <- as(g1, "graphNEL")
#edgeWeights(g2)



