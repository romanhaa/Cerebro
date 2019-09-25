### R code from vignette source 'MultiGraphClass.Rnw'

###################################################
### code chunk number 1: loadGraph
###################################################
library(graph)


###################################################
### code chunk number 2: creategraphBAM
###################################################
df <- data.frame(from = c("SEA", "SFO", "SEA", "LAX", "SEA"),
                   to = c("SFO", "LAX", "LAX", "SEA", "DEN"),
               weight = c( 90, 96, 124, 115, 259))
g <- graphBAM(df, edgemode = "directed")
g 


###################################################
### code chunk number 3: nodeAndWeights
###################################################
nodes(g)
edgeWeights(g, index = c("SEA", "LAX"))


###################################################
### code chunk number 4: addNodeEdge
###################################################
g <- addNode("IAH", g)
g <- addEdge(from = "DEN", to = "IAH", graph = g, weight = 120)
g


###################################################
### code chunk number 5: removeEdge
###################################################
g <- removeEdge(from ="DEN", to = "IAH", g)
g <- removeNode(node = "IAH", g)
g


###################################################
### code chunk number 6: subGraph
###################################################
g <- subGraph(snodes = c("DEN","LAX", "SEA"), g)
g


###################################################
### code chunk number 7: fromTo
###################################################
extractFromTo(g)


###################################################
### code chunk number 8: loadData1
###################################################
data("esetsFemale")
data("esetsMale")


###################################################
### code chunk number 9: dataFrames
###################################################
dfMale <- esetsMale[["brain"]]
dfFemale <- esetsFemale[["brain"]]
head(dfMale)


###################################################
### code chunk number 10: creategraphBAMs
###################################################
male <- graphBAM(dfMale, edgemode = "directed")
female <- graphBAM(dfFemale, edgemode = "directed")


###################################################
### code chunk number 11: bamIntersect
###################################################
intrsct <- graphIntersect(male, female, edgeFun=list(weight =  sum))
intrsct


###################################################
### code chunk number 12: removeEdges
###################################################
resWt <- removeEdgesByWeight(intrsct, lessThan = 1.5)


###################################################
### code chunk number 13: updateColor
###################################################
ftSub <- extractFromTo(resWt)
edgeDataDefaults(male, attr = "color") <- "white"
edgeDataDefaults(female, attr = "color") <- "white"
edgeData(male, from = as.character(ftSub[,"from"]), 
        to = as.character(ftSub[,"to"]), attr = "color") <- "red"

edgeData(female, from = as.character(ftSub[,"from"]), 
        to = as.character(ftSub[,"to"]), attr = "color") <- "red"



###################################################
### code chunk number 14: loadRBGL
###################################################
library(graph)
library(RBGL)


###################################################
### code chunk number 15: createDataFrames
###################################################
ft1 <- data.frame(
        from = c("SEA", "SFO", "SEA", "LAX", "SEA"),
          to = c("SFO", "LAX", "LAX", "SEA", "DEN"),
      weight = c( 90, 96, 124, 115, 259))

ft2 <- data.frame(
        from = c("SEA", "SFO", "SEA", "LAX", "SEA", "DEN", "SEA", "IAH", "DEN"),
          to = c("SFO", "LAX", "LAX", "SEA", "DEN", "IAH", "IAH", "DEN", "BWI"),
       weight= c(169, 65, 110, 110, 269, 256, 304, 256, 271))

ft3 <- data.frame( 
   from = c("SEA", "SFO", "SEA", "LAX", "SEA", "DEN", "SEA", "IAH", "DEN", "BWI"),
    to  = c("SFO", "LAX", "LAX", "SEA", "DEN", "IAH", "IAH", "DEN", "BWI", "SFO"),
 weight = c(237, 65, 156, 139, 281, 161, 282, 265, 298, 244))

ft4 <- data.frame( 
    from = c("SEA", "SFO", "SEA", "SEA", "DEN", "SEA", "BWI"),
     to  = c("SFO", "LAX", "LAX", "DEN", "IAH", "IAH", "SFO"),
  weight = c(237, 60, 125, 259, 265, 349, 191))


###################################################
### code chunk number 16: createMG
###################################################
esets <- list(Alaska = ft1, United = ft2, Delta = ft3, American = ft4)
mg <- MultiGraph(esets, directed = TRUE)
mg


###################################################
### code chunk number 17: cities
###################################################
nodes(mg)


###################################################
### code chunk number 18: DeltafromSeattle
###################################################
mgEdgeData(mg, "Delta", from = "SEA", attr = "weight")


###################################################
### code chunk number 19: nodeData
###################################################
nodeDataDefaults(mg, attr="shape") <- "square"
nodeData(mg, n = c("SEA", "DEN", "IAH", "LAX", "SFO"), attr = "shape")  <- 
    c("triangle", "circle", "circle", "circle", "circle")


###################################################
### code chunk number 20: nodeDataVal
###################################################
nodeData(mg,  attr = "shape")  


###################################################
### code chunk number 21: edgeDataVal
###################################################
mgEdgeDataDefaults(mg, "Delta", attr = "color")  <- "red"
mgEdgeData(mg, "Delta", from = c("SEA", "SEA", "SEA", "SEA"), 
        to = c("DEN", "IAH", "LAX", "SFO"), attr = "color") <- "green"


###################################################
### code chunk number 22: mgEdgeDataVal
###################################################
mgEdgeData(mg, "Delta", attr = "color")


###################################################
### code chunk number 23: subsetMG
###################################################
g <- subsetEdgeSets(mg, edgeSets = c("Alaska", "United", "Delta"))


###################################################
### code chunk number 24: intersecmg
###################################################
edgeFun <- list( weight = min)
gInt <- edgeSetIntersect0(g, edgeFun = edgeFun)
gInt


###################################################
### code chunk number 25: intersectWeights
###################################################
mgEdgeData(gInt, "Alaska_United_Delta", attr= "weight")


###################################################
### code chunk number 26: loadData
###################################################
data("esetsFemale")
data("esetsMale")
names(esetsFemale)
head(esetsFemale$brain)


###################################################
### code chunk number 27: createMultiGraphs
###################################################
female  <- MultiGraph(edgeSets = esetsFemale, directed = TRUE)
male  <- MultiGraph(edgeSets = esetsMale, directed = TRUE )
male
female


###################################################
### code chunk number 28: graphBAMs
###################################################
maleBrain <- extractGraphBAM(male, "brain")[["brain"]]
maleBrain
femaleBrain <-  extractGraphBAM(female, "brain")[["brain"]]


###################################################
### code chunk number 29: edgeDistance
###################################################
maleWt <- bellman.ford.sp(maleBrain, start = c("10024416717"))$distance
maleWt <- maleWt[maleWt != Inf & maleWt != 0]
maleWt

femaleWt <- bellman.ford.sp(femaleBrain, start = c("10024416717"))$distance
femaleWt <- femaleWt[femaleWt != Inf & femaleWt != 0]
femaleWt


###################################################
### code chunk number 30: nodeAttr
###################################################
nodeDataDefaults(male, attr = "color") <- "gray"
nodeData(male , n = c("10024416717", names(maleWt)), attr = "color" ) <- c("red")

nodeDataDefaults(female, attr = "color") <- "gray"
nodeData(female , n = c("10024416717", names(femaleWt)), attr = "color" ) <- c("red")


###################################################
### code chunk number 31: nodeSub
###################################################
resInt <- graphIntersect(male, female)
resInt


