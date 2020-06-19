########################################
##EXAMPLES

set.seed(123)

x1 <- matrix(rnorm(500), nc=10)
y1 <- matrix(100+runif(500), nc=10)

eL1 <- new.env()
assign("exprs", x1, eL1)
assign("y1", y1, eL1)

eL2 <- list(exprs=x1, y1=y1)

var1 <- sample(1:3, 10, replace=TRUE)
var2 <- sample(c("M","F"), 10, replace=TRUE)

pD1 <- data.frame(var1, var2)
vL <- list(var1="Variable 1", var2="Variable 2")

pheno1 <- new("phenoData", pData=pD1, varLabels=vL)


##now how do we get the eList
v1<-new("eSet", eList=eL1, phenoData=pheno1)

v2<-new("eSet", eList=eL2,  phenoData=pheno1)

esApply(v1, 1, mean)

esApply(v2, 1, mean)

#c0 <- try( combine(v1,v2) ) # error
c1 <- combine(v2,v2)

v3 <- v2
pData(v3)$var3 <- 1:10
v3@phenoData@varLabels$var3 <- "var 3 label"
c2 <- combine(v2,v3)
