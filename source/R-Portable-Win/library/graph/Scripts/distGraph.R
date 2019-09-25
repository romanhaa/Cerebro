
 ##some code for distGraphs

 ## a concrete example:
 library(Biobase)

  library(mva)

 setwd("c:/cygwin/home/rgentlem/Software/graph/R")

 source("clustergraph.R")

  data(eset)
  d1 <- dist(exprs(eset))
  length(d1)
  ##should be 124750
  ##500*499/2
 ##  124750

 ## we could take deciles and 

   deciles <- quantile(unclass(d1), probs=seq(0.1,0.9,0.1))

   dG <- new("distGraph", Dist=d1)

   dG2 <- threshold(dG, deciles[8])
 
   dG3 <- threshold(dG, deciles[4])

   xx<- d1[1:10, c(3,5,11)]

   ad1 <- adj(dG3, 10)

   ac1 <- acc(dG3, 10)

  cc <- connComp(dG3)

   dG4 <- threshold(dG, deciles[1]/2.1)

   cc4 <- connComp(dG4)

  J<- sapply(cc4, function(x) length(x))
  table(J)

