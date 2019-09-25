
  library(methods)
  
  .initpTreeClass(globalenv())

  pT1 <- new("pTree")

  npT<-pTreeInsert(pT1, "a", 10)

  objs <- letters[1:12]
  vals<-c(2, 14, 7, 21, 16, 18, 11, 1, 25, 3, 8, 5)

  for( i in 1:12)
    npT <- pTreeInsert(npT, objs[i], vals[i])


 unlist(npT@values)

  xx<-pTreeDelete(npT, "d")
