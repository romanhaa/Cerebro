library(IRanges)  # for CharacterList

test_expand <- function(){
  ## setup
  aa <- CharacterList("a", paste0("d", 1:2), paste0("b", 1:3), c(), "c")
  bb <- CharacterList(paste0("sna", 1:2),"foo", paste0("bar", 1:3),c(),"hica")
  df <- DataFrame(aa=aa, bb=bb, cc=11:15)

  ## tests
  ## test one col without dropping
  res1 <- expand(df, colnames="aa", keepEmptyRows=TRUE)
  checkTrue(dim(res1)[1]==8)
  checkTrue(dim(res1)[2]==3)
  checkIdentical(res1$aa,c("a","d1","d2","b1","b2","b3",NA,"c"))
  checkIdentical(res1$bb[[4]],c("bar1","bar2","bar3"))
  ## test one col with dropping
  res2 <- expand(df, colnames="aa", keepEmptyRows=FALSE)
  checkTrue(dim(res2)[1]==7)
  checkTrue(dim(res2)[2]==3)
  checkIdentical(res2$aa,c("a","d1","d2","b1","b2","b3","c"))
  checkIdentical(res2$bb[[4]],c("bar1","bar2","bar3"))
  ## test two columns no dropping
  res3 <- expand(df, colnames=c("aa","bb"),  keepEmptyRows=TRUE)
  checkTrue(dim(res3)[1]==15)
  checkTrue(dim(res3)[2]==3)
  checkIdentical(res3$aa,
    c("a","a","d1","d2","b1","b1","b1","b2","b2","b2","b3","b3","b3",NA,"c"))
  checkIdentical(as.character(as.data.frame(res3[14,])),
                 c(NA, NA, "14"))
  ## test two columns with dropping
  res4 <- expand(df, colnames=c("aa","bb"),  keepEmptyRows=FALSE)
  checkTrue(dim(res4)[1]==14)
  checkTrue(dim(res4)[2]==3)
  checkIdentical(res4$aa,
    c("a","a","d1","d2","b1","b1","b1","b2","b2","b2","b3","b3","b3","c"))
  ## inverted order (different sorting of 2 cols, no dropping
  res5 <- expand(df, colnames=c("bb","aa"),  keepEmptyRows=TRUE)
  checkTrue(dim(res5)[1]==15)
  checkTrue(dim(res5)[2]==3)
  checkIdentical(res5$aa,
    c("a","a","d1","d2","b1","b2","b3","b1","b2","b3","b1","b2","b3",NA,"c"))
  ## inverted order (different sorting of 2 cols, with dropping  
  res6 <- expand(df, colnames=c("bb","aa"),  keepEmptyRows=FALSE)
  checkTrue(dim(res6)[1]==14)
  checkTrue(dim(res6)[2]==3)
  checkIdentical(res6$aa,
    c("a","a","d1","d2","b1","b2","b3","b1","b2","b3","b1","b2","b3","c"))
}
