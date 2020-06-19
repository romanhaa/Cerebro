### NOTE: List is an abstract type, so we just test with IntegerList

library(IRanges)

test_List_funprog <- function() {
  int1 <- c(1L,2L,3L,5L,2L,8L)
  int2 <- c(15L,45L,20L,1L,15L,100L)
  for (compress in c(TRUE, FALSE)) {
    collection <- IntegerList(int1, int2, int1, compress=compress)
    addcollect <- IntegerList(int2, int1, int1, compress=compress)
    checkIdentical(Reduce("+", collection), Reduce("+", list(int1, int2, int1)))
    checkIdentical(as.list(Filter(function(x) mean(x) > 10, collection)),
                   Filter(function(x) mean(x) > 10, list(int1, int2, int1)))
    checkIdentical(Find(function(x) mean(x) > 10, collection),
                   Find(function(x) mean(x) > 10, list(int1, int2, int1)))
    checkIdentical(Map("+", collection, addcollect),
                   Map("+", list(int1, int2, int1), list(int2, int1, int1)))
    checkIdentical(mapply("+", collection, addcollect),
                   mapply("+", list(int1, int2, int1), list(int2, int1, int1)))
    checkIdentical(Position(function(x) mean(x) > 10, collection),
                   Position(function(x) mean(x) > 10, list(int1, int2, int1)))
  }
}

