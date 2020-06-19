test_unsplit <- function() {
  ir <- IRanges(1:5, 11:15)
  f <- factor(c("a", "b", "a", "b", "b"), c("b", "a", "c"))

  rl <- split(ir, f)
  checkIdentical(unsplit(rl, f), ir)  

  rl <- split(ir, f, drop=TRUE)
  checkIdentical(unsplit(rl, Rle(f), drop=TRUE), ir)
  checkException(unsplit(rl, f, drop=FALSE), silent=TRUE)

  v <- 1:5
  l <- splitAsList(v, f)
  checkIdentical(unsplit(l, Rle(f)), v)    
}

