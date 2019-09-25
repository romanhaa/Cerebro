test_IRanges_coverage <- function() {
  ir <- IRanges(c(1, 8, 14, 15, 19, 34, 40), width = c(12, 6, 6, 15, 6, 2, 7))
  checkIdentical(as.vector(coverage(ir)),
                 rep(c(1L, 2L, 1L, 2L, 3L, 2L, 1L, 0L, 1L, 0L, 1L),
                     c(7, 5, 2, 4, 1, 5, 5, 4, 2, 4, 7)))
  ir <- IRanges(start=c(-2L, 6L, 9L, -4L, 1L, 0L, -6L, 10L),
                width=c( 5L, 0L, 6L,  1L, 4L, 3L,  2L,  3L))
  checkIdentical(as.vector(coverage(ir)),
                 rep(c(3L, 1L, 0L, 1L, 2L, 1L), c(2, 2, 4, 1, 3, 2)))
  checkIdentical(as.vector(coverage(ir, shift=7)),
                 rep(c(1L, 0L, 1L, 2L, 3L, 1L, 0L, 1L, 2L, 1L),
                     c(3, 1, 2, 1, 2, 2, 4, 1, 3, 2)))
  checkIdentical(as.vector(coverage(ir, shift=7, width=27)),
                 rep(c(1L, 0L, 1L, 2L, 3L, 1L, 0L, 1L, 2L, 1L, 0L),
                     c(3, 1, 2, 1, 2, 2, 4, 1, 3, 2, 6)))
}

