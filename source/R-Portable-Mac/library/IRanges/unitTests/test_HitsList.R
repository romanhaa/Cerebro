test_HitsList_as_matrix <- function() {
  x <- IRangesList(chr1=IRanges(1, 5), chr2=IRanges(6, 10))
  y <- IRangesList(chr2=IRanges(8, 10))
  checkIdentical(as.matrix(findOverlaps(x, y)),
                 cbind(queryHits = 2L, subjectHits = 1L))
}
