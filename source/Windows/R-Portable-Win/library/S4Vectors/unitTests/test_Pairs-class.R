test_Pairs <- function() {
    score <- rnorm(10)
    p <- Pairs(1:10, Rle(1L, 10), score=score, names=letters[1:10])
    checkIdentical(first(p), 1:10)
    checkIdentical(mcols(p)$score, score)
    checkIdentical(p %in% p[1:5], c(rep(TRUE, 5), rep(FALSE, 5)))
    checkIdentical(as.data.frame(p),
                   data.frame(first=first(p), second=second(p), score,
                              names=names(p), stringsAsFactors=FALSE))
    z <- zipup(p)
    first(p) <- Rle(1:10)
    checkIdentical(zipdown(z), p)
}
