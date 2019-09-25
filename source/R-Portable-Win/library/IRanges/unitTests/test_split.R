test_splitAsList <- function() {
    ir <- IRanges(sample(100),sample(100)+100)
    ir2 <- unlist(split(ir, ceiling(1:100 / 10)))
    checkTrue(all(ir==ir2))
}
