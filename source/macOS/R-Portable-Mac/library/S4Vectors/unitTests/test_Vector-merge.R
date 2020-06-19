test_Vector_merge <- function() {
    library(GenomicRanges)

    ## Binary merge

    gr <- GRanges(c("chr1:1-1000", "chr1:2000-3000"), a=1:2, b=2:1)
    gr2 <- GRanges(c("chr1:1-1000", "chr1:2000-3000"), c=c(1,3), d=c(3,1))

    target <- granges(gr)
    mcols(target) <- DataFrame(a=1:2, b=2:1, c=c(1,3), d=c(3,1))
    current <- merge(gr, gr2)
    checkIdentical(target, current)
        
    gr <- GRanges(c("chr1:1-1000", "chr1:2000-3000", "chr1:1-10"),
                  a=1:3, b=c(2,1,3))

    target <- granges(gr2)
    mcols(target) <- DataFrame(a=1:2, b=c(2,1), c=c(1,3), d=c(3,1))
    current <- merge(gr, gr2)
    checkIdentical(target, current)
    current <- merge(gr, gr2, all.y=TRUE)
    checkIdentical(target, current)

    target <- granges(gr)
    mcols(target) <- DataFrame(a=1:3, b=c(2,1,3), c=c(1,3,NA), d=c(3,1,NA))
    current <- merge(gr, gr2, all=TRUE, sort=FALSE)
    checkIdentical(target, current)
    current <- merge(gr, gr2, all.x=TRUE, sort=FALSE)
    checkIdentical(target, current)
    target <- sort(target)
    current <- merge(gr, gr2, all=TRUE, sort=TRUE)
    checkIdentical(target, current)

    x <- GRanges(c("chr1:1-1000", "chr2:2000-3000"),
                 score=c(0.45, 0.1), a1=c(5L, 7L), a2=c(6, 8))
    y <- GRanges(c("chr2:150-151", "chr1:1-10",
                   "chr2:2000-3000", "chr2:2000-3000"),
                 score=c(0.7, 0.82, 0.1, 0.2),
                 b1=c(0L, 5L, 1L, 7L), b2=c(1, -2, 1, 1.5))

    checkException(merge(x, y[-3]))

    target0 <- c(granges(x), granges(y[-4]))[c(4, 1, 3, 2)]
    mcols(target0) <- DataFrame(score=c(0.82, 0.45, 0.7, 0.1),
                                a1=c(NA, 5L, NA, 7L),
                                a2=c(NA, 6, NA,  8),
                                b1=c(5L, NA, 0L, 1L),
                                b2=c(-2, NA, 1, 1))

    current <- merge(x, y[-4], all=TRUE)
    checkIdentical(target0, current)

    current <- merge(x, y[-4], all.x=TRUE)
    checkIdentical(target0[c(2, 4)], current)

    current <- merge(x, y[-4], all.y=TRUE)
    target <- target0[c(3, 4, 1)]
    seqlevels(target) <- seqlevels(current)
    checkIdentical(target, current)

    current <- merge(x, y[-4])
    checkIdentical(target0[4], current)

    ## Self merge is a no-op if 'sort=FALSE' (or object already sorted) and
    ## if the object has no duplicates

    checkIdentical(x, merge(x, x))

    ## N-ary merge

    current <- merge(x, y[-4], x, all=TRUE)
    checkIdentical(target0, current)
}

