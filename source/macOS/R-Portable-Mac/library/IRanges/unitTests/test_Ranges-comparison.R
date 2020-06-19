test_pcompare_IntegerRanges <- function()
{
    x1 <- IRanges(6:16, width=4)
    y <- IRanges(11, 14)
    target <- c(-6:-4, -4L, -4L, 0L, 4L, 4L, 4:6)
    checkIdentical(target, pcompare(x1, y))
    checkIdentical(-target, pcompare(y, x1))

    x2 <- IRanges(4:16, width=6)
    target <- c(-6:-4, -4L, -4L, -3L, -2L, 1L, 4L, 4L, 4:6)
    checkIdentical(target, pcompare(x2, y))
    checkIdentical(-target, pcompare(y, x2))

    x3 <- IRanges(8:16, width=2)
    target <- c(-6:-4, -1L, 2L, 3L, 4:6)
    checkIdentical(target, pcompare(x3, y))
    checkIdentical(-target, pcompare(y, x3))

    ## Moving a 0-width range over a non 0-width range.
    ## Note that when the end of the 0-width range is equal to the start of
    ## the non 0-width range minus 1, returning code -5 (which describes
    ## a situation of adjacent ranges) seems appropriate.
    ## However, one could argue that returning code -1 (which describes a
    ## situation where one range is inside the other) would also be
    ## appropriate, because, in that case, the two ranges have the same start.
    ## So the question really is whether the 0-width range should be considered
    ## *outside* or *inside* the non 0-width range.
    ## It's an arbitrary choice and we chose the former.
    x0 <- IRanges(10:16, width=0)
    target <- c(-6:-5, 2L, 2L, 2L, 5:6)
    checkIdentical(target, pcompare(x0, y))
    checkIdentical(-target, pcompare(y, x0))

    ## Moving a 0-width range over a 0-width range.
    y0 <- IRanges(13, 12)
    target <- c(-6L, -6L, -6L, 0L, 6L, 6L, 6L)
    checkIdentical(target, pcompare(x0, y0))
    checkIdentical(-target, pcompare(y0, x0))
}

test_order_IntegerRanges <- function()
{
    ir1 <- IRanges(c(2,5,1,5), c(3,7,3,6))
    ir1.sort <- IRanges(c(1,2,5,5), c(3,3,6,7))
    ir1.rev <- IRanges(c(5,5,2,1), c(7,6,3,3))
    checkIdentical(sort(ir1), ir1.sort)
    checkIdentical(sort(ir1, decreasing=TRUE), ir1.rev)
    checkException(sort(ir1, decreasing=NA), silent = TRUE)
}

