###

test_PartitioningByEnd <- function()
{
    ## on a numeric vector, NG not supplied
    current0 <- PartitioningByEnd()
    checkTrue(validObject(current0))

    target <- new("PartitioningByEnd")
    checkIdentical(target, current0)

    breakpoints <- c(0, 5, 5, 8)
    current1 <- PartitioningByEnd(breakpoints)
    checkTrue(validObject(current1))
    checkIdentical(4L, length(current1))
    checkIdentical(as.integer(breakpoints), end(current1))
    checkIdentical(end(current1), cumsum(width(current1)))
    checkIdentical(NULL, names(current1))

    checkException(PartitioningByEnd(breakpoints, names=letters), silent=TRUE)

    current2 <- PartitioningByEnd(breakpoints, names=letters[1:4])
    checkTrue(validObject(current2))
    checkIdentical(letters[1:4], names(current2))

    names(breakpoints) <- names(current2)
    current3 <- PartitioningByEnd(breakpoints)
    checkIdentical(current2, current3)

    current4 <- PartitioningByEnd(breakpoints, names=LETTERS[4:1])
    checkIdentical(LETTERS[4:1], names(current4))

    breakpoints <- rep.int(0, 1000)
    current5 <- PartitioningByEnd(breakpoints)
    checkTrue(validObject(current5))
    checkIdentical(as.integer(breakpoints), end(current5))
    checkIdentical(end(current5), cumsum(width(current5)))

    ## on a PartitioningByEnd object
    checkIdentical(current1, PartitioningByEnd(current1))  # no-op
    checkIdentical(current2, PartitioningByEnd(current2))  # no-op
    checkException(PartitioningByEnd(current2, names=LETTERS), silent=TRUE)
    current6 <- PartitioningByEnd(current2, names=names(current4))
    checkTrue(validObject(current6))
    checkIdentical(names(current4), names(current6))

    ## on CompressedList, SimpleList, IRanges, and list objects
    do_checks <- function(x) {
        checkIdentical(current1, PartitioningByEnd(x))
        checkException(PartitioningByEnd(x, names=letters), silent=TRUE)
        checkIdentical(current2, PartitioningByEnd(x, names=names(current2)))
        names(x) <- names(current2)
        checkIdentical(current2, PartitioningByEnd(x))
        checkIdentical(current4, PartitioningByEnd(x, names=names(current4)))
    }
    x <- RleList(Rle(), Rle(-3, 5), Rle(), Rle(1:0, c(2,1)), compress=TRUE)
    do_checks(x)
    do_checks(as(x, "SimpleList"))
    do_checks(as.list(x))
    x <- IRanges(seq(148, by=-50, length.out=4), width=width(current1))
    do_checks(x)
    ## TODO: Uncomment this when as.list() works again on IRanges objects
    #do_checks(as.list(x))
    do_checks(list(NULL, integer(5), complex(0), raw(3)))
}

test_PartitioningByWidth <- function()
{
    ## on a numeric vector, NG not supplied
    current0 <- PartitioningByWidth()
    checkTrue(validObject(current0))

    target <- new("PartitioningByWidth")
    checkIdentical(target, current0)

    widths <- c(0, 5, 0, 3)
    current1 <- PartitioningByWidth(widths)
    checkTrue(validObject(current1))
    checkIdentical(4L, length(current1))
    checkIdentical(as.integer(widths), width(current1))
    checkIdentical(end(current1), cumsum(width(current1)))
    checkIdentical(NULL, names(current1))

    checkException(PartitioningByWidth(widths, names=letters), silent=TRUE)

    current2 <- PartitioningByWidth(widths, names=letters[1:4])
    checkTrue(validObject(current2))
    checkIdentical(letters[1:4], names(current2))

    names(widths) <- names(current2)
    current3 <- PartitioningByWidth(widths)
    checkIdentical(current2, current3)

    current4 <- PartitioningByWidth(widths, names=LETTERS[4:1])
    checkIdentical(LETTERS[4:1], names(current4))

    widths <- rep.int(0, 1000)
    current5 <- PartitioningByWidth(widths)
    checkTrue(validObject(current5))
    checkIdentical(as.integer(widths), width(current5))
    checkIdentical(end(current5), cumsum(width(current5)))

    ## on a PartitioningByWidth object
    checkIdentical(current1, PartitioningByWidth(current1))  # no-op
    checkIdentical(current2, PartitioningByWidth(current2))  # no-op
    checkException(PartitioningByWidth(current2, names=LETTERS), silent=TRUE)
    current6 <- PartitioningByWidth(current2, names=names(current4))
    checkTrue(validObject(current6))
    checkIdentical(names(current4), names(current6))

    ## on CompressedList, SimpleList, IRanges, and list objects
    do_checks <- function(x) {
        checkIdentical(current1, PartitioningByWidth(x))
        checkException(PartitioningByWidth(x, names=letters), silent=TRUE)
        checkIdentical(current2, PartitioningByWidth(x, names=names(current2)))
        names(x) <- names(current2)
        checkIdentical(current2, PartitioningByWidth(x))
        checkIdentical(current4, PartitioningByWidth(x, names=names(current4)))
    }
    x <- RleList(Rle(), Rle(-3, 5), Rle(), Rle(1:0, c(2,1)), compress=TRUE)
    do_checks(x)
    do_checks(as(x, "SimpleList"))
    do_checks(as.list(x))
    x <- IRanges(seq(148, by=-50, length.out=4), width=width(current1))
    do_checks(x)
    ## TODO: Uncomment this when as.list() works again on IRanges objects
    #do_checks(as.list(x))
    do_checks(list(NULL, integer(5), complex(0), raw(3)))
}

test_PartitioningByEndOrWidth_NG_supplied <- function()
{
    for (class in c("PartitioningByEnd", "PartitioningByWidth")) {
        CONSTRUCTOR <- get(class)
        x <- c(3, 3, 4, 6)
        NG <- 8
        current1 <- CONSTRUCTOR(x, NG)
        checkTrue(is(current1, class))
        checkTrue(validObject(current1))
        checkIdentical(8L, length(current1))
        checkIdentical(tabulate(x, nbins=NG), width(current1))

        checkException(CONSTRUCTOR(x, NG, names=letters[1:4]), silent=TRUE)

        current2 <- CONSTRUCTOR(x, NG, names=letters[1:8])
        checkTrue(validObject(current2))
        checkIdentical(letters[1:8], names(current2))

        names(x) <- letters[1:4]
        current3 <- CONSTRUCTOR(x, NG)
        checkIdentical(current1, current3)
    }
}

