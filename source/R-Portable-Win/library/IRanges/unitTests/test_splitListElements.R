###

test_regroupBySupergroup <- function()
{
    regroupBySupergroup <- IRanges:::regroupBySupergroup

    .do_checks <- function(x, breakpoints, target)
    {
        supergroups <- PartitioningByEnd(breakpoints)

        current <- regroupBySupergroup(x, supergroups)
        checkIdentical(target, current)
        checkIdentical(target, regroupBySupergroup(x, breakpoints))

        x_partitioning <- PartitioningByEnd(x)
        current2 <- regroupBySupergroup(x_partitioning, supergroups)
        checkIdentical(PartitioningByEnd(target), current2)
    }

    x <- CharacterList(
           x1=NULL,
           x2=LETTERS[1:3],
           x3=LETTERS[4:5],
           x4=letters[1:5],
           x5=NULL,
           x6=letters[6:7]
    )

    breakpoints <- c(SG1=3, SG2=6)
    target <- CharacterList(SG1=LETTERS[1:5], SG2=letters[1:7],
                            compress=TRUE)
    .do_checks(as(x, "CompressedList"), breakpoints, target)
    .do_checks(as(x, "SimpleList"), breakpoints, target)

    breakpoints <- c(SG1=2, SG2=5, SG3=6)
    target <- CharacterList(SG1=LETTERS[1:3],
                            SG2=c(LETTERS[4:5], letters[1:5]),
                            SG3=letters[6:7],
                            compress=TRUE)
    .do_checks(as(x, "CompressedList"), breakpoints, target)
    .do_checks(as(x, "SimpleList"), breakpoints, target)

    breakpoints <- c(SG1=2, 2, SG2=5, SG3=6)
    target <- CharacterList(SG1=LETTERS[1:3],
                            NULL,
                            SG2=c(LETTERS[4:5], letters[1:5]),
                            SG3=letters[6:7],
                            compress=TRUE)
    .do_checks(as(x, "CompressedList"), breakpoints, target)
    .do_checks(as(x, "SimpleList"), breakpoints, target)

    breakpoints <- 6
    target <- CharacterList(unlist(x, use.names=FALSE), compress=TRUE)
    .do_checks(as(x, "CompressedList"), breakpoints, target)
    .do_checks(as(x, "SimpleList"), breakpoints, target)

    breakpoints <- c(SG1=6, SG2=6, SG3=6)
    target <- CharacterList(SG1=unlist(x, use.names=FALSE),
                            SG2=NULL,
                            SG3=NULL,
                            compress=TRUE)
    .do_checks(as(x, "CompressedList"), breakpoints, target)
    .do_checks(as(x, "SimpleList"), breakpoints, target)

    breakpoints <- c(0, 0, 0, 6, 6)
    target <- CharacterList(NULL, NULL, NULL,
                            unlist(x, use.names=FALSE),
                            NULL,
                            compress=TRUE)
    .do_checks(as(x, "CompressedList"), breakpoints, target)
    .do_checks(as(x, "SimpleList"), breakpoints, target)

    breakpoints <- seq_along(x)
    target <- unname(as(x, "CompressedList"))
    .do_checks(as(x, "CompressedList"), breakpoints, target)
    .do_checks(as(x, "SimpleList"), breakpoints, target)

    names(breakpoints) <- names(x)
    target <- as(x, "CompressedList")
    .do_checks(as(x, "CompressedList"), breakpoints, target)  # no-op
    .do_checks(as(x, "SimpleList"), breakpoints, target)

    x0 <- CharacterList()
    breakpoints <- setNames(integer(0), character(0))
    target <- setNames(CharacterList(compress=TRUE), character(0))
    .do_checks(as(x0, "CompressedList"), breakpoints, target)
    # Fails at the moment because unlist() is doesn't work properly on
    # SimpleCharacterList
    #.do_checks(as(x0, "SimpleList"), breakpoints, target)

    x2 <- RleList(Rle(44:45, 2:1), Rle(45), Rle(-2, 3))
    breakpoints <- c(SG1=2, SG2=3)
    target <- RleList(SG1=Rle(44:45, c(2,2)), SG2=Rle(-2, 3),
                      compress=TRUE)
    .do_checks(as(x2, "CompressedList"), breakpoints, target)
    .do_checks(as(x2, "SimpleList"), breakpoints, target)

    x3 <- Views(unlist(x2, use.names=FALSE),
                start=c(3, 1, 1), end=c(6, 1, 3))
    breakpoints <- c(SG1=2, SG2=3)
    target <- RleList(SG1=Rle(c(45,-2,44), c(2,2,1)), SG2=Rle(44:45, 2:1),
                      compress=TRUE)
    .do_checks(x3, breakpoints, target)
}

