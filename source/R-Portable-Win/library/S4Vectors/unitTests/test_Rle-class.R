test_Rle_construction <- function() {
    empty <- Rle()
    checkTrue(validObject(empty))
    checkIdentical(Rle(), new("Rle"))
    checkIdentical(length(empty), 0L)
    x <- Rle(rep(6:10, 1:5))
    checkTrue(validObject(x))
    checkIdentical(x, Rle(6:10, 1:5))
    y <- Rle(factor(rep(letters, 1:26)))
    checkTrue(validObject(y))
    checkIdentical(y, Rle(factor(letters), 1:26))

    checkIdentical(Rle(c(TRUE, TRUE, FALSE, FALSE, FALSE, NA, NA, NA)),
                   Rle(c(TRUE, FALSE, NA), c(2, 3, 3)))
    checkIdentical(Rle(c(1L, 1L, 1L, 2L, 2L, NA, NA, NA)),
                   Rle(c(1L, 2L, NA), c(3, 2, 3)))
    checkIdentical(Rle(c(1, 1, 1, 2, 2, NA, NA, NA)),
                   Rle(c(1, 2, NA), c(3, 2, 3)))
    checkIdentical(Rle(c("a", "a", "b", "b", "b", NA, NA, NA)),
                   Rle(c("a", "b", NA), c(2, 3, 3)))
}

test_Rle_replace <- function() {
    x <- Rle(1:26, 1:26)
    runValue(x) <- letters
    checkTrue(validObject(x))
    checkIdentical(x, Rle(letters, 1:26))
    runLength(x) <- 26:1
    checkTrue(validObject(x))
    checkIdentical(x, Rle(letters, 26:1))
}

test_Rle_coercion <- function() {
    x <- rep(6:10, 1:5)
    xRle <- Rle(x)
    y <- c(TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE)
    yRle <- Rle(y)
    checkIdentical(x, as.vector(xRle))
    checkIdentical(as.integer(x), as.integer(xRle))
    checkIdentical(as.numeric(x), as.numeric(xRle))
    checkIdentical(as.complex(x), as.complex(xRle))
    checkIdentical(as.factor(x), as.factor(xRle))
    checkIdentical(y, as.vector(yRle))
    checkIdentical(as.logical(y), as.logical(yRle))
    checkIdentical(as.character(y), as.character(yRle))
    checkIdentical(as.raw(y), as.raw(yRle))
    checkIdentical(as.factor(y), as.factor(yRle))
}

test_extract_ranges_from_Rle <- function() {
    extract_ranges_from_Rle <- S4Vectors:::extract_ranges_from_Rle

    # Extract single range.
    x <- Rle()
    for (method in 0:3) {
        current <- extract_ranges_from_Rle(x, 1L, 0L, method)
        checkIdentical(x, current)
        checkException(extract_ranges_from_Rle(x, 1L, 1L, method), silent=TRUE)
        checkException(extract_ranges_from_Rle(x, 0L, 0L, method), silent=TRUE)
        checkException(extract_ranges_from_Rle(x, 0L, 1L, method), silent=TRUE)
    }

    x <- Rle(0.8, 10L)
    for (method in 0:3) {
        target <- Rle(numeric(0))
        for (start in 1:11) {
            current <- extract_ranges_from_Rle(x, start, 0L, method)
            checkIdentical(target, current)
        }
        checkException(extract_ranges_from_Rle(x, 0L, 0L, method), silent=TRUE)
        checkException(extract_ranges_from_Rle(x, 12L, 1L, method), silent=TRUE)

        target <- Rle(0.8)
        for (start in 1:10) {
            current <- extract_ranges_from_Rle(x, start, 1L, method)
            checkIdentical(target, current)
        }
        checkException(extract_ranges_from_Rle(x, 0L, 1L, method), silent=TRUE)
        checkException(extract_ranges_from_Rle(x, 11L, 1L, method), silent=TRUE)
    }

    # Extract multiple ranges.
    x <- Rle(factor(letters[1:3], levels=rev(letters)), 7:5)

    start <- 1L
    width <- length(x)
    for (method in 0:3) {
        current <- extract_ranges_from_Rle(x, start, width, method)
        checkIdentical(x, current)
    }

    start <- seq_along(x)
    width <- rep(1L, length(start))
    for (method in 0:3) {
        current <- extract_ranges_from_Rle(x, start, width, method)
        checkIdentical(x, current)
    }

    start <- seq_len(length(x) + 1L)
    width <- rep(0L, length(start))
    target <- Rle(factor(levels=rev(letters)))
    for (method in 0:3) {
        current <- extract_ranges_from_Rle(x, start, width, method)
        checkIdentical(target, current)
    }

    start <- seq_len(length(x) - 5L)
    width <- rep(c(6L, 2L, 7L), length.out=length(start))
    target <- S4Vectors:::extract_ranges_from_vector_OR_factor(
                                  S4Vectors:::decodeRle(x), start, width)
    for (method in 0:3) {
        current <- extract_ranges_from_Rle(x, start, width, method)
        checkIdentical(target, S4Vectors:::decodeRle(current))
    }

    start <- rev(start)
    width <- rev(width)
    target <- S4Vectors:::extract_ranges_from_vector_OR_factor(
                                  S4Vectors:::decodeRle(x), start, width)
    for (method in 0:3) {
        current <- extract_ranges_from_Rle(x, start, width, method)
        checkIdentical(target, S4Vectors:::decodeRle(current))
    }
}

test_Rle_general <- function() {
    x <- rep(6:10, 1:5)
    xRle <- Rle(x)
    checkIdentical(unique(x), unique(xRle))
    checkIdentical(x[c(3,2,4,6)], as.vector(xRle[c(3,2,4,6)]))
    checkIdentical(append(x,x), as.vector(append(xRle,xRle)))
    checkIdentical(append(x,x,3), as.vector(append(xRle,xRle,3)))
    checkIdentical(c(x,x) %in% c(7:9), as.vector(c(xRle,xRle)) %in% c(7:9))
    checkIdentical(c(x, x), as.vector(c(xRle, xRle)))
    checkIdentical(is.na(c(NA, x, NA, NA, NA, x, NA)),
                   as.vector(is.na(c(Rle(NA), xRle, Rle(NA, 3), xRle, Rle(NA)))))
    checkIdentical(is.unsorted(c(1,2,2,3)), is.unsorted(Rle(c(1,2,2,3))))
    checkIdentical(is.unsorted(c(1,2,2,3), strictly = TRUE),
                   is.unsorted(Rle(c(1,2,2,3)), strictly = TRUE))
    checkIdentical(length(x), length(xRle))
    checkIdentical(match(c(x,x), c(7:9)), as.vector(match(c(xRle,xRle), c(7:9))))
    checkIdentical(rep(x, times = 2), as.vector(rep(xRle, times = 2)))
    checkIdentical(rep(x, times = x), as.vector(rep(xRle, times = x)))
    checkIdentical(rep(x, length.out = 20), as.vector(rep(xRle, length.out = 20)))
    checkIdentical(rep(x, each = 2), as.vector(rep(xRle, each = 2)))
    checkIdentical(rep(x, x, 20), as.vector(rep(xRle, x, 20)))
    checkException(rep(xRle, x, each = 2), silent = TRUE)
    checkIdentical(rep(x, 2, each = 2), as.vector(rep(xRle, 2, each = 2)))
    checkIdentical(rep(x, length.out = 20, each = 2),
                   as.vector(rep(xRle, length.out = 20, each = 2)))
    checkIdentical(rep(x, x, 20, 2), as.vector(rep(xRle, x, 20, 2)))
    checkIdentical(rep.int(x, times = 2), as.vector(rep.int(xRle, times = 2)))
    checkIdentical(rev(x), as.vector(rev(xRle)))

    library(IRanges)
    checkIdentical(as.vector(xRle[IRanges(start=1:3, width=1:3)]),
                   x[c(1,2,3,3,4,5)])
    z <- x
    z[] <- rev(z)
    zRle <- xRle
    zRle[] <- rev(zRle)
    checkIdentical(z, as.vector(zRle))
    z <- x
    z[c(1,5,3)] <- 3:1
    zRle <- xRle
    zRle[c(1,5,3)] <- 3:1
    checkIdentical(z, as.vector(zRle))
    z <- x
    z[1:5] <- 0L
    zRle <- xRle
    zRle[IRanges(start=1:3, width=1:3)] <- 0L
    checkIdentical(z, as.vector(zRle))
    checkIdentical(sort(c(x,x)), as.vector(sort(c(xRle,xRle))))

    checkIdentical(as.vector(subset(xRle, rep(c(TRUE, FALSE), length.out = length(.(x))))),
                   subset(x, rep(c(TRUE, FALSE), length.out = length(x))))
    checkIdentical(as.vector(window(x, start = 3, end = 13)),
                   as.vector(window(xRle, start = 3, end = 13)))
    z <- x
    z[3:13] <- 0L
    zRle <- xRle
    window(zRle, start = 3, end = 13) <- 0L
    checkIdentical(z, as.vector(zRle))
}

## ---------------------------------------------
## table() and sort()
## ---------------------------------------------

test_Rle_sort <- function()
{
    ## atomic
    ix <- c(NA, 3L, NA)
    nx <- c(2, 5, 1, 2, NA, 5, NA)
    cx <- c("c", "B", NA, "a")
    lx <- c(FALSE, FALSE, NA, TRUE, NA)
    checkIdentical(sort(nx), as.numeric(sort(Rle(nx))))
    checkIdentical(sort(nx, na.last=TRUE), 
        as.numeric(sort(Rle(nx), na.last=TRUE)))
    checkIdentical(sort(nx, na.last=FALSE), 
        as.numeric(sort(Rle(nx), na.last=FALSE)))
    checkIdentical(sort(ix), as.integer(sort(Rle(ix))))
    checkIdentical(sort(cx), as.character(sort(Rle(cx))))
    checkIdentical(sort(lx), as.logical(sort(Rle(lx))))
    checkIdentical(sort(numeric()), as.numeric(sort(Rle(numeric()))))
    checkIdentical(sort(character()), as.character(sort(Rle(character()))))
    
    ## factor 
    nf <- factor(nx)
    checkIdentical(sort(nf), as.factor(sort(Rle(nf))))
    checkIdentical(sort(nf, decreasing=TRUE, na.last=TRUE), 
        as.factor(sort(Rle(nf), decreasing=TRUE, na.last=TRUE)))
    checkIdentical(sort(nf, na.last=FALSE), 
        as.factor(sort(Rle(nf), na.last=FALSE)))
    checkIdentical(sort(factor()), as.factor(sort(Rle(factor()))))

    ## factor, unused levels
    nf <- factor(nx, levels=1:6)
    checkIdentical(levels(sort(nf)), levels(sort(Rle(nf))))
}

test_Rle_table <- function()
{
    ## atomic
    ix <- c(NA, 3L, NA)
    nx <- c(2, 5, 1, 2, NA, 5, NA)
    cx <- c("c", "B", NA, "a")
    lx <- c(FALSE, FALSE, NA, TRUE, NA)
    checkIdentical(table(ix), table("ix"=Rle(ix)))
    checkIdentical(table(nx), table("nx"=Rle(nx)))
    checkIdentical(table(cx), table("cx"=Rle(cx)))
    checkIdentical(table(lx), table("lx"=Rle(lx)))
    checkIdentical(table(numeric()), table(Rle(numeric())))
    checkIdentical(table(character()), table(Rle(character())))
    
    ## factor
    nf <- factor(nx)
    checkIdentical(table("nx"=nx), table("nx"=Rle(nx)))
    checkIdentical(table(factor()), table(Rle(factor())))
    
    ## factor, unused levels
    nf <- factor(nx, levels=1:6)
    cf <- factor(cx, levels=c("a", "c", "B", "b"))
    checkIdentical(as.factor(table(nf)), as.factor(table(Rle(nf))))
    checkIdentical(as.factor(table(cf)), as.factor(table(Rle(cf))))
}

test_Rle_Integer_overflow <- function() {
    v <- as.integer(c(1,(2^31)-1,1))
    x0 <- Rle(v)
    checkIdentical(sum(v), sum(x0))

    x <- Rle(c(1,(2^31)-1,1))
    checkIdentical(mean(x0), mean(x))
}

