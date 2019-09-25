library(IRanges)  # many tests in this file use functionalities defined
                  # in IRanges

test_Rle_groupGeneric <- function() {
    set.seed(0)
    x <- sample(0:3, 50, replace = TRUE)
    xRle <- Rle(x)
    checkIdentical(numeric(0) + 1, as.vector(Rle(numeric(0)) + 1))
    checkIdentical(x + 1, as.vector(xRle + 1))
    checkIdentical(2 * x + 3, as.vector(2 * xRle + 3))    
    checkIdentical(x[(x > 0) & (x < 3)], as.vector(xRle[(xRle > 0) & (xRle < 3)]))
    checkIdentical(log(x), as.vector(log(xRle)))
    checkIdentical(range(x), range(xRle))
    checkIdentical(sum(x), sum(xRle))
    checkIdentical(prod(x), prod(xRle))
    checkIdentical(cumsum(x), as.vector(cumsum(xRle)))
    checkIdentical(cumprod(x), as.vector(cumprod(xRle)))
    checkIdentical(round(x + .25), as.vector(round(xRle + .25)))
    checkIdentical(signif(x + .25), as.vector(signif(xRle + .25)))
    checkIdentical(Im(x + 5i), as.vector(Im(xRle + 5i)))
}

test_Rle_general <- function() {
    x <- rep(6:10, 1:5)
    xRle <- Rle(x)
    checkIdentical(aggregate(xRle, IRanges(start = 3:6, end = 13:10), FUN = mean),
                   aggregate(xRle, FUN = mean, start = 3:6, width = seq(11, 5, by = -2)))
    exp <- c(mean(x[3:13]), mean(x[4:12]), mean(x[5:11]), mean(x[6:10]))
    agg <- aggregate(xRle, FUN = function(x) x, start = 3:6, end = 13:10)
    checkEquals(exp, aggregate(xRle, FUN = mean, start = 3:6, end = 13:10))
    checkEquals(as.vector(aggregate.ts(ts(x, frequency = 5), FUN = mean)),
                aggregate(xRle, FUN = mean, start = c(1, 6, 11), end = c(5, 10, 15)))
    
    #checkIdentical(findRange(c(1, 3, 5), xRle), IRanges(start = c(1,2,4), width = 1:3))
    #checkIdentical(head(x, 8), as.vector(head(xRle, 8)))
    #checkIdentical(head(x, -3), as.vector(head(xRle, -3)))

    #checkException(split(Rle(1:26), integer()), silent = TRUE)
    #checkException(split(Rle(1:26), Rle()), silent = TRUE)
    #checkIdentical(lapply(as.list(split(Rle(1:26), letters)), as.vector),
    #               split(1:26, letters))
    #checkIdentical(lapply(as.list(split(Rle(1:26), Rle(letters))), as.vector),
    #               split(1:26, letters))
    #checkIdentical(lapply(as.list(split(Rle(1:26), letters[1:2])), as.vector),
    #               split(1:26, letters[1:2]))
    #checkIdentical(lapply(as.list(split(Rle(1:26), Rle(letters[1:2]))), as.vector),
    #               split(1:26, letters[1:2]))
    #checkIdentical(lapply(as.list(split(Rle(integer()), letters)), as.vector),
    #               split(integer(), letters))
    #checkIdentical(lapply(as.list(split(Rle(integer()), Rle(letters))), as.vector),
    #               split(integer(), letters))

    #checkIdentical(splitRanges(Rle(letters, 1:26)),
    #               split(IRanges(end = cumsum(1:26), width = 1:26), letters))

    checkIdentical(summary(x), summary(xRle))
    #checkIdentical(tail(x, 8), as.vector(tail(xRle, 8)))
    #checkIdentical(tail(x, -3), as.vector(tail(xRle, -3)))
    #checkException(tapply(xRle), silent = TRUE)
    #checkIdentical(tapply(x, x), tapply(xRle, xRle))
    #checkIdentical(tapply(x, x, mean), tapply(xRle, xRle, mean))
    #checkIdentical(tapply(xRle, x, mean), tapply(xRle, xRle, mean))
    #checkIdentical(tapply(x, x, mean, simplify = FALSE),
    #               tapply(xRle, xRle, mean, simplify = FALSE))
    #checkIdentical(tapply(xRle, x, mean, simplify = FALSE),
    #               tapply(xRle, xRle, mean, simplify = FALSE))
}

test_Rle_logical <- function() {
    checkIdentical(logical(), as.vector(Rle(logical())))

    x <- c(TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE)
    xRle <- Rle(x)
    checkIdentical(!x, as.vector(!x))
    checkIdentical(which(x), as.vector(which(x)))
    checkIdentical(as(xRle, "IRanges"),
                   IRanges(start = c(1,5,7), width = c(2, 1, 3)))
}

test_Rle_numerical <- function() {
    checkIdentical(numeric(), as.vector(Rle(numeric())))

    x <- cumsum(cumsum(1:10))
    xRle <- Rle(x)
    checkIdentical(pmax(x, rev(x)), as.vector(pmax(xRle, rev(xRle))))
    checkIdentical(pmin(x, rev(x)), as.vector(pmin(xRle, rev(xRle))))
    checkIdentical(pmax.int(x, rev(x)), as.vector(pmax.int(xRle, rev(xRle))))
    checkIdentical(pmin.int(x, rev(x)), as.vector(pmin.int(xRle, rev(xRle))))
    checkIdentical(diff(x), as.vector(diff(xRle)))
    checkIdentical(diff(x, lag = 2), as.vector(diff(xRle, lag = 2)))
    checkIdentical(diff(x, differences = 2), as.vector(diff(xRle, differences = 2)))
    checkIdentical(diff(x, lag = 2, differences = 2), 
                   as.vector(diff(xRle, lag = 2, differences = 2)))

    x <- rep(c(1.2, 3.4, NA, 7.8, 9.0), 1:5)
    y <- x - rev(x)
    xRle <- Rle(x)
    yRle <- Rle(y)
    checkIdentical(mean(x), mean(xRle))
    checkIdentical(mean(x, na.rm = TRUE), mean(xRle, na.rm = TRUE))
    checkIdentical(var(x), var(xRle))
    checkEqualsNumeric(var(x, na.rm = TRUE), var(xRle, na.rm = TRUE))
    checkIdentical(var(x, y), var(xRle, yRle))
    checkEqualsNumeric(var(x, y, na.rm = TRUE), var(xRle, yRle, na.rm = TRUE))
    checkIdentical(cov(x, y), cov(xRle, yRle))
    checkEqualsNumeric(cov(x, y, use = "complete"), cov(xRle, yRle, use = "complete"))
    checkIdentical(cor(x, y), cor(xRle, yRle))
    checkEqualsNumeric(cor(x, y, use = "complete"), cor(xRle, yRle, use = "complete"))
    checkIdentical(sd(x), sd(xRle))
    checkEqualsNumeric(sd(x, na.rm = TRUE), sd(xRle, na.rm = TRUE))
    checkIdentical(median(x), median(xRle))
    checkIdentical(median(x, na.rm = TRUE), median(xRle, na.rm = TRUE))
    checkIdentical(quantile(x, na.rm = TRUE), quantile(xRle, na.rm = TRUE))
    checkIdentical(mad(x), mad(xRle))
    checkIdentical(mad(x, na.rm = TRUE), mad(xRle, na.rm = TRUE))
    checkIdentical(IQR(x, na.rm = TRUE), IQR(xRle, na.rm = TRUE))

    y <- (-20:20)^2
    y[c(1,10,21,41)] <- c(100L, 30L, 400L, 470L)
    checkEqualsNumeric(smoothEnds(y), as.vector(smoothEnds(Rle(y))))
    checkEqualsNumeric(runmed(y, 7), as.vector(runmed(Rle(y), 7)))
    checkEqualsNumeric(runmed(y, 11), as.vector(runmed(Rle(y), 11)))
    checkEqualsNumeric(runmed(y, 7, "keep"),
                       as.vector(runmed(Rle(y), 7, "keep")))
    checkEqualsNumeric(runmed(y, 11, "keep"),
                       as.vector(runmed(Rle(y), 11, "keep")))
    checkEqualsNumeric(runmed(y, 7, "constant"),
                       as.vector(runmed(Rle(y), 7, "constant")))
    checkEqualsNumeric(runmed(y, 11, "constant"),
                       as.vector(runmed(Rle(y), 11, "constant")))

    x <- rep(c(1.2, 3.4, 5.6, 7.8, 9.0), 1:5)
    y <- rep(1:5, c(4, 2, 5, 1, 3))
    xRle <- Rle(x)
    yRle <- Rle(y)
    checkEqualsNumeric(sapply(1:13, function(i) sum(window(x, i, i + 2))),
                       as.numeric(runsum(xRle, k = 3)))
#    checkEqualsNumeric(sapply(1:13, function(i) sum(window(rev(x), i, i + 2))),
#                       as.numeric(runsum(rev(xRle), k = 3)))
    checkEqualsNumeric(sapply(1:13, function(i) sum(window(y, i, i + 2))),
                       as.integer(runsum(yRle, k = 3)))
    checkEqualsNumeric(sapply(1:13, function(i) sum(window(rev(y), i, i + 2))),
                       as.integer(runsum(rev(yRle), k = 3)))
    checkEqualsNumeric(sapply(1:13, function(i) mean(window(x, i, i + 2))),
                       as.numeric(runmean(xRle, k = 3)))
    checkEqualsNumeric(sapply(1:13, function(i) mean(window(rev(x), i, i + 2))),
                       as.numeric(runmean(rev(xRle), k = 3)))
    checkEqualsNumeric(sapply(1:13, function(i) mean(window(y, i, i + 2))),
                       as.numeric(runmean(yRle, k = 3)))
    checkEqualsNumeric(sapply(1:13, function(i) mean(window(rev(y), i, i + 2))),
                       as.numeric(runmean(rev(yRle), k = 3)))
    checkEqualsNumeric(sapply(1:13, function(i) sum(window(x, i, i + 2))),
                       as.numeric(runwtsum(xRle, k = 3, wt = rep(1,3))))
    checkEqualsNumeric(sapply(1:13, function(i) sum(window(x, i, i + 2)/3)),
                       as.numeric(runwtsum(xRle, k = 3, wt = rep(1/3,3))))
    checkEqualsNumeric(sapply(1:13, function(i) sum(window(y, i, i + 2))),
                       as.numeric(runwtsum(yRle, k = 3, wt = rep(1,3))))
    checkEqualsNumeric(sapply(1:13, function(i) sum(window(y, i, i + 2)/3)),
                       as.numeric(runwtsum(yRle, k = 3, wt = rep(1/3,3))))
    checkEqualsNumeric(sapply(1:13, function(i) min(window(x, i, i + 2))),
                       as.numeric(runq(xRle, k = 3, i = 1)))
    checkEqualsNumeric(sapply(1:13, function(i) median(window(x, i, i + 2))),
                       as.numeric(runq(xRle, k = 3, i = 2)))
    checkEqualsNumeric(sapply(1:13, function(i) max(window(x, i, i + 2))),
                       as.numeric(runq(xRle, k = 3, i = 3)))
    checkIdentical(runq(xRle, k = 3, i = 2),
                   rev(runq(rev(xRle), k = 3, i = 2)))
    checkEqualsNumeric(sapply(1:13, function(i) min(window(y, i, i + 2))),
                       as.numeric(runq(yRle, k = 3, i = 1)))
    checkEqualsNumeric(sapply(1:13, function(i) median(window(y, i, i + 2))),
                       as.numeric(runq(yRle, k = 3, i = 2)))
    checkEqualsNumeric(sapply(1:13, function(i) max(window(y, i, i + 2))),
                       as.numeric(runq(yRle, k = 3, i = 3)))
    checkIdentical(runq(yRle, k = 3, i = 2),
                   rev(runq(rev(yRle), k = 3, i = 2)))
}

test_Rle_character <- function() {
    checkIdentical(character(), as.vector(Rle(character())))

    txt <-
      c("The", "licenses", "for", "most", "software", "are", "designed",
        "to", "take", "away", "your", "freedom", "to", "share", "and",
        "change", "it.", "", "By", "contrast,", "the", "GNU", "General",
        "Public", "License", "is", "intended", "to", "guarantee", "your",
        "freedom", "to", "share", "and", "change", "free", "software",
        "--", "to", "make", "sure", "the", "software", "is", "free", "for",
        "all", "its", "users")
     txt <- rep(txt, seq_len(length(txt)))
     txtRle <- Rle(txt)
     checkIdentical(nchar(txt), as.vector(nchar(txtRle)))
     checkIdentical(substr(txt, 3, 7), as.vector(substr(txtRle, 3, 7)))
     checkIdentical(substring(txt, 4, 9), as.vector(substring(txtRle, 4, 9)))
     checkIdentical(chartr("@!*", "alo", txt),
                    as.vector(chartr("@!*", "alo", txtRle)))
     checkIdentical(tolower(txt), as.vector(tolower(txtRle)))
     checkIdentical(toupper(txt), as.vector(toupper(txtRle)))
     checkIdentical(sub("[b-e]",".", txt), as.vector(sub("[b-e]",".", txtRle)))
     checkIdentical(gsub("[b-e]",".", txt), as.vector(gsub("[b-e]",".", txtRle)))
     checkIdentical(paste(txt, rev(txt), sep = "|"),
                    as.vector(paste(txtRle, rev(txtRle), sep = "|")))

     modifyFactor <- function(x, FUN, ...) {
         levels(x) <- FUN(levels(x), ...)
         x
     }
     fac <- factor(txt)
     facRle <- Rle(fac)
     checkIdentical(modifyFactor(fac, substr, 3, 7),
                    as.factor(substr(facRle, 3, 7)))
     checkIdentical(modifyFactor(fac, substring, 4, 9),
                    as.factor(substring(facRle, 4, 9)))
     checkIdentical(modifyFactor(fac, chartr, old = "@!*", new = "alo"),
                    as.factor(chartr("@!*", "alo", facRle)))
     checkIdentical(modifyFactor(fac, tolower), as.factor(tolower(facRle)))
     checkIdentical(modifyFactor(fac, toupper), as.factor(toupper(facRle)))
     checkIdentical(modifyFactor(fac, sub, pattern = "[b-e]",
                                 replacement = "."),
                    as.factor(sub("[b-e]",".", facRle)))
     checkIdentical(modifyFactor(fac, gsub, pattern = "[b-e]",
                                 replacement = "."),
                    as.factor(gsub("[b-e]",".", facRle)))
     checkTrue(is.factor(runValue(paste(facRle, rev(facRle), sep = "|"))))
}

test_Rle_factor <- function() {
    checkIdentical(factor(character()),
                   as.factor(Rle(factor(character()))))

    x <- factor(rep(letters, 1:26))
    xRle <- Rle(x)
    checkIdentical(levels(x), levels(xRle))
    levels(x) <- LETTERS
    levels(xRle) <- LETTERS
    checkIdentical(levels(x), levels(xRle))
    checkIdentical(nlevels(x), 26L)
    xRle[] <- xRle
    checkIdentical(Rle(x), xRle)
    checkIdentical(x, xRle[TRUE,drop=TRUE])
}

## ---------------------------------------------
## runsum(), runmean(), runwtsum()
## ---------------------------------------------

.naive_runsum <- function(x, k, na.rm=FALSE)
    sapply(0:(length(x)-k),
        function(offset) sum(x[1:k + offset], na.rm=na.rm)) 

checkIdenticalIfNaNsWereNAs <- function(x, y)
{
    x[is.nan(x)] <- NA_real_
    y[is.nan(y)] <- NA_real_
    checkIdentical(x, y)
}

test_Rle_runsum_real <- function() {

    x0 <- c(NA, NaN, Inf, -Inf) 
    x <- Rle(x0)
    ## na.rm = TRUE 
    target1 <- .naive_runsum(x0, 4, na.rm=TRUE)
    target2 <- .naive_runsum(x, 4, na.rm=TRUE)
    checkIdenticalIfNaNsWereNAs(target1, target2) 
    current <- as.vector(runsum(x, 4, na.rm=TRUE))
    checkIdenticalIfNaNsWereNAs(target1, current)
    ## na.rm = FALSE 
    target1 <- .naive_runsum(x0, 4, na.rm=FALSE)
    target2 <- .naive_runsum(x, 4, na.rm=FALSE)
    checkIdenticalIfNaNsWereNAs(target1, target2) 
    current <- as.vector(runsum(x, 4, na.rm=FALSE))
    checkIdenticalIfNaNsWereNAs(target1, current) 

    x0 <- c(NA, Inf, NA, -Inf, Inf, -Inf, NaN, Inf, NaN, -Inf)
    x <- Rle(x0)
    for (k in 1:2) {
        target1 <- .naive_runsum(x0, k, na.rm=TRUE)
        target2 <- .naive_runsum(x, k, na.rm=TRUE)
        checkIdenticalIfNaNsWereNAs(target1, target2)
        current <- as.vector(runsum(x, k, na.rm=TRUE))
        checkIdenticalIfNaNsWereNAs(target1, current) 

        target1 <- .naive_runsum(x0, k, na.rm=FALSE)
        target2 <- .naive_runsum(x, k, na.rm=FALSE)
        checkIdenticalIfNaNsWereNAs(target1, target2)
        current <- as.vector(runsum(x, k, na.rm=FALSE))
        checkIdenticalIfNaNsWereNAs(target1, current)
    }
 
    ## NOTE : Inconsistent behavior in base::sum()
    ## sum(x, y) and x + y:
    ## > sum(NaN, NA)
    ##   [1] NA
    ## > NaN + NA
    ##   [1] NaN
    ## also between sum(c(x, y)) and sum(x, y):
    ## This inconsistency only exists on linux, not Mac or Windows
    ##  > sum(c(NaN, NA))
    ##  [1] NaN
    ##  > sum(NaN, NA)
    ##  [1] NA 
    ## x0 <- c(NA, NaN, NA)
    ## x <- Rle(x0)
    ## target1 <- c(x0[1] + x0[2], x0[2] + x0[3]) 
    ## target2 <- as.vector(c(x[1] + x[2], x[2] + x[3]))
    ## checkIdentical(target1, target2)
    ## current <- as.vector(runsum(x, k=2, na.rm=FALSE))
    ## checkIdentical(target1, current)
}

test_Rle_runsum_integer <- function() {

    x0 <- c(NA_integer_, 1L, 1L)
    x <- Rle(x0)
    for (k in 1:3) {
        target1 <- .naive_runsum(x0, k, na.rm=TRUE)
        target2 <- .naive_runsum(x, k, na.rm=TRUE)
        checkIdentical(target1, target2) 
        current <- as.vector(runsum(x, k, na.rm=TRUE))
        checkIdentical(target1, current)

        target1 <- .naive_runsum(x0, k, na.rm=FALSE)
        target2 <- .naive_runsum(x, k, na.rm=FALSE)
        checkIdentical(target1, target2) 
        current <- as.vector(runsum(x, k, na.rm=FALSE))
        checkIdentical(target1, current)
    }

    x0 <- c(1L, NA_integer_, 1L)
    x <- Rle(x0)
    for (k in 1:3) {
        target1 <- .naive_runsum(x0, k, na.rm=TRUE)
        target2 <- .naive_runsum(x, k, na.rm=TRUE)
        checkIdentical(target1, target2) 
        current <- as.vector(runsum(x, k, na.rm=TRUE))
        checkIdentical(target1, current)

        target1 <- .naive_runsum(x0, k, na.rm=FALSE)
        target2 <- .naive_runsum(x, k, na.rm=FALSE)
        checkIdentical(target1, target2) 
        current <- as.vector(runsum(x, k, na.rm=FALSE))
        checkIdentical(target1, current)
    }
}

.naive_runmean <- function(x, k, na.rm=FALSE)
    sapply(0:(length(x)-k),
        function(offset) mean(x[1:k + offset], na.rm=na.rm)) 

test_Rle_runmean <- function() {

    x0 <- c(NA, 1, 1)
    x <- Rle(x0)
    for (k in 1:3) {
        target1 <- .naive_runmean(x0, k, na.rm=TRUE)
        target2 <- .naive_runmean(x, k, na.rm=TRUE)
        checkIdentical(target1, target2) 
        current <- as.vector(runmean(x, k, na.rm=TRUE))
        checkIdentical(target1, current)

        target1 <- .naive_runmean(x0, k, na.rm=FALSE)
        target2 <- .naive_runmean(x, k, na.rm=FALSE)
        checkIdentical(target1, target2) 
        current <- as.vector(runmean(x, k, na.rm=FALSE))
        checkIdentical(target1, current)
    }

    x0 <- c(0, NA, NaN, 0, NA, Inf, 0, NA, -Inf, 0, Inf, -Inf)
    x <- Rle(x0)
    for (k in 1:2) {
        target1 <- .naive_runmean(x0, k, na.rm=TRUE)
        target2 <- .naive_runmean(x, k, na.rm=TRUE)
        checkIdentical(target1, target2)
        current <- as.vector(runmean(x, k, na.rm=TRUE))
        checkIdentical(target1, current) 
 
        target1 <- .naive_runmean(x0, k, na.rm=FALSE)
        target2 <- .naive_runmean(x, k, na.rm=FALSE)
        checkIdentical(target1, target2)
        #current <- as.vector(runmean(x, k, na.rm=FALSE))
        #checkIdentical(target1, current)
    } 
}

.naive_runwtsum <- function(x, k, wt, na.rm=FALSE)
    sapply(0:(length(x)-k),
        function(offset) {
            xwt <- x[1:k + offset] * wt 
            sum(xwt, na.rm=na.rm)}) 

test_Rle_runwtsum_real <- function() {

    x0 <- c(NA, NaN, Inf, -Inf) 
    x <- Rle(x0)
    wt <- rep(1, 4)
    target1 <- .naive_runwtsum(x0, 4, wt, na.rm=TRUE)
    target2 <- .naive_runwtsum(x, 4, wt, na.rm=TRUE)
    checkIdentical(target1, target2) 
    current <- as.vector(runwtsum(x, 4, wt, na.rm=TRUE))
    checkIdentical(target1, current)
    target1 <- .naive_runwtsum(x0, 4, wt, na.rm=FALSE)
    target2 <- .naive_runwtsum(x, 4, wt, na.rm=FALSE)
    checkIdentical(target1, target2) 
    #current <- as.vector(runwtsum(x, 4, wt, na.rm=FALSE))
    #checkIdentical(target1, current) 

    x0 <- c(NA, Inf, NA, -Inf, Inf, -Inf, NaN, Inf, NaN, -Inf)
    x <- Rle(x0)
    for (k in 1:2) {
        if (k==1)
            wt <- 1
        else
            wt <- c(1, 1) 
        target1 <- .naive_runwtsum(x0, k, wt, na.rm=TRUE)
        target2 <- .naive_runwtsum(x, k, wt, na.rm=TRUE)
        checkIdentical(target1, target2)
        current <- as.vector(runwtsum(x, k, wt, na.rm=TRUE))
        checkIdentical(target1, current) 

        target1 <- .naive_runwtsum(x0, k, wt, na.rm=FALSE)
        target2 <- .naive_runwtsum(x, k, wt, na.rm=FALSE)
        checkIdentical(target1, target2)
        current <- as.vector(runwtsum(x, k, wt, na.rm=FALSE))
        checkIdentical(target1, current)
    }
 
    x0 <- c(1, NA, 1, NaN, 1, NA)
    x <- Rle(x0)
    for (k in 1:2) {
        if (k==1)
            wt <- 2 
        else
            wt <- c(1, 1)
        target1 <- .naive_runwtsum(x0, k, wt, na.rm=FALSE)
        target2 <- .naive_runwtsum(x, k, wt, na.rm=FALSE)
        checkIdentical(target1, target2)
        current <- as.vector(runwtsum(x, k, wt, na.rm=FALSE))
        checkIdentical(target1, current)
    }
}

test_Rle_runwtsum_integer <- function() {

    x0 <- c(NA_integer_, 1L, 1L)
    x <- Rle(x0)
    iwt <- rep(2L, 3)
    for (k in 1:3) {
        wt <- iwt[1:k]
        target1 <- .naive_runwtsum(x0, k, wt, na.rm=TRUE)
        target2 <- .naive_runwtsum(x, k, wt, na.rm=TRUE)
        checkIdentical(target1, target2) 
        current <- as.vector(runwtsum(x, k, wt, na.rm=TRUE))
        checkIdentical(as.numeric(target1), current)

        target1 <- .naive_runwtsum(x0, k, wt, na.rm=FALSE)
        target2 <- .naive_runwtsum(x, k, wt, na.rm=FALSE)
        checkIdentical(target1, target2) 
        current <- as.vector(runwtsum(x, k, wt, na.rm=FALSE))
        checkIdentical(as.numeric(target1), current)
    }

    x0 <- c(1L, NA_integer_, 1L)
    x <- Rle(x0)
    iwt <- rep(2L, 3)
    for (k in 1:3) {
        wt <- iwt[1:k]
        target1 <- .naive_runwtsum(x0, k, wt, na.rm=TRUE)
        target2 <- .naive_runwtsum(x, k, wt, na.rm=TRUE)
        checkIdentical(target1, target2) 
        current <- as.vector(runwtsum(x, k, wt, na.rm=TRUE))
        checkIdentical(as.numeric(target1), current)

        target1 <- .naive_runwtsum(x0, k, wt, na.rm=FALSE)
        target2 <- .naive_runwtsum(x, k, wt, na.rm=FALSE)
        checkIdentical(target1, target2) 
        current <- as.vector(runwtsum(x, k, wt, na.rm=FALSE))
        checkIdentical(as.numeric(target1), current)
    }
}

.naive_runq <- function(x, k, i, na.rm=FALSE)
    sapply(0:(length(x)-k),
        function(offset) {
            xsub <- x[1:k + offset]
            if (!na.rm) { 
                ## Manually handle NA's because they are not allowed
                ## in 'x' of quantile(x, ...) when na.rm=FALSE.
                if (any(is.na(xsub)))
                    NA 
                else
                    quantile(xsub, probs=i/k, na.rm=na.rm, names=FALSE, type=3)
            } else {
                ## If all NA's, just return first NA value.
                ## Not handled in quantile().
                if (all(is.na(xsub))) {
                    xsub[1]
                } else {
                    xsub <- xsub[!is.na(xsub)]
                    quantile(xsub, probs=i/k, na.rm=na.rm, names=FALSE, type=3)
                }
            }
        }, USE.NAMES=FALSE)

test_Rle_runq_real <- function() {

    x0 <- c(NA_real_)
    x <- Rle(x0)
    k <- length(x); i <- 1
    target1 <- as.numeric(.naive_runq(x0, k, i, na.rm=TRUE))
    current <- as.numeric(runq(x, k, i, na.rm=TRUE))
    checkIdentical(target1, current)

    x0 <- c(3, NA, 1, NaN, 4, Inf, 2, -Inf)
    x <- Rle(x0)
    k <- length(x)
    for (i in c(1, length(x))) {
        target1 <- as.numeric(.naive_runq(x0, k, i, na.rm=TRUE))
        current <- as.numeric(runq(x, k, i, na.rm=TRUE))
        checkIdentical(target1, current)
 
        target1 <- as.numeric(.naive_runq(x0, k, i, na.rm=FALSE))
        current <- as.numeric(runq(x, k, i, na.rm=FALSE))
        checkIdentical(target1, current)
    }

    x0 <- c(3, NA, 1, NaN, 4, Inf, 2, -Inf)
    x <- Rle(x0)
    i <- 1 
    ## NOTE : special case k=1, returns NA not NaN
    target1 <- c(3, NA, 1, NA, 4, Inf, 2, -Inf)
    current <- as.numeric(runq(x, k=1, i=1, na.rm=TRUE))
    checkIdentical(target1, current)
    for (k in c(2:length(x))) {
        target1 <- as.numeric(.naive_runq(x0, k, i, na.rm=TRUE))
        current <- as.numeric(runq(x, k, i, na.rm=TRUE))
        checkIdentical(target1, current)

        target1 <- as.numeric(.naive_runq(x0, k, i, na.rm=FALSE))
        current <- as.numeric(runq(x, k, i, na.rm=FALSE))
        checkIdentical(target1, current)
    }

    x0 <- c(1, 2, 3, 4, 5)
    x <- Rle(x0)
    k <- length(x); i <- 4 
    target1 <- .naive_runq(x0, k, i, na.rm=TRUE)
    current <- as.vector(runq(x, k, i, na.rm=TRUE))
    checkIdentical(target1, current)

    x0 <- c(1, 2, 3, NA, NA)
    x <- Rle(x0)
    k <- length(x); i <- 4 
    target1 <- .naive_runq(x0, k, i, na.rm=TRUE)
    current <- as.vector(runq(x, k, i, na.rm=TRUE))
    checkIdentical(target1, current)
}

test_Rle_runq_integer <- function() {

    x0 <- c(NA_integer_)
    x <- Rle(x0)
    k <- length(x); i <- 1
    target1 <- as.numeric(.naive_runq(x0, k, i, na.rm=TRUE))
    current <- as.numeric(runq(x, k, i, na.rm=TRUE))
    checkIdentical(target1, current)

    x0 <- NA_integer_
    x <- Rle(x0)
    k <- i <- 1 
    target1 <- unlist(.naive_runq(x0, k, i, na.rm=TRUE))
    target2 <- as.vector(do.call(c, (.naive_runq(x, k, i, na.rm=TRUE))))
    checkIdentical(target1, target2) 
    current <- as.vector(runq(x, k, i, na.rm=TRUE))
    checkIdentical(as.integer(unname(target1)), current)

    x0 <- c(NA_integer_, 2L, 1L)
    x <- Rle(x0)
    k <- 3 
    for (i in 1:3) {
        target1 <- unlist(.naive_runq(x0, k, i, na.rm=TRUE))
        current <- as.vector(runq(x, k, i, na.rm=TRUE))
        checkIdentical(unname(target1), current)

        target1 <- unlist(.naive_runq(x0, k, i, na.rm=FALSE))
        current <- as.integer(runq(x, k, i, na.rm=FALSE))
        checkIdentical(as.integer(target1), current)
    }

    x0 <- c(3L, 2L, NA_integer_, NA_integer_, 1L, 2L)
    x <- Rle(x0)
    i <- 1
    for (k in 1:6) {
        target1 <- unlist(.naive_runq(x0, k, i, na.rm=TRUE))
        current <- as.vector(runq(x, k, i, na.rm=TRUE))
        checkIdentical(target1, current)

        target1 <- unlist(.naive_runq(x0, k, i, na.rm=FALSE))
        current <- as.integer(runq(x, k, i, na.rm=FALSE))
        checkIdentical(as.integer(target1), current)
    }
}

