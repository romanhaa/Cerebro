
test_ellipsis_forwarding_for_Extremes <- function()
{
    for (FUN in c("pmax", "pmin", "pmax.int", "pmin.int")) {
        FUN <- match.fun(FUN)
        FUN_wrapper <- function(x, ...) FUN(x, ...)
        x <- c(1:3, NA)
        y <- c(NA, 3:1)
        checkIdentical(FUN(x, y), FUN_wrapper(x, y))
        checkIdentical(FUN(x, y, na.rm=FALSE), FUN_wrapper(x, y, na.rm=FALSE))
        checkIdentical(FUN(x, y, na.rm=TRUE), FUN_wrapper(x, y, na.rm=TRUE))
    }
}

