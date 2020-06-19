
test_ellipsis_forwarding_for_mapply <- function()
{
    mapply_wrapper <- function(FUN, x, ...) mapply(FUN, x, ...)
    x <- list(a=1:3, 1:2)
    y <- list(104:105, B=103)

    target <- mapply(append, x, y)
    checkIdentical(target, mapply_wrapper(append, x, y))

    MoreArgs <- list(after=0)
    target <- mapply(append, x, y, MoreArgs=MoreArgs)
    current <- mapply_wrapper(append, x, y, MoreArgs=MoreArgs)
    checkIdentical(target, current)

    MoreArgs <- list(after=2)
    target <- mapply(append, x, y, MoreArgs=MoreArgs, USE.NAMES=FALSE)
    current <- mapply_wrapper(append, x, y, MoreArgs=MoreArgs, USE.NAMES=FALSE)
    checkIdentical(target, current)
}

