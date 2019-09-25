
test_ellipsis_forwarding_for_order <- function()
{
    x <- list(c(NA,11:13), c(21:22,NA))

    target <- lapply(x, base::order)
    checkIdentical(target, lapply(x, order))

    target <- lapply(x, base::order, na.last=TRUE)
    checkIdentical(target, lapply(x, order, na.last=TRUE))

    target <- lapply(x, base::order, na.last=FALSE)
    checkIdentical(target, lapply(x, order, na.last=FALSE))

    target <- lapply(x, base::order, na.last=NA)
    checkIdentical(target, lapply(x, order, na.last=NA))
}

