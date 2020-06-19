
test_ellipsis_forwarding_for_paste <- function()
{
    x <- list(letters, LETTERS)

    target <- sapply(x, base::paste)
    checkIdentical(target, sapply(x, paste))

    target <- sapply(x, base::paste, collapse="")
    checkIdentical(target, sapply(x, paste, collapse=""))
}

