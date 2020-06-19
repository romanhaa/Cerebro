test_stopf_warningf <- function() {
    .stopf <- GSEABase:::.stopf
    .warningf <- GSEABase:::.warningf

    current <- tryCatch(.stopf(letters[1]), error=conditionMessage)
    checkIdentical("a", current)

    current <- tryCatch(.warningf(letters[1]), warning=conditionMessage)
    checkIdentical("a", current)

    ## support strings of length more than 1
    current <- tryCatch(.stopf(letters[1:2]), error=conditionMessage)
    checkIdentical("a\nb", current)

    current <- tryCatch(.warningf(letters[1:2]), warning=conditionMessage)
    checkIdentical("a\nb", current)
}
