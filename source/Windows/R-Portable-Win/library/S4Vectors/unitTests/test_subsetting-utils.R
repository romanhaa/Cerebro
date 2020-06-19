.NAMES0 <- c("C", "AA", "BB", "A", "", "A", "AA", "BB", "DD")

test_normalizeDoubleBracketSubscript <- function()
{
    ## These "core tests" don't even look at 'x'.
    do_core_tests <- function(x, exact=TRUE) {
        for (i in list(TRUE, FALSE, 1i, as.raw(1),
                       integer(0), 1:3, character(0), c("A", "b"))) {
            checkException(normalizeDoubleBracketSubscript(i, x,
                                                           exact=exact))
            checkException(normalizeDoubleBracketSubscript(Rle(i), x,
                                                           exact=exact))
        }

        for (i in list(NA, NA_integer_, NA_real_, NA_character_, NA_complex_)) {
            checkException(normalizeDoubleBracketSubscript(i, x, exact=exact))
            current <- normalizeDoubleBracketSubscript(i, x, exact=exact,
                                                       allow.NA=TRUE)
            checkIdentical(NA, current)
            checkException(normalizeDoubleBracketSubscript(Rle(i), x,
                                                           exact=exact))
            current <- normalizeDoubleBracketSubscript(Rle(i), x, exact=exact,
                                                       allow.NA=TRUE)
            checkIdentical(NA, current)
        }

        ## Error: [[ subscript must be >= 1
        for (i in list(0L, 0.99, -1)) {
            checkException(normalizeDoubleBracketSubscript(i, x,
                                                           exact=exact))
            checkException(normalizeDoubleBracketSubscript(Rle(i), x,
                                                           exact=exact))
            checkException(normalizeDoubleBracketSubscript(i, x,
                                                           exact=exact,
                                                           allow.append=TRUE))
            checkException(normalizeDoubleBracketSubscript(Rle(i), x,
                                                           exact=exact,
                                                           allow.append=TRUE))
        }
    }

    test_invalid_position <- function(i, x, allow.append=FALSE) {
        for (exact in list(TRUE, FALSE)) {
            for (allow.NA in list(FALSE, TRUE)) {
                for (allow.nomatch in list(FALSE, TRUE)) {
                    checkException(normalizeDoubleBracketSubscript(i, x,
                                              exact=exact,
                                              allow.append=allow.append,
                                              allow.NA=allow.NA,
                                              allow.nomatch=allow.nomatch))
                    checkException(normalizeDoubleBracketSubscript(Rle(i), x,
                                              exact=exact,
                                              allow.append=allow.append,
                                              allow.NA=allow.NA,
                                              allow.nomatch=allow.nomatch))
                }
            }
        }
    }

    test_valid_position <- function(i, x, target, allow.append=FALSE) {
        for (exact in list(TRUE, FALSE)) {
            for (allow.NA in list(FALSE, TRUE)) {
                for (allow.nomatch in list(FALSE, TRUE)) {
                    current <- normalizeDoubleBracketSubscript(i, x,
                                              exact=exact,
                                              allow.append=allow.append,
                                              allow.NA=allow.NA,
                                              allow.nomatch=allow.nomatch)
                    checkIdentical(target, current)
                    current <- normalizeDoubleBracketSubscript(Rle(i), x,
                                              exact=exact,
                                              allow.append=allow.append,
                                              allow.NA=allow.NA,
                                              allow.nomatch=allow.nomatch)
                    checkIdentical(target, current)
                }
            }
        }
    }

    test_invalid_name <- function(name, x, exact=TRUE) {
        for (i in list(name, Rle(name), factor(name), Rle(factor(name)))) {
            for (allow.append in list(FALSE, TRUE)) {
                for (allow.NA in list(FALSE, TRUE)) {
                    checkException(normalizeDoubleBracketSubscript(i, x,
                                              exact=exact,
                                              allow.append=allow.append,
                                              allow.NA=allow.NA))
                    checkException(normalizeDoubleBracketSubscript(i, x,
                                              exact=exact,
                                              allow.append=allow.append,
                                              allow.NA=allow.NA,
                                              allow.nomatch=FALSE))
                    current <- normalizeDoubleBracketSubscript(i, x,
                                              exact=exact,
                                              allow.append=allow.append,
                                              allow.NA=allow.NA,
                                              allow.nomatch=TRUE)
                    checkIdentical(NA, current)
                }
            }
        }
    }

    test_valid_name <- function(name, x, target, exact=TRUE) {
        for (i in list(name, Rle(name), factor(name), Rle(factor(name)))) {
            for (allow.append in list(FALSE, TRUE)) {
                for (allow.NA in list(FALSE, TRUE)) {
                    for (allow.nomatch in list(FALSE, TRUE)) {
                        current <- normalizeDoubleBracketSubscript(i, x,
                                              exact=exact,
                                              allow.append=allow.append,
                                              allow.NA=allow.NA,
                                              allow.nomatch=allow.nomatch)
                        checkIdentical(target, current)
                    }
                }
            }
        }
    }

    ## ----------------------------------------------------------------- ##

    do_basic_tests_on_empty_object <- function(x) {
        do_core_tests(x, exact=TRUE)
        do_core_tests(x, exact=FALSE)

        ## (1) With a single non-NA number.

        ## Error: subscript is out of bounds
        test_invalid_position(1L, x, allow.append=FALSE)
        test_invalid_position(1, x, allow.append=FALSE)

        test_valid_position(1L, x, 1L, allow.append=TRUE)
        test_valid_position(1.99, x, 1L, allow.append=TRUE)

        ## Error: [[ subscript must be <= length(x) + 1
        test_invalid_position(2L, x, allow.append=TRUE)
        test_invalid_position(2, x, allow.append=TRUE)

        ## (2) With a single non-NA string.

        test_invalid_name("A", x, exact=TRUE)
        test_invalid_name("A", x, exact=FALSE)
    }

    x <- list()
    do_basic_tests_on_empty_object(x)

    ## ----------------------------------------------------------------- ##

    names(x) <- character(0)
    do_basic_tests_on_empty_object(x)

    ## ----------------------------------------------------------------- ##

    do_basic_tests_on_full_object <- function(x) {
        do_core_tests(x, exact=TRUE)
        do_core_tests(x, exact=FALSE)

        ## (1) With a single non-NA number.

        test_valid_position(1L, x, 1L, allow.append=FALSE)
        test_valid_position(1L, x, 1L, allow.append=TRUE)
        test_valid_position(1.99, x, 1L, allow.append=FALSE)
        test_valid_position(1.99, x, 1L, allow.append=TRUE)

        test_valid_position(9L, x, 9L, allow.append=FALSE)
        test_valid_position(9L, x, 9L, allow.append=TRUE)
        test_valid_position(9.99, x, 9L, allow.append=FALSE)
        test_valid_position(9.99, x, 9L, allow.append=TRUE)

        ## Error: subscript is out of bounds
        test_invalid_position(10L, x, allow.append=FALSE)
        test_invalid_position(10.99, x, allow.append=FALSE)

        test_valid_position(10L, x, 10L, allow.append=TRUE)
        test_valid_position(10.99, x, 10L, allow.append=TRUE)

        ## Error: [[ subscript must be <= length(x) + 1
        test_invalid_position(11L, x, allow.append=TRUE)
        test_invalid_position(11, x, allow.append=TRUE)
    }

    x <- as.list(letters[1:9])
    do_basic_tests_on_full_object(x)

    ## (2) With a single non-NA string.

    test_invalid_name("A", x, exact=TRUE)
    test_invalid_name("A", x, exact=FALSE)

    ## ----------------------------------------------------------------- ##

    names(x) <- .NAMES0
    do_basic_tests_on_full_object(x)

    ## (2) With a single non-NA string.

    ## Exact matching.

    test_invalid_name("Z", x, exact=TRUE)
    test_invalid_name("B", x, exact=TRUE)
    test_invalid_name("D", x, exact=TRUE)

    test_valid_name("C", x, 1L, exact=TRUE)
    test_valid_name("BB", x, 3L, exact=TRUE)
    test_valid_name("A", x, 4L, exact=TRUE)
    test_valid_name("AA", x, 2L, exact=TRUE)
    test_valid_name("DD", x, 9L, exact=TRUE)

    ## Partial matching.

    test_invalid_name("Z", x, exact=FALSE)
    test_invalid_name("B", x, exact=FALSE)  # ambiguous partial matching

    test_valid_name("C", x, 1L, exact=FALSE)
    test_valid_name("BB", x, 3L, exact=FALSE)
    test_valid_name("A", x, 4L, exact=FALSE)
    test_valid_name("AA", x, 2L, exact=FALSE)
    test_valid_name("DD", x, 9L, exact=FALSE)
    test_valid_name("D", x, 9L, exact=FALSE)
}

.do_test_getListElement_list_or_data.frame <- function(x0)
{
    ## These "core tests" don't even look at 'x'.
    do_core_tests <- function(x, exact=TRUE) {
        for (i in list(TRUE, FALSE, 1i, as.raw(1),
                       integer(0), 1:3, character(0), c("A", "b"))) {
            checkException(getListElement(x, i, exact=exact))
            checkException(getListElement(x, Rle(i), exact=exact))
        }

        for (i in list(NA, NA_integer_, NA_real_, NA_character_, NA_complex_)) {
            current <- getListElement(x, i, exact=exact)
            checkIdentical(NULL, current)
            current <- getListElement(x, Rle(i), exact=exact)
            checkIdentical(NULL, current)
        }

        ## Error: [[ subscript must be >= 1
        for (i in list(0L, 0.99, -1)) {
            checkException(getListElement(x, i, exact=exact))
            checkException(getListElement(x, Rle(i), exact=exact))
        }
    }

    test_invalid_position <- function(x, i) {
        for (exact in list(TRUE, FALSE)) {
            checkException(getListElement(x, i, exact=exact))
            checkException(getListElement(x, Rle(i), exact=exact))
        }
    }

    test_valid_position <- function(x, i) {
        target <- `[[`(x, i)
        for (exact in list(TRUE, FALSE)) {
            current <- getListElement(x, i, exact=exact)
            checkIdentical(target, current)
            current <- getListElement(x, Rle(i), exact=exact)
            checkIdentical(target, current)
        }
    }

    test_valid_name <- function(x, name, exact=TRUE) {
        target <- `[[`(x, name, exact=exact)
        for (i in list(name, Rle(name), factor(name), Rle(factor(name)))) {
            current <- getListElement(x, i, exact=exact)
            checkIdentical(target, current)
        }
    }

    ## ----------------------------------------------------------------- ##

    stopifnot(identical(names(x0), .NAMES0))

    do_basic_tests_on_empty_object <- function(x) {
        do_core_tests(x, exact=TRUE)
        do_core_tests(x, exact=FALSE)

        ## (1) With a single non-NA number.

        ## Error: subscript is out of bounds
        test_invalid_position(x, 1L)
        test_invalid_position(x, 1)

        ## (2) With a single non-NA string.

        ## No match
        test_valid_name(x, "A", exact=TRUE)
        test_valid_name(x, "A", exact=FALSE)
    }

    if (!(is.data.frame(x0) || is(x0, "DataFrame"))) {
        ## Test on empty unnamed object.
        x <- x0[0]
        names(x) <- NULL
        do_basic_tests_on_empty_object(x)
    }

    ## ----------------------------------------------------------------- ##

    ## Test on empty named object.
    x <- x0[0]
    do_basic_tests_on_empty_object(x)

    ## ----------------------------------------------------------------- ##

    do_basic_tests_on_full_object <- function(x) {
        do_core_tests(x, exact=TRUE)
        do_core_tests(x, exact=FALSE)

        ## (1) With a single non-NA number.

        test_valid_position(x, 1L)
        test_valid_position(x, 1.99)

        test_valid_position(x, 9L)
        test_valid_position(x, 9.99)

        test_invalid_position(x, 10L)
        test_invalid_position(x, 10)
        test_invalid_position(x, 10.99)
    }

    if (!(is.data.frame(x0) || is(x0, "DataFrame"))) {
        ## Test on full unnamed object.
        x <- x0
        names(x) <- NULL
        do_basic_tests_on_full_object(x)

        ## (2) With a single non-NA string.

        ## No match
        test_valid_name(x, "A", exact=TRUE)
        test_valid_name(x, "A", exact=FALSE)
    }

    ## ----------------------------------------------------------------- ##

    ## Test on full named object.
    x <- x0
    do_basic_tests_on_full_object(x)

    ## (2) With a single non-NA string.

    ## Exact matching.

    ## No match
    test_valid_name(x, "Z", exact=TRUE)
    test_valid_name(x, "B", exact=TRUE)
    test_valid_name(x, "D", exact=TRUE)

    ## Match
    test_valid_name(x, "C", exact=TRUE)
    test_valid_name(x, "BB", exact=TRUE)
    test_valid_name(x, "A", exact=TRUE)
    test_valid_name(x, "AA", exact=TRUE)
    test_valid_name(x, "DD", exact=TRUE)

    ## Partial matching.

    ## No match
    test_valid_name(x, "Z", exact=FALSE)
    test_valid_name(x, "B", exact=FALSE)  # ambiguous partial matching

    ## Match
    test_valid_name(x, "C", exact=FALSE)
    test_valid_name(x, "BB", exact=FALSE)
    test_valid_name(x, "A", exact=FALSE)
    test_valid_name(x, "AA", exact=FALSE)
    test_valid_name(x, "DD", exact=FALSE)
    test_valid_name(x, "D", exact=FALSE)
}

test_getListElement_list <- function()
{
    x <- setNames(as.list(letters[1:9]), .NAMES0)
    .do_test_getListElement_list_or_data.frame(x)
    x <- as.data.frame(lapply(1:9, function(i) {10L*i + 1:4} ))
    colnames(x) <- .NAMES0
    .do_test_getListElement_list_or_data.frame(x)
}

.do_test_setListElement_list_or_data.frame <- function(x0, value0)
{
    ## These "core tests" don't even look at 'x' or 'value'.
    do_core_tests <- function(x, value) {
        for (i in list(TRUE, FALSE, 1i, as.raw(1),
                       integer(0), 1:3, character(0), c("A", "b"))) {
            checkException(setListElement(x, i, value))
            checkException(setListElement(x, Rle(i), value))
        }

        for (i in list(NA, NA_integer_, NA_real_, NA_character_, NA_complex_)) {
            checkException(setListElement(x, i, value))
            checkException(setListElement(x, Rle(i), value))
        }

        ## Error: [[ subscript must be >= 1
        for (i in list(0L, 0.99, -1)) {
            checkException(setListElement(x, i, value))
            checkException(setListElement(x, Rle(i), value))
        }
    }

    ## Does not look at 'value'.
    test_invalid_position <- function(x, i, value) {
        checkException(setListElement(x, i, value))
        checkException(setListElement(x, Rle(i), value))
    }

    test_valid_position <- function(x, i, value) {
        target <- `[[<-`(x, i, value=value)
        ## `[[<-.data.frame` does some terrible mangling of the colnames when
        ## appending a column to 'x' if 'colnames(x)' contains duplicates.
        ## We fix this.
        if (is.data.frame(x) && ncol(target) > ncol(x))
            colnames(target) <- c(colnames(x), "")
        current <- setListElement(x, i, value)
        checkIdentical(target, current)
        current <- setListElement(x, Rle(i), value)
        checkIdentical(target, current)
    }

    test_valid_name <- function(x, name, value) {
        target <- `[[<-`(x, name, value=value)
        ## `[[<-.data.frame` does some terrible mangling of the colnames when
        ## appending a column to 'x' if 'colnames(x)' contains duplicates.
        ## We fix this.
        if (is.data.frame(x) && ncol(target) > ncol(x))
            colnames(target) <- c(colnames(x), name)
        for (i in list(name, Rle(name), factor(name), Rle(factor(name)))) {
            current <- setListElement(x, i, value)
            checkIdentical(target, current)
        }
    }

    ## ----------------------------------------------------------------- ##

    stopifnot(identical(names(x0), .NAMES0))

    do_basic_tests_on_empty_object <- function(x) {
        do_core_tests(x, NULL)
        do_core_tests(x, value0)

        ## (1) With a single non-NA number.

        ## No-op
        test_valid_position(x, 1L, NULL)
        test_valid_position(x, 1, NULL)
        test_valid_position(x, 1.99, NULL)

        ## Append naked 'value0' to 'x'.
        test_valid_position(x, 1L, value0)
        test_valid_position(x, 1, value0)
        test_valid_position(x, 1.99, value0)

        ## Error: [[ subscript must be <= length(x) + 1
        test_invalid_position(x, 2L, NULL)
        test_invalid_position(x, 2, value0)

        ## (2) With a single non-NA string.

        ## No match
        test_valid_name(x, "A", NULL)  # no-op
        test_valid_name(x, "A", value0)  # append
    }

    if (!(is.data.frame(x0) || is(x0, "DataFrame"))) {
        ## Test on empty unnamed object.
        x <- x0[0]
        names(x) <- NULL
        do_basic_tests_on_empty_object(x)
    }

    ## ----------------------------------------------------------------- ##

    ## Test on empty named object.
    x <- x0[0]
    do_basic_tests_on_empty_object(x)

    ## ----------------------------------------------------------------- ##

    do_basic_tests_on_full_object <- function(x) {
        do_core_tests(x, NULL)
        do_core_tests(x, value0)

        ## (1) With a single non-NA number.

        ## Remove 1st list element
        test_valid_position(x, 1L, NULL)
        test_valid_position(x, 1.99, NULL)

        ## Replace 1st list element
        test_valid_position(x, 1L, value0)
        test_valid_position(x, 1.99, value0)

        ## Remove last list element
        test_valid_position(x, 9L, NULL)
        test_valid_position(x, 9.99, NULL)

        ## Replace last list element
        test_valid_position(x, 9L, value0)
        test_valid_position(x, 9.99, value0)

        ## No-op
        test_valid_position(x, 10L, NULL)
        test_valid_position(x, 10, NULL)
        test_valid_position(x, 10.99, NULL)

        ## Append naked 'value0' to 'x'
        test_valid_position(x, 10L, value0)
        test_valid_position(x, 10, value0)
        test_valid_position(x, 10.99, value0)

        ## Error: [[ subscript must be <= length(x) + 1
        test_invalid_position(x, 11L, NULL)
        test_invalid_position(x, 11, value0)
    }

    if (!(is.data.frame(x0) || is(x0, "DataFrame"))) {
        ## Test on full unnamed object.
        x <- x0
        names(x) <- NULL
        do_basic_tests_on_full_object(x)

        ## (2) With a single non-NA string.

        ## No match
        test_valid_name(x, "A", NULL)  # no-op
        test_valid_name(x, "A", value0)  # append
    }

    ## ----------------------------------------------------------------- ##

    ## Test on full named object.
    x <- x0
    do_basic_tests_on_full_object(x)

    ## (2) With a single non-NA string.

    ## No match.

    ## No-op
    test_valid_name(x, "Z", NULL)
    test_valid_name(x, "B", NULL)
    test_valid_name(x, "D", NULL)

    ## Append named 'value0' to 'x'
    test_valid_name(x, "Z", value0)
    test_valid_name(x, "B", value0)
    test_valid_name(x, "D", value0)

    ## Match.

    ## Remove named list element
    test_valid_name(x, "C", NULL)
    test_valid_name(x, "BB", NULL)
    test_valid_name(x, "A", NULL)
    test_valid_name(x, "AA", NULL)
    test_valid_name(x, "DD", NULL)

    ## Replace named list element
    test_valid_name(x, "C", value0)
    test_valid_name(x, "BB", value0)
    test_valid_name(x, "A", value0)
    test_valid_name(x, "AA", value0)
    test_valid_name(x, "DD", value0)
}

test_setListElement_list <- function()
{
    x <- setNames(as.list(letters[1:9]), .NAMES0)
    .do_test_setListElement_list_or_data.frame(x, 9:6)
    x <- as.data.frame(lapply(1:9, function(i) {10L*i + 1:4} ))
    colnames(x) <- .NAMES0
    .do_test_setListElement_list_or_data.frame(x, 9:6)
    .do_test_setListElement_list_or_data.frame(x, letters[1:4])
}

