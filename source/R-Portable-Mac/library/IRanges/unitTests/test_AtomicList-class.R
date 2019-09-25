test_AtomicList_constructors <- function() {
    subclasses <- c(logical="LogicalList",
                    integer="IntegerList",
                    #double="NumericList",
                    numeric="NumericList",
                    complex="ComplexList",
                    character="CharacterList",
                    raw="RawList",
                    Rle="RleList")
    for (elt_type in names(subclasses)) {
        subclass <- subclasses[[elt_type]]
        constructor <- get(subclass)
        vec1 <- get(elt_type)(6)
        vec2 <- get(elt_type)(8)
        target <- list(A=vec1, B=vec2)
        for (compress in c(TRUE, FALSE)) {
            current <- constructor(A=vec1, B=vec2, compress=compress)
            checkTrue(is(current, subclass))
            checkIdentical(compress, is(current, "CompressedList"))
            checkIdentical(elt_type, elementType(current))
            checkIdentical(target, as.list(current))
            checkIdentical(unname(target), as.list(current, use.names=FALSE))
        }
    }
}

test_AtomicList_general <- function() {
    vec1 <- c(1L,2L,NA,3L,NA,5L,2L,8L)
    vec2 <- c(NA,15L,45L,20L,NA,1L,15L,100L,80L,5L,NA)
    for (compress in c(TRUE, FALSE)) {
        for (type in c("IntegerList", "RleList")) {
            list1 <- do.call(type, list(one = vec1, vec2, compress = compress))
            checkIdentical(as.list(list1 %in% c(1L, 5L)),
                           lapply(list1, "%in%", c(1L, 5L)))
            checkIdentical(lapply(list1 %in%
                                  IntegerList(one = vec1, vec2,
                                              compress = compress),
                                  as.vector),
                           mapply("%in%", lapply(list1, as.vector),
                                  list(one = vec1, vec2)))
            checkIdentical(as.list(is.na(list1)), lapply(list1, is.na))
            checkIdentical(as.list(match(list1, c(1L, 5L))),
                           lapply(list1, match, c(1L, 5L)))
            checkIdentical(lapply(match(list1,
                                        IntegerList(one = vec1, vec2,
                                                    compress = compress)),
                                  as.vector),
                           mapply(match, lapply(list1, as.vector),
                                  list(one = vec1, vec2)))
            checkIdentical(as.list(sort(list1)), lapply(list1, sort))
            checkIdentical(as.list(unique(list1)), lapply(list1, unique))
        }
    }
}

test_RleList_methods <- function() {
    x1 <- RleList(11:15, 15L, integer(0), 15:16, compress=FALSE)
    x2 <- RleList(11:15, 15L, integer(0), 15:16, compress=TRUE)
    checkIdentical(as(runValue(x1), "CompressedIntegerList"), runValue(x2))
    checkIdentical(as(runLength(x1), "CompressedIntegerList"), runLength(x2))
    checkIdentical(as(ranges(x1), "CompressedIRangesList"), ranges(x2))

    a1 <- Rle(1, 999722111)
    a2 <- 20 * a1
    a <- RleList(a1, a2, compress=TRUE)
    b1 <- c(a1, a1)
    b2 <- 20 * b1
    b <- RleList(b1, b2, compress=FALSE)
    ## FIXME: 'a1 <= 19:21' is taking forever and eats up all the memory in
    ## BioC <= 2.12! Seems like 'a1' is expanded to integer vector first, which
    ## is not good :-/
    #for (y in list(8L, 8, 19:21)) {
    for (y in list(8L, 8)) {
        ## With a CompressedRleList
        target <- RleList(a1 <= y, a2 <= y, compress=TRUE)
        current <- a <= y
        checkIdentical(target, current)
        ## With a SimpleRleList
        target <- RleList(b1 <= y, b2 <= y, compress=FALSE)
        current <- b <= y
        checkIdentical(target, current)
    }
}

