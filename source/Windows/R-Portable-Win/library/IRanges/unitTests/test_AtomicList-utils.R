test_AtomicList_GroupGenerics <- function() {
    vec1 <- c(1L,2L,3L,5L,2L,8L)
    vec2 <- c(15L,45L,20L,1L,15L,100L,80L,5L)
    for (compress in c(TRUE, FALSE)) {
        for (type in c("IntegerList", "RleList")) {
            list1 <- do.call(type, list(one = vec1, vec2, compress = compress))
            checkIdentical(as.list(list1 + list1), Map("+", list1, list1))
            checkIdentical(as.list(log(list1)), lapply(list1, log))
            checkIdentical(as.list(round(sqrt(list1))),
                           lapply(list1, function(x) round(sqrt(x))))
            checkIdentical(sum(list1), sapply(list1, sum))
        }
    }
}

test_AtomicList_logical <- function() {
    vec1 <- c(TRUE,NA,FALSE, NA)
    vec2 <- c(TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE)
    for (compress in c(TRUE, FALSE)) {
        for (type in c("LogicalList", "RleList")) {
            list1 <- do.call(type, list(one = vec1, vec2, compress = compress))
            checkIdentical(as.list(!list1), lapply(list1, "!"))
            checkIdentical(as.list(which(list1)), lapply(list1, which))
        }
    }
}

test_AtomicList_numerical <- function() {
    vec1 <- c(1L,2L,NA,3L,NA,5L,2L,8L)
    vec2 <- c(NA,15L,45L,20L,NA,1L,15L,100L,80L,5L,NA)
    for (compress in c(TRUE, FALSE)) {
        for (type in c("IntegerList", "RleList")) {
            list1 <- do.call(type, list(one = vec1, vec2, compress = compress))
            list2 <- endoapply(list1, rev)
            checkIdentical(as.list(diff(list1)), lapply(list1, diff))
            checkIdentical(as.list(pmax(list1, list2)),
                           mapply(pmax, list1, list2))
            checkIdentical(as.list(pmin(list1, list2)),
                           mapply(pmin, list1, list2))
            checkIdentical(as.list(pmax.int(list1, list2)),
                           mapply(pmax.int, list1, list2))
            checkIdentical(as.list(pmin.int(list1, list2)),
                           mapply(pmin.int, list1, list2))
            checkIdentical(mean(list1, na.rm=TRUE),
                           sapply(list1, mean, na.rm=TRUE))
            checkIdentical(var(list1, na.rm=TRUE),
                           sapply(list1, var, na.rm=TRUE))
            checkIdentical(cov(list1, list2, use="complete.obs"),
                           mapply(cov, list1, list2,
                                  MoreArgs = list(use="complete.obs")))
            checkIdentical(cor(list1, list2, use="complete.obs"),
                           mapply(cor, list1, list2,
                                  MoreArgs = list(use="complete.obs")))
            checkIdentical(sd(list1, na.rm=TRUE),
                           sapply(list1, sd, na.rm=TRUE))
            checkIdentical(median(list1, na.rm=TRUE),
                           sapply(list1, median, na.rm=TRUE))
            checkIdentical(quantile(list1, na.rm=TRUE),
                           sapply(list1, quantile, na.rm=TRUE))
            checkIdentical(mad(list1, na.rm=TRUE),
                           sapply(list1, mad, na.rm=TRUE))
            checkIdentical(IQR(list1, na.rm=TRUE),
                           sapply(list1, IQR, na.rm=TRUE))

            vec3 <- (-20:20)^2
            vec3[c(1,10,21,41)] <- c(100L, 30L, 400L, 470L)
            list3 <- do.call(type, list(one = vec3, rev(vec3), compress = compress))
            checkIdentical(as.list(smoothEnds(list3)), lapply(list3, smoothEnds))
            checkIdentical(as.list(runmed(list3, 7)),
                           lapply(list3, function(x) {
                                      y <- runmed(x, 7)
                                      if (type != "RleList")
                                          y <- as.vector(y)
                                      y
                                  }))
        }
    }
}

test_AtomicList_character <- function() {
    txt <- c("The", "licenses", "for", "most", "software", "are",
             "designed", "to", "take", "away", "your", "freedom",
             "to", "share", "and", "change", "it.",
             "", "By", "contrast,", "the", "GNU", "General", "Public", "License",
             "is", "intended", "to", "guarantee", "your", "freedom", "to",
             "share", "and", "change", "free", "software", "--",
             "to", "make", "sure", "the", "software", "is",
             "free", "for", "all", "its", "users")
     for (compress in c(TRUE, FALSE)) {
         for (type in c("CharacterList", "RleList")) {
             list1 <- do.call(type, list(one = txt, rev(txt), compress = compress))
             checkIdentical(as.list(nchar(list1)), lapply(list1, nchar))
             checkIdentical(as.list(chartr("@!*", "alo", list1)),
                            lapply(list1, chartr, old="@!*", new="alo"))
             checkIdentical(as.list(tolower(list1)), lapply(list1, tolower))
             checkIdentical(as.list(toupper(list1)), lapply(list1, toupper))
             checkIdentical(as.list(sub("[b-e]",".", list1)),
                            lapply(list1, sub, pattern="[b-e]", replacement="."))
             checkIdentical(as.list(gsub("[b-e]",".", list1)),
                            lapply(list1, gsub, pattern="[b-e]", replacement="."))
        }
    }
}

test_RleList_methods <- function() {
    ## na.rm
    x <- RleList(c(NA,1,1), 
                 c(1L,NA_integer_,1L), 
                 c(1,Inf,1,-Inf),compress=TRUE)

    target <- RleList(c(1,2), c(1L,1L), c(Inf,Inf,-Inf))
    current <- runsum(x,2, na.rm = TRUE)
    checkIdentical(target, current)
    target <- RleList(c(NA,2), c(NA_integer_,NA_integer_), c(Inf,Inf,-Inf))
    current <- runsum(x,2, na.rm = FALSE)
    checkIdentical(target, current)

    target <- RleList(c(2,4), c(2,2), c(Inf, Inf, -Inf))
    current <- runwtsum(x,2, c(2,2), na.rm = TRUE)
    checkIdentical(target, current)
    target <- RleList(c(NA,4), c(NA_real_,NA_real_), c(Inf,Inf,-Inf))
    current <- runwtsum(x,2, c(2,2), na.rm = FALSE)
    checkIdentical(target, current)

    target <- RleList(c(1,1), c(1,1), c(Inf,Inf,-Inf))
    current <- runmean(x, 2, na.rm = TRUE)
    checkIdentical(target, current)
    target <- RleList(c(NA,1), c(NA_real_, NA_real_), c(Inf, Inf, -Inf))
    current <- runmean(x, 2, na.rm = FALSE)
    checkIdentical(target, current)

    x <- RleList(c(NA,1,2), 
                 c(2L,NA_integer_,1L), 
                 c(1,Inf,1,-Inf),compress=TRUE)
    target <- RleList(c(1,2), c(2L,1L), c(Inf,Inf,1))
    current <- runq(x, 2, 2, na.rm = TRUE)
    checkIdentical(target, current)
    target <- RleList(c(NA,2), c(NA_integer_, NA_integer_), c(Inf, Inf, 1))
    current <- runq(x, 2, 2, na.rm = FALSE)
    checkIdentical(target, current)

    ## Binary operations between an RleList and an atomic vector:
    a1 <- Rle(1, 999722111)
    a2 <- 20 * a1
    a <- RleList(a1, a2, compress=TRUE)
    b1 <- c(a1, a1)
    b2 <- 20 * b1
    b <- RleList(b1, b2, compress=FALSE)
    for (y in list(8L, 8)) {
        ## With a CompressedRleList
        target <- RleList(a1 + y, a2 + y, compress=TRUE)
        current <- a + y
        checkIdentical(target, current)
        target <- RleList(a1 * y, a2 * y, compress=TRUE)
        current <- a * y
        checkIdentical(target, current)
        target <- RleList(a1 / y, a2 / y, compress=TRUE)
        current <- a / y
        checkIdentical(target, current)
        ## With a SimpleRleList
        target <- RleList(b1 + y, b2 + y, compress=FALSE)
        current <- b + y
        checkIdentical(target, current)
        target <- RleList(b1 * y, b2 * y, compress=FALSE)
        current <- b * y
        checkIdentical(target, current)
        target <- RleList(b1 / y, b2 / y, compress=FALSE)
        current <- b / y
        checkIdentical(target, current)
    }
}

