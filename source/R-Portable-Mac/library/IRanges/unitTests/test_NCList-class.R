###

findOverlaps_NCList <- IRanges:::findOverlaps_NCList
findOverlaps_NCLists <- IRanges:::findOverlaps_NCLists

.transpose_hits <- function(hits)
{
    if (is.list(hits))
        return(lapply(hits, .transpose_hits))
    t(hits)
}

### Used in the unit tests for GNCList located in GenomicRanges.
.compare_hits <- function(target, current)
{
    if (is.list(target) || is(target, "List")
     && is.list(current) || is(current, "List"))
        return(all(mapply(.compare_hits, target, current)))
    identical(.transpose_hits(target), .transpose_hits(current))
}

### Used in the unit tests for GNCList located in GenomicRanges.
.make_Hits_from_q2s <- function(q2s, s_len)
{
    q_hits <- rep.int(seq_along(q2s), elementNROWS(q2s))
    s_hits <- as.integer(unlist(q2s, use.names=FALSE))
    Hits(q_hits, s_hits, length(q2s), s_len, sort.by.query=TRUE)
}

.make_Hits_from_s2q <- function(s2q, q_len)
    .transpose_hits(.make_Hits_from_q2s(s2q, q_len))

.select_hits <- function(x, select)
{
    if (is.list(x))
        return(lapply(x, .select_hits, select))
    selectHits(x, select)
}

### Vectorized. Return -1 if the query and subject overlap (i.e. if
### end(query) < start(subject) and end(subject) < start(query) are both
### false). Otherwise (i.e. if they are disjoint), return the width of the
### gap between them. Note that a gap width of 0 means that they are adjacent.
### TODO: Rename this pgapWidth(), make it a generic with various methods
### (at least one for IntegerRanges and one for GenomicRanges objects), and
### export it.
.gapwidth <- function(query, subject)
{
    ifelse(end(query) < start(subject),
           start(subject) - end(query),
           ifelse(end(subject) < start(query),
                  start(query) - end(subject),
                  0L)) - 1L
}

### Vectorized.
### TODO: Rename this poverlapWidth(), make it a generic with various methods
### (at least one for IntegerRanges and one for GenomicRanges objects), and
### export it.
.overlapwidth <- function(query, subject)
{
    score <- pmin.int(end(query), end(subject)) -
             pmax.int(start(query), start(subject)) + 1L
    pmax.int(score, 0L)
}

### Used in the unit tests for GNCList located in GenomicRanges.
.get_query_overlaps <- function(query, subject,
            maxgap=-1L, minoverlap=0L,
            type=c("any", "start", "end", "within", "extend", "equal"))
{
    type <- match.arg(type)
    if (type == "any" && maxgap != -1L && minoverlap != 0L)
        stop("when 'type' is \"any\", at least one of 'maxgap' ",
             "and 'minoverlap' must be set to its default value")
    overlapwidth <- .overlapwidth(query, subject)
    ok <- overlapwidth >= minoverlap
    if (type == "any") {
        gapwidth <- .gapwidth(query, subject)
        ok <- ok & gapwidth <= maxgap
        return(ok)
    }
    if (maxgap == -1L)
        maxgap <- 0L
    if (type != "end") 
        d1 <- abs(start(subject) - start(query))
    if (type != "start")
        d2 <- abs(end(subject) - end(query))
    if (type == "start")
        return(ok & d1 <= maxgap)
    if (type == "end")
        return(ok & d2 <= maxgap)
    if (type == "equal")
        return(ok & d1 <= maxgap & d2 <= maxgap)
    if (type == "within") {
        ok2 <- start(query) >= start(subject) & end(query) <= end(subject)
    } else {  # type == "extend"
        ok2 <- start(query) <= start(subject) & end(query) >= end(subject)
    }
    ok <- ok & ok2
    if (maxgap > 0L)
        ok <- ok & (d1 + d2) <= maxgap
    ok
}

.findOverlaps_naive <- function(query, subject,
                                maxgap=-1L, minoverlap=0L,
                                type=c("any", "start", "end",
                                       "within", "extend", "equal"),
                                select=c("all", "first", "last", "arbitrary",
                                         "count"))
{
    type <- match.arg(type)
    select <- match.arg(select)
    hits_per_query <- lapply(seq_along(query),
        function(i)
            which(.get_query_overlaps(query[i], subject,
                                      maxgap=maxgap, minoverlap=minoverlap,
                                      type=type)))
    hits <- .make_Hits_from_q2s(hits_per_query, length(subject))
    selectHits(hits, select=select)
}

test_NCList <- function()
{
    x <- IRanges(rep.int(1:6, 6:1), c(0:5, 1:5, 2:5, 3:5, 4:5, 5),
                 names=LETTERS[1:21])
    mcols(x) <- DataFrame(score=seq(0.7, by=0.045, length.out=21))
    nclist <- NCList(x)

    checkTrue(is(nclist, "NCList"))
    checkTrue(validObject(nclist, complete=TRUE))
    checkIdentical(x, ranges(nclist, use.mcols=TRUE))
    checkIdentical(length(x), length(nclist))
    checkIdentical(names(x), names(nclist))
    checkIdentical(start(x), start(nclist))
    checkIdentical(end(x), end(nclist))
    checkIdentical(width(x), width(nclist))
    checkIdentical(x, as(nclist, "IRanges"))
    checkIdentical(x[-6], as(nclist[-6], "IRanges"))
}

### Test findOverlaps_NCList() *default* behavior, that is, with all optional
### arguments (i.e. 'maxgap', 'minoverlap', 'type', 'select', and
### 'circle.length') set to their default value.
test_findOverlaps_NCList <- function()
{
    query <- IRanges(-3:7, width=3)
    subject <- IRanges(rep.int(1:6, 6:1), c(0:5, 1:5, 2:5, 3:5, 4:5, 5))

    target0 <- .findOverlaps_naive(query, subject)
    current <- findOverlaps_NCList(query, NCList(subject))
    checkTrue(.compare_hits(target0, current))
    current <- findOverlaps_NCList(NCList(query), subject)
    checkTrue(.compare_hits(target0, current))
    current <- findOverlaps_NCList(query, subject)
    checkTrue(.compare_hits(target0, current))

    ## Shuffle query and/or subject elements.
    permute_input <- function(q_perm, s_perm) {
        q_revperm <- integer(length(q_perm))
        q_revperm[q_perm] <- seq_along(q_perm)
        s_revperm <- integer(length(s_perm))
        s_revperm[s_perm] <- seq_along(s_perm)
        target <- remapHits(target0, Lnodes.remapping=q_revperm,
                                     new.nLnode=length(q_perm),
                                     Rnodes.remapping=s_revperm,
                                     new.nRnode=length(s_perm))
        current <- findOverlaps_NCList(query[q_perm], NCList(subject[s_perm]))
        checkTrue(.compare_hits(target, current))
        current <- findOverlaps_NCList(NCList(query[q_perm]), subject[s_perm])
        checkTrue(.compare_hits(target, current))
        current <- findOverlaps_NCList(query[q_perm], subject[s_perm])
        checkTrue(.compare_hits(target, current))
    }

    q_perm <- rev(seq_along(query))
    s_perm <- rev(seq_along(subject))
    permute_input(q_perm, seq_along(subject))  # reverse query
    permute_input(seq_along(query), s_perm)    # reverse subject
    permute_input(q_perm, s_perm)              # reverse both

    set.seed(97)
    for (i in 1:33) {
        ## random permutations
        q_perm <- sample(length(query))
        s_perm <- sample(length(subject))
        permute_input(q_perm, seq_along(subject))
        permute_input(seq_along(query), s_perm)
        permute_input(q_perm, s_perm)
    }
}

test_findOverlaps_NCList_with_filtering <- function()
{
    query <- IRanges(-3:7, width=3)
    subject <- IRanges(rep.int(1:6, 6:1), c(0:5, 1:5, 2:5, 3:5, 4:5, 5))

    pp_query <- NCList(query)
    pp_subject <- NCList(subject)
    for (type in c("any", "start", "end", "within", "extend", "equal")) {
      for (maxgap in -1:3) {
        if (type != "any" || maxgap == -1L)
            max_minoverlap <- 4L
        else
            max_minoverlap <- 0L
        for (minoverlap in 0:max_minoverlap) {
          for (select in c("all", "first", "last", "count")) {
            ## query - subject
            target <- .findOverlaps_naive(query, subject,
                                          maxgap=maxgap, minoverlap=minoverlap,
                                          type=type, select=select)
            current <- findOverlaps_NCList(query, pp_subject,
                                           maxgap=maxgap, minoverlap=minoverlap,
                                           type=type, select=select)
            checkTrue(.compare_hits(target, current))
            current <- findOverlaps_NCList(pp_query, subject,
                                           maxgap=maxgap, minoverlap=minoverlap,
                                           type=type, select=select)
            checkTrue(.compare_hits(target, current))
            current <- findOverlaps_NCList(query, subject,
                                           maxgap=maxgap, minoverlap=minoverlap,
                                           type=type, select=select)
            checkTrue(.compare_hits(target, current))
            ## subject - query
            target <- .findOverlaps_naive(subject, query,
                                          maxgap=maxgap, minoverlap=minoverlap,
                                          type=type, select=select)
            current <- findOverlaps_NCList(pp_subject, query,
                                           maxgap=maxgap, minoverlap=minoverlap,
                                           type=type, select=select)
            checkTrue(.compare_hits(target, current))
            current <- findOverlaps_NCList(subject, pp_query,
                                           maxgap=maxgap, minoverlap=minoverlap,
                                           type=type, select=select)
            checkTrue(.compare_hits(target, current))
            current <- findOverlaps_NCList(subject, query,
                                           maxgap=maxgap, minoverlap=minoverlap,
                                           type=type, select=select)
            checkTrue(.compare_hits(target, current))
            ## subject - subject
            target <- .findOverlaps_naive(subject, subject,
                                          maxgap=maxgap, minoverlap=minoverlap,
                                          type=type, select=select)
            current <- findOverlaps_NCList(pp_subject, subject,
                                           maxgap=maxgap, minoverlap=minoverlap,
                                           type=type, select=select)
            checkTrue(.compare_hits(target, current))
            current <- findOverlaps_NCList(subject, pp_subject,
                                           maxgap=maxgap, minoverlap=minoverlap,
                                           type=type, select=select)
            checkTrue(.compare_hits(target, current))
            current <- findOverlaps_NCList(subject, subject,
                                           maxgap=maxgap, minoverlap=minoverlap,
                                           type=type, select=select)
            checkTrue(.compare_hits(target, current))
          }
        }
      }
    }
}

### Only test "start" and "end" types at the moment.
test_findOverlaps_NCList_special_types <- function()
{
    x <- IRanges(10, 10)

    x1 <- IRanges(10, 9)
    y1 <- IRanges(start=c(7, 7, 13, 13), width=c(2, 0, 2, 0))
    stopifnot(all(abs(start(x) - start(y1)) == 3L))
    stopifnot(all(abs(start(x1) - start(y1)) == 3L))

    x2 <- IRanges(11, 10)
    y2 <- IRanges(end=c(7, 7, 13, 13), width=c(2, 0, 2, 0))
    stopifnot(all(abs(end(x) - end(y2)) == 3L))
    stopifnot(all(abs(end(x2) - end(y2)) == 3L))

    test_maxgap_and_type <- function(maxgap, minoverlap, nhit)
    {
        hits <- findOverlaps(x, y1, maxgap=maxgap, minoverlap=minoverlap,
                             type="start")
        checkEquals(nhit, length(hits))
        hits <- findOverlaps(y1, x, maxgap=maxgap, minoverlap=minoverlap,
                             type="start")
        checkEquals(nhit, length(hits))
        hits <- findOverlaps(x1, y1, maxgap=maxgap, minoverlap=minoverlap,
                             type="start")
        checkEquals(nhit, length(hits))
        hits <- findOverlaps(y1, x1, maxgap=maxgap, minoverlap=minoverlap,
                             type="start")
        checkEquals(nhit, length(hits))
        hits <- findOverlaps(x, y2, maxgap=maxgap, minoverlap=minoverlap,
                             type="end")
        checkEquals(nhit, length(hits))
        hits <- findOverlaps(y2, x, maxgap=maxgap, minoverlap=minoverlap,
                             type="end")
        checkEquals(nhit, length(hits))
        hits <- findOverlaps(x2, y2, maxgap=maxgap, minoverlap=minoverlap,
                             type="end")
        checkEquals(nhit, length(hits))
        hits <- findOverlaps(y2, x2, maxgap=maxgap, minoverlap=minoverlap,
                             type="end")
        checkEquals(nhit, length(hits))
    }
    ## no hits
    for (maxgap in -1:2) {
        test_maxgap_and_type(maxgap, minoverlap=1L, 0L)
        test_maxgap_and_type(maxgap, minoverlap=0L, 0L)
    }
    for (maxgap in 3:5) {
        ## no hits
        test_maxgap_and_type(maxgap, minoverlap=1L, 0L)
        ## 4 hits
        test_maxgap_and_type(maxgap, minoverlap=0L, 4L)
    }
}

.test_arbitrary_selection <- function(query, subject)
{
    pp_query <- NCList(query)
    pp_subject <- NCList(subject)
    for (type in c("any", "start", "end", "within", "extend", "equal")) {
      for (maxgap in -1:3) {
        if (type != "any" || maxgap == -1L)
            max_minoverlap <- 4L
        else
            max_minoverlap <- 0L
        for (minoverlap in 0:max_minoverlap) {
          target <- as(.findOverlaps_naive(query, subject,
                                           maxgap=maxgap, minoverlap=minoverlap,
                                           type=type, select="all"),
                       "CompressedIntegerList")
          target_idx0 <- elementNROWS(target) == 0L
          check_arbitrary_hits <- function(current) {
              current_idx0 <- is.na(current)
              checkIdentical(target_idx0, current_idx0)
              current <- as(current, "CompressedIntegerList")
              checkTrue(all(current_idx0 | as.logical(current %in% target)))
          }
          current <- findOverlaps_NCList(query, pp_subject,
                                         maxgap=maxgap, minoverlap=minoverlap,
                                         type=type, select="arbitrary")
          check_arbitrary_hits(current)
          current <- findOverlaps_NCList(pp_query, subject,
                                         maxgap=maxgap, minoverlap=minoverlap,
                                         type=type, select="arbitrary")
          check_arbitrary_hits(current)
          current <- findOverlaps_NCList(query, subject,
                                         maxgap=maxgap, minoverlap=minoverlap,
                                         type=type, select="arbitrary")
          check_arbitrary_hits(current)
        }
      }
    }
}

test_findOverlaps_NCList_arbitrary <- function()
{
    query <- IRanges(4:3, 6)
    subject <- IRanges(2:4, 10)
    .test_arbitrary_selection(query, subject)
    query <- IRanges(-3:7, width=3)
    subject <- IRanges(rep.int(1:6, 6:1), c(0:5, 1:5, 2:5, 3:5, 4:5, 5))
    .test_arbitrary_selection(query, subject)
}

.test_circularity <- function(query0, subject0, circle_length, target0,
                              pp, findOverlaps_pp, type)
{
    for (i in -2:2) {
        query <- shift(query0, shift=i*circle_length)
        pp_query <- pp(query, circle.length=circle_length)
        for (j in -2:2) {
            subject <- shift(subject0, shift=j*circle_length)
            pp_subject <- pp(subject, circle.length=circle_length)
            for (select in c("all", "first", "last", "count")) {
                target <- .select_hits(target0, select=select)
                current <- findOverlaps_pp(query, pp_subject,
                                           type=type, select=select,
                                           circle.length=circle_length)
                checkTrue(.compare_hits(target, current))
                current <- findOverlaps_pp(pp_query, subject,
                                           type=type, select=select,
                                           circle.length=circle_length)
                checkTrue(.compare_hits(target, current))
                current <- findOverlaps_pp(query, subject,
                                           type=type, select=select,
                                           circle.length=circle_length)
                checkTrue(.compare_hits(target, current))

                target <- .select_hits(.transpose_hits(target0), select=select)
                current <- findOverlaps_pp(pp_subject, query,
                                           type=type, select=select,
                                           circle.length=circle_length)
                checkTrue(.compare_hits(target, current))
                current <- findOverlaps_pp(subject, pp_query,
                                           type=type, select=select,
                                           circle.length=circle_length)
                checkTrue(.compare_hits(target, current))
                current <- findOverlaps_pp(subject, query,
                                           type=type, select=select,
                                           circle.length=circle_length)
                checkTrue(.compare_hits(target, current))
            }
        }
    }
}

test_findOverlaps_NCList_with_circular_space <- function()
{
    query <- IRanges(-2:17, width=3)
    subject <- IRanges(c(4, -1, 599), c(7, 0, 999))
    circle_length <- 10L

    ## type "any"
    s2q <- list(c(5:10, 15:20L), c(1:3, 10:13, 20L), 1:20)
    target <- .make_Hits_from_s2q(s2q, length(query))
    .test_circularity(query, subject, circle_length, target,
                      NCList, findOverlaps_NCList, "any")

    ## type "start"
    s2q <- lapply(start(subject),
                  function(s)
                      which((start(query) - s) %% circle_length == 0L))
    target <- .make_Hits_from_s2q(s2q, length(query))
    .test_circularity(query, subject, circle_length, target,
                      NCList, findOverlaps_NCList, "start")

    ## type "end"
    s2q <- lapply(end(subject),
                  function(e)
                      which((end(query) - e) %% circle_length == 0L))
    target <- .make_Hits_from_s2q(s2q, length(query))
    .test_circularity(query, subject, circle_length, target,
                      NCList, findOverlaps_NCList, "end")
}

test_NCLists <- function()
{
    x1 <- IRanges(-3:7, width=3)
    x2 <- IRanges()
    x3 <- IRanges(rep.int(1:6, 6:1), c(0:5, 1:5, 2:5, 3:5, 4:5, 5))
    x <- IRangesList(x1=x1, x2=x2, x3=x3)
    mcols(x) <- DataFrame(label=c("first", "second", "third"))
    nclists <- NCLists(x)

    checkTrue(is(nclists, "NCLists"))
    checkTrue(validObject(nclists, complete=TRUE))
    checkIdentical(x, ranges(nclists, use.mcols=TRUE))
    checkIdentical(length(x), length(nclists))
    checkIdentical(names(x), names(nclists))
    checkIdentical(start(x), start(nclists))
    checkIdentical(end(x), end(nclists))
    checkIdentical(width(x), width(nclists))
    checkIdentical(x, as(nclists, "IRangesList"))
    checkIdentical(x[-1], as(nclists[-1], "IRangesList"))
    checkIdentical(elementNROWS(x), elementNROWS(nclists))

    nclist <- nclists[[3]]
    checkTrue(is(nclist, "NCList"))
    checkTrue(validObject(nclist, complete=TRUE))
    checkIdentical(x3, as(nclist, "IRanges"))
}

test_findOverlaps_NCLists <- function()
{
    ir1 <- IRanges(-3:7, width=3)
    ir2 <- IRanges(rep.int(1:6, 6:1), c(0:5, 1:5, 2:5, 3:5, 4:5, 5))
    target0 <- mapply(findOverlaps_NCList, list(ir1, ir2), list(ir2, ir1))
    for (compress in c(TRUE, FALSE)) {
        query <- IRangesList(ir1, ir2, IRanges(2, 7), compress=compress)
        pp_query <- NCLists(query)
        subject <- IRangesList(ir2, ir1, compress=compress)
        pp_subject <- NCLists(subject)
        for (select in c("all", "first", "last", "count")) {
            target <- .select_hits(target0, select=select)
            current <- findOverlaps_NCLists(query, pp_subject, select=select)
            checkTrue(.compare_hits(target, current))
            current <- findOverlaps_NCLists(pp_query, subject, select=select)
            checkTrue(.compare_hits(target, current))
            current <- findOverlaps_NCLists(query, subject, select=select)
            checkTrue(.compare_hits(target, current))
        }
    }
}

test_findOverlaps_NCLists_with_circular_space <- function()
{
    query1 <- IRanges(-2:17, width=3)
    subject1 <- IRanges(c(4, -1, 599), c(7, 0, 999))
    query <- IRangesList(query1, IRanges(), subject1)
    subject <- IRangesList(subject1, IRanges(), query1)
    circle_length <- c(10L, NA_integer_, 10L)

    s2q <- list(c(5:10, 15:20L), c(1:3, 10:13, 20L), 1:20)
    target1 <- .make_Hits_from_s2q(s2q, length(query1))
    target2 <- .make_Hits_from_s2q(list(), 0)
    target3 <- .transpose_hits(target1)
    target <- list(target1, target2, target3)
    .test_circularity(query, subject, circle_length, target,
                      NCLists, findOverlaps_NCLists, "any")
}

