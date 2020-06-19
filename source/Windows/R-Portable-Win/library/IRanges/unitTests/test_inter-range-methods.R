test_range_IntegerRanges <- function() {
  ir1 <- IRanges(c(2,5,1), c(3,7,3))
  ir2 <- IRanges(c(5,2,0), c(6,3,1))
  checkIdentical(range(ir1), IRanges(1, 7))
  checkIdentical(range(ir1, ir2), IRanges(0, 7))
  checkIdentical(range(IRanges()), IRanges())
  checkException(range(ir1, c(2,3)), silent = TRUE)
  # check with.revmap
  rng1 <- range(ir1, with.revmap=TRUE)
  rng2 <- range(ir2, with.revmap=TRUE)
  rng3 <- range(ir1,ir2, with.revmap=TRUE)
  checkIdentical(mcols(rng1)$revmap, IntegerList(seq(3)))
  checkIdentical(mcols(rng2)$revmap, IntegerList(seq(3)))
  checkIdentical(mcols(rng3)$revmap, IntegerList(seq(6)))
  ir3 <- IRanges()
  checkIdentical(mcols(range(ir3, with.revmap=TRUE))$revmap, IntegerList())
}

test_range_IntegerRangesList <- function() {
  for (compress in c(TRUE, FALSE)) {
    rl1 <- IRangesList(a = IRanges(c(1,2),c(4,3)), b = IRanges(c(4,6),c(10,7)),
                       compress = compress)
    rl2 <- IRangesList(c = IRanges(c(0,2),c(4,5)), a = IRanges(c(4,5),c(6,7)),
                       compress = compress)
    ans <- IRangesList(a = IRanges(1,7), b = IRanges(4,10), c = IRanges(0,5),
                       compress = compress)
    checkIdentical(range(rl1, rl2), ans)
    names(rl2) <- NULL
    ans <- IRangesList(IRanges(0,5), IRanges(4,10), compress = compress)
    checkIdentical(range(rl1, rl2), ans)
    ## must be same length
    checkException(range(rl2, rep.int(rl2, 2L)), silent=TRUE)
  }
  # check with.revmap
  revmap1 <- mcols(range(rl1,rl2, with.revmap=TRUE)[[1]])$revmap
  revmap2 <- mcols(range(rl1,rl2, with.revmap=TRUE)[[2]])$revmap
  ans <- IntegerList(seq(4))
  checkIdentical(revmap1, ans)
  checkIdentical(revmap2, ans)
  range1 <- IRanges(start=c(1, 2, 3), end=c(5, 2, 8))
  range2 <- IRanges(start=c(15, 45, 20, 1), end=c(15, 100, 80, 5))
  range3 <- IRanges()
  range4 <- IRanges(start=c(-2, 6, 7), width=c(8, 0, 0))
  collection <- IRangesList(range1, range2, range3, range4)
  rng <- range(collection, with.revmap=TRUE)
  checkIdentical(mcols(rng[[1]])$revmap, IntegerList(1:3))
  checkIdentical(mcols(rng[[2]])$revmap, IntegerList(1:4))
  checkIdentical(mcols(rng[[3]])$revmap, IntegerList())
  checkIdentical(mcols(rng[[4]])$revmap, IntegerList(1:3))
  rng <- range(IRangesList(IRanges(), IRanges()), with.revmap=TRUE)
  checkIdentical(mcols(rng[[1]])$revmap, IntegerList())
  checkIdentical(mcols(rng[[2]])$revmap, IntegerList())
}

test_reduce_IntegerRanges <- function() {
  x <- IRanges()
  current <- reduce(x)
  checkIdentical(x, current)

  x <- IRanges(1:3, width=0)

  current <- reduce(x, with.revmap=TRUE)
  target <- x
  mcols(target) <- DataFrame(revmap=as(seq_along(target), "IntegerList"))
  checkIdentical(target, current)

  current <- reduce(x, drop.empty.ranges=TRUE, with.revmap=TRUE)
  target <- IRanges()
  mcols(target) <- DataFrame(revmap=IntegerList(seq_along(target)))
  checkIdentical(target, current)

  x <- IRanges(c(1:4, 10:11, 11), width=c(0,1,1,0,0,0,1))

  current <- reduce(x, with.revmap=TRUE)
  target <- IRanges(c(1:2, 10:11), width=c(0,2,0,1))
  mcols(target) <- DataFrame(revmap=IntegerList(1,2:4,5,6:7))
  checkIdentical(target, current)

  current <- reduce(x, drop.empty.ranges=TRUE, with.revmap=TRUE)
  target <- IRanges(c(2, 11), width=c(2,1))
  mcols(target) <- DataFrame(revmap=IntegerList(2:3,7))
  checkIdentical(target, current)

  x <- IRanges(start=c(1,2,3), end=c(5,2,8))
  y <- reduce(x, with.revmap=TRUE)

  target <- IRanges(start=1, end=8)
  mcols(target) <- DataFrame(revmap=IntegerList(1:3))
  checkIdentical(target, y)

  mcols(target)$revmap <- as(seq_along(target), "IntegerList")
  checkIdentical(target, reduce(y, with.revmap=TRUE))

  x <- IRanges(start=c(15,45,20,1), end=c(15,100,80,5))
  y <- reduce(x, with.revmap=TRUE)

  target <- IRanges(start=c(1,15,20), end=c(5,15,100))
  mcols(target) <- DataFrame(revmap=IntegerList(4, 1, 3:2))
  checkIdentical(target, y)

  mcols(target)$revmap <- as(seq_along(target), "IntegerList")
  checkIdentical(target, reduce(y, with.revmap=TRUE))

  x <- IRanges(start=c(7,3,-2,6,7,-10,-2,3), width=c(3,1,0,0,0,0,8,0))
  ## Before reduction:
  ##     start end width  ==-10===-5====0===+5==+10===
  ## [1]     7   9     3  ....:....:....:....:.xxx:...
  ## [2]     3   3     1  ....:....:....:..x.:....:...
  ## [3]    -2  -3     0  ....:....:..[.:....:....:...
  ## [4]     6   5     0  ....:....:....:....:[...:...
  ## [5]     7   6     0  ....:....:....:....:.[..:...
  ## [6]   -10 -11     0  ....[....:....:....:....:...
  ## [7]    -2   5     8  ....:....:..xxxxxxxx....:...
  ## [8]     3   2     0  ....:....:....:..[.:....:...
  ## ---------------------==-10===-5====0===+5==+10===
  ## After reduction:
  ##                  y1: ....[....:..xxxxxxxx.xxx:...
  ##                  y3: ....:....:..xxxxxxxx....:...
  y1 <- reduce(x)
  checkIdentical(y1, IRanges(start=c(-10,-2,7), end=c(-11,5,9)))
  checkIdentical(reduce(y1), y1)
  y2 <- reduce(x, with.inframe.attrib=TRUE)
  checkIdentical(start(attr(y2, "inframe")), c(9L,6L,1L,9L,9L,1L,1L,6L))
  checkIdentical(width(attr(y2, "inframe")), width(x))
  y3 <- reduce(x, drop.empty.ranges=TRUE)
  checkIdentical(y3, y1[width(y1) != 0L])
  checkIdentical(reduce(y3), y3)
  y4 <- reduce(x, drop.empty.ranges=TRUE, with.inframe.attrib=TRUE)
  checkIdentical(attr(y4, "inframe"), attr(y2, "inframe"))
  y5 <- reduce(x, min.gapwidth=0)
  checkIdentical(y5, IRanges(start=c(-10,-2,-2,6,7,7), end=c(-11,-3,5,5,6,9)))
  y6 <- reduce(x, drop.empty.ranges=TRUE, min.gapwidth=0)
  checkIdentical(y6, y5[width(y5) != 0L])
  y7 <- reduce(x, min.gapwidth=2)
  checkIdentical(y7, IRanges(start=c(-10,-2), end=c(-11,9)))
  y8 <- reduce(x, min.gapwidth=8)
  checkIdentical(y8, y7)
  y9 <- reduce(x, min.gapwidth=9)
  checkIdentical(y9, IRanges(start=-10, end=9))
}

test_reduce_IntegerRangesList <- function() {
  range1 <- IRanges(start=c(1,2,3), end=c(5,2,8))
  range2 <- IRanges(start=c(15,45,20,1), end=c(15,100,80,5))
  range3 <- IRanges(start=c(3,-2,6,7,-10,-2,3), width=c(1,0,0,0,0,8,0))
  range4 <- IRanges()
  for (compress in c(TRUE, FALSE)) {
    collection <- IRangesList(one=range1,
                              range2,
                              range3,
                              range4,
                              compress=compress)
    for (with.revmap in c(FALSE, TRUE)) {
      for (drop.empty.ranges in c(FALSE, TRUE)) {
        current <- reduce(collection, drop.empty.ranges=drop.empty.ranges,
                                      with.revmap=with.revmap)
        target <- IRangesList(one=reduce(range1,
                                         drop.empty.ranges=drop.empty.ranges,
                                         with.revmap=with.revmap),
                              reduce(range2,
                                     drop.empty.ranges=drop.empty.ranges,
                                     with.revmap=with.revmap),
                              reduce(range3,
                                     drop.empty.ranges=drop.empty.ranges,
                                     with.revmap=with.revmap),
                              reduce(range4,
                                     drop.empty.ranges=drop.empty.ranges,
                                     with.revmap=with.revmap),
                              compress=compress)
        checkIdentical(target, current)
      }
    }
  }
}

test_gaps_IntegerRanges <- function() {
  checkIdentical(gaps(IRanges()), IRanges())
  checkIdentical(gaps(IRanges(), start=1, end=4),
                 IRanges(start=1, end=4))

  x <- IRanges(start=2, end=3)
  checkIdentical(gaps(x), IRanges())
  checkIdentical(gaps(x, start=2), IRanges())
  checkIdentical(gaps(x, start=4), IRanges())
  checkIdentical(gaps(x, start=0), IRanges(start=0, end=1))
  checkIdentical(gaps(x, end=3), IRanges())
  checkIdentical(gaps(x, end=1), IRanges())
  checkIdentical(gaps(x, end=5), IRanges(start=4, end=5))
  checkIdentical(gaps(x, start=0, end=5), IRanges(start=c(0,4), end=c(1,5)))
}

test_gaps_IntegerRangesList <- function() {
  range1 <- IRanges(start=c(1,2,3), end=c(5,2,8))
  range2 <- IRanges(start=c(15,45,20,1), end=c(15,100,80,5))
  for (compress in c(TRUE, FALSE)) {
    collection <- IRangesList(one = range1, range2, compress = compress)
    checkIdentical(gaps(collection),
                   IRangesList(one = gaps(range1), gaps(range2),
                               compress = compress))
  }
}

test_disjoin_IntegerRanges <- function()
{
  checkIdentical(disjoin(IRanges()), IRanges())
  ir <- IRanges(c(1, 21, 10, 1, 15, 5, 20, 20), c(6, 20, 9, 3, 14, 11, 20, 19))
  current <- disjoin(ir)
  checkTrue(validObject(current, complete=TRUE))

  ## The result of disjoin(x) must verify the following properties:
  check_disjoin_general_properties <- function(y, x) {
      checkTrue(isDisjoint(y))
      checkTrue(isStrictlySorted(y))
      checkIdentical(reduce(x, drop.empty.ranges=TRUE), reduce(y))
      checkTrue(all(start(y) %in% c(start(x), end(x) + 1L)))
      checkTrue(all(end(y) %in% c(end(x), start(x) - 1L)))
  }

  check_disjoin_general_properties(current, ir)

  target <- IRanges(c(1, 4, 5, 7, 10, 20), c(3, 4, 6, 9, 11, 20))
  checkIdentical(target, current)

  ## Check 'revmap'.
  mcols(ir)$label <- LETTERS[seq_along(ir)]
  current <- disjoin(ir, with.revmap=TRUE)
  revmap <- IntegerList(c(1, 4), 1, c(1, 6), 6, 6, 7)
  mcols(target)$revmap <- revmap
  checkIdentical(target, current)

  ## With many randomly generated ranges.
  set.seed(2009L)
  ir <- IRanges(start=sample(580L, 500L, replace=TRUE),
                width=sample(10L, 500L, replace=TRUE) - 1L)
  check_disjoin_general_properties(disjoin(ir), ir)
  ir <- IRanges(start=sample(4900L, 500L, replace=TRUE),
                width=sample(35L, 500L, replace=TRUE) - 1L)
  check_disjoin_general_properties(disjoin(ir), ir)
}

test_disjoin_IntegerRangesList <- function()
{
    ir0 <- IRanges(10, 20)
    checkTrue(validObject(disjoin(IRangesList())))
    ## unnamed; incl. 0-length
    irl <- IRangesList(IRanges())
    checkIdentical(irl, disjoin(irl))
    irl <- IRangesList(ir0, IRanges(), ir0)
    checkIdentical(irl, disjoin(irl))
    irl <- IRangesList(ir0, IRanges(), IRanges(), ir0)
    checkIdentical(irl, disjoin(irl))
    ## named; incl. 0-length
    irl <- IRangesList(a=IRanges())
    checkIdentical(irl, disjoin(irl))
    irl <- IRangesList(a=ir0, b=IRanges(), c=ir0)
    checkIdentical(irl, disjoin(irl))
    irl <- IRangesList(a=ir0, b=IRanges(), c=IRanges(), d=ir0)
    checkIdentical(irl, disjoin(irl))
    ## no interference between separate elements
    ir0 <- IRanges(10, c(15, 20))
    dr0 <- disjoin(ir0)
    irl <- IRangesList(ir0, ir0)
    checkIdentical(IRangesList(dr0, dr0), disjoin(irl))
    irl <- IRangesList(ir0, IRanges(), ir0)
    checkIdentical(IRangesList(dr0, IRanges(), dr0), disjoin(irl))
    ## 0-width
    ## 1-width
    ir0 <- IRanges(c(1, 10), 10)
    irl <- IRangesList(ir0, IRanges())
    checkIdentical(disjoin(ir0), disjoin(irl)[[1]])
    irl <- IRangesList(IRanges(), ir0)
    checkIdentical(disjoin(ir0), disjoin(irl)[[2]])

    ## check don't collapse levels
    irl <- IRangesList(IRanges(1, 5), IRanges(3, 7))
    names(irl) <- character(2)
    checkIdentical(irl, disjoin(irl))

    ## check 'revmap' on many randomly generated ranges
    set.seed(2009L)
    ir1 <- IRanges(start=sample(580L, 500L, replace=TRUE),
                   width=sample(10L, 500L, replace=TRUE) - 1L)
    ir2 <- IRanges(start=sample(4900L, 500L, replace=TRUE),
                   width=sample(35L, 500L, replace=TRUE) - 1L)
    for (compress in c(TRUE, FALSE)) {
      collection <- IRangesList(one=ir1,
                                IRanges(),
                                ir0,
                                ir0,
                                ir2,
                                IRanges(),
                                compress=compress)
      for (with.revmap in c(FALSE, TRUE)) {
        current <- disjoin(collection, with.revmap=with.revmap)
        target <- IRangesList(one=disjoin(ir1, with.revmap=with.revmap),
                              disjoin(IRanges(), with.revmap=with.revmap),
                              disjoin(ir0, with.revmap=with.revmap),
                              disjoin(ir0, with.revmap=with.revmap),
                              disjoin(ir2, with.revmap=with.revmap),
                              disjoin(IRanges(), with.revmap=with.revmap),
                              compress=compress)
        checkIdentical(target, current)
      }
    }
}

test_disjointBins_IntegerRanges <- function()
{
  checkIdentical(disjointBins(IRanges()), integer())
  checkIdentical(disjointBins(IRanges(1, 5)), 1L)
  checkIdentical(disjointBins(IRanges(c(1, 3), c(5, 12))), c(1L, 2L))
  checkIdentical(disjointBins(IRanges(c(1, 3, 10), c(5, 12, 13))),
                 c(1L, 2L, 1L))
  checkIdentical(disjointBins(IRanges(c(3, 1, 10), c(5, 12, 13))),
                 c(2L, 1L, 2L))
}

