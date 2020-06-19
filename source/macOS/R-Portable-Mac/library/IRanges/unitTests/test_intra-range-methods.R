test_shift_IntegerRanges <- function() {
  ir0 <- IRanges(0, 0)

  ir1 <- shift(ir0, .Machine$integer.max)
  checkTrue(validObject(ir1))
  checkIdentical(.Machine$integer.max, start(ir1))
  checkIdentical(.Machine$integer.max, end(ir1))
  checkIdentical(1L, width(ir1))

  checkIdentical(ir1, shift(ir1))
  checkIdentical(ir1, shift(shift(ir1, -10), 10))

  ir2 <- shift(ir0, -.Machine$integer.max)
  checkTrue(validObject(ir2))
  checkIdentical(-.Machine$integer.max, start(ir2))
  checkIdentical(-.Machine$integer.max, end(ir2))
  checkIdentical(1L, width(ir2))

  checkIdentical(ir2, shift(ir2))
  checkIdentical(ir2, shift(shift(ir2, 10), -10))

  ## shift() would produce an object with ranges that are not within the
  ## [-.Machine$integer.max, .Machine$integer.max] range.
  checkException(suppressWarnings(shift(ir1, 1)), silent=TRUE)
  checkException(suppressWarnings(shift(ir2, -1)), silent=TRUE)

  ir3 <- IRanges(1999222000, width=1000)
  checkException(suppressWarnings(shift(ir3, 188222000)), silent=TRUE)

  ir4 <- IRanges(1:20, width=222000000)
  checkException(suppressWarnings(shift(ir4, 1:20 * 99000000L)), silent=TRUE)
}

test_narrow_IntegerRanges <- function() {
  ir1 <- IRanges(c(2,5,1), c(3,7,3))
  checkIdentical(narrow(ir1, start=1, end=2),
                 IRanges(c(2, 5, 1), c(3, 6, 2)))
  checkException(narrow(ir1, start=10, end=20), silent = TRUE)
}

test_narrow_IRangesList <- function() {
  range1 <- IRanges(start=c(2,5), end=c(3,7))
  range2 <- IRanges(start=1, end=3)
  for (compress in c(TRUE, FALSE)) {
    collection <- IRangesList(range1, range2, compress = compress)
    checkIdentical(narrow(collection, start=1, end=2),
                   IRangesList(IRanges(c(2, 5), c(3, 6)), IRanges(1, 2),
                               compress = compress))
    checkException(narrow(collection, start=10, end=20), silent = TRUE)
  }
}

test_resize_IntegerRanges <- function() {
  ir1 <- IRanges(c(2,5,1), c(3,7,3))
  checkIdentical(resize(ir1, width=10),
                 IRanges(c(2, 5, 1), width=10))
  checkIdentical(resize(ir1, width=10, fix="end"),
                 IRanges(c(-6, -2, -6), width=10))
  checkIdentical(resize(ir1, width=10, fix="center"),
                 IRanges(c(-2, 1, -3), width=10))
  checkIdentical(resize(ir1, width=10, fix=c("start", "end", "center")),
                 IRanges(c(2, -2, -3), width=10))
  checkException(resize(ir1, -1), silent = TRUE)
}

test_resize_IRangesList <- function() {
  range1 <- IRanges(start=c(2,5), end=c(3,7))
  range2 <- IRanges(start=1, end=3)
  for (compress in c(TRUE, FALSE)) {
    collection <- IRangesList(range1, range2, compress = compress)
    checkIdentical(resize(collection, width=10),
                   IRangesList(IRanges(c(2, 5), width=10), IRanges(1, width=10),
                               compress = compress))
    checkIdentical(resize(collection, width=10, fix="end"),
                   IRangesList(IRanges(c(-6, -2), width=10), IRanges(-6, width=10),
                               compress = compress))
    checkIdentical(resize(collection, width=10, fix="center"),
                   IRangesList(IRanges(c(-2, 1), width=10), IRanges(-3, width=10),
                               compress = compress))
    checkIdentical(resize(collection, width=10,
                          fix=CharacterList(c("start", "end"), "center")),
                   IRangesList(IRanges(c(2, -2), width=10), IRanges(-3, width=10),
                               compress = compress))
    checkException(resize(collection, -1), silent = TRUE)
  }
}

test_flank_IntegerRanges <- function() {
  checkIdentical(flank(IRanges(), 2), IRanges())

  ir1 <- IRanges(c(2, 5, 1), c(3, 7, 3))
  checkIdentical(flank(ir1, 2), IRanges(c(0, 3, -1), c(1, 4, 0)))
  checkIdentical(flank(ir1, 2, FALSE), IRanges(c(4, 8, 4), c(5, 9, 5)))
  checkIdentical(flank(ir1, 2, c(FALSE, TRUE, FALSE)),
                 IRanges(c(4, 3, 4), c(5, 4, 5)))
  checkIdentical(flank(ir1, c(2, -2, 2)), IRanges(c(0, 5, -1), c(1, 6, 0)))
  checkIdentical(flank(ir1, 2, both = TRUE), IRanges(c(0, 3, -1), c(3, 6, 2)))
  checkIdentical(flank(ir1, 2, FALSE, TRUE), IRanges(c(2, 6, 2), c(5, 9, 5)))
  checkIdentical(flank(ir1, -2, FALSE, TRUE), IRanges(c(2, 6, 2), c(5, 9, 5)))
  checkException(flank(ir1, 2, both = c(TRUE, FALSE, TRUE)),
                 silent = TRUE) # not vectorized
  checkException(flank(ir1, 2, c(FALSE, TRUE, NA)), silent = TRUE)
  checkException(flank(ir1, NA), silent = TRUE)
}

test_flank_IRangesList <- function() {
  range1 <- IRanges(start=c(2,5), end=c(3,7))
  range2 <- IRanges(start=1, end=3)
  for (compress in c(TRUE, FALSE)) {
    collection <- IRangesList(range1, range2, compress = compress)
    checkIdentical(flank(collection, 2),
                   IRangesList(IRanges(c(0, 3), c(1, 4)), IRanges(-1, 0),
                               compress = compress))
    checkIdentical(flank(collection, 2, FALSE),
                   IRangesList(IRanges(c(4, 8), c(5, 9)), IRanges(4, 5),
                               compress = compress))
    checkIdentical(flank(collection, 2, LogicalList(c(FALSE, TRUE), FALSE)),
                   IRangesList(IRanges(c(4, 3), c(5, 4)), IRanges(4, 5),
                               compress = compress))
    checkIdentical(flank(collection, IntegerList(c(2, -2), 2)),
                   IRangesList(IRanges(c(0, 5), c(1, 6)), IRanges(-1, 0),
                               compress = compress))
    checkIdentical(flank(collection, 2, both = TRUE),
                   IRangesList(IRanges(c(0, 3), c(3, 6)), IRanges(-1, 2),
                               compress = compress))
    checkIdentical(flank(collection, 2, FALSE, TRUE),
                   IRangesList(IRanges(c(2, 6), c(5, 9)), IRanges(2, 5),
                               compress = compress))
    checkIdentical(flank(collection, -2, FALSE, TRUE),
                   IRangesList(IRanges(c(2, 6), c(5, 9)),
                               IRanges(2, 5), compress = compress))
    checkException(flank(collection, 2, both = c(TRUE, FALSE, TRUE)),
                   silent = TRUE) # not vectorized
    checkException(flank(collection, 2, LogicalList(c(FALSE, TRUE), NA)),
                   silent = TRUE)
    checkException(flank(collection, NA), silent = TRUE)
  }
}

test_promoters <- function() {
  ir <- IRanges(c(10, 10), width=c(0, 1))
  checkIdentical(width(promoters(ir, 0, 0)), c(0L, 0L))
  checkIdentical(width(promoters(ir, 1, 0)), c(1L, 1L))
  checkIdentical(start(promoters(ir, 1, 0)), c(9L, 9L))
  checkIdentical(width(promoters(ir, 0, 1)), c(1L, 1L))
  checkIdentical(start(promoters(ir, 0, 1)), c(10L, 10L))
  ir <- IRanges(c(5, 2, 20), width=1)
  checkIdentical(start(promoters(ir, 5, 2)), c(0L, -3L, 15L))

  rl <- IRangesList("A"=IRanges(5:7, width=1), "B"=IRanges(10:12, width=5))
  current <- promoters(rl, 0, 0)
  checkIdentical(names(current), names(rl))
  checkIdentical(start(current), start(rl))
  current <- promoters(rl, 2, 0)
  checkIdentical(unique(unlist(width(current))), 2L)

  library(XVector)
  subject <- XInteger(10, 3:-6)
  view <- Views(subject, start=4:2, end=4:6)
  current <- promoters(view, 0, 0)
  checkIdentical(start(current), start(view))
  current <- promoters(view, 2, 0)
  checkIdentical(unique(width(current)), 2L)

  cmp <- IRangesList("A"=IRanges(5:7, width=1), "B"=IRanges(10:12, width=5))
  current <- promoters(rl, 0, 0)
  checkIdentical(names(current), names(rl))
  checkIdentical(start(current), start(rl))
  current <- promoters(rl, 2, 0)
  checkIdentical(unique(unlist(width(current))), 2L)
}

test_reflect_IntegerRanges <- function() {
  ir1 <- IRanges(c(2,5,1), c(3,7,3))
  bounds <- IRanges(c(0, 5, 3), c(10, 6, 9))
  checkIdentical(reflect(ir1, bounds),
                 IRanges(c(7, 4, 9), c(8, 6, 11)))
  checkException(reflect(ir1, IRanges()), silent = TRUE)
}

test_restrict_IntegerRanges <- function() {
  ir1 <- IRanges(c(2,5,1), c(3,7,3))
  checkIdentical(restrict(ir1, start=2, end=5),
                 IRanges(c(2, 5, 2), c(3, 5, 3)))
  checkIdentical(restrict(ir1, start=1, end=2),
                 IRanges(c(2, 1), c(2, 2)))
  checkIdentical(restrict(ir1, start=1, end=2, keep.all.ranges=TRUE),
                 IRanges(c(2, 3, 1), c(2, 2, 2)))
}

test_restrict_IRangesList <- function() {
  range1 <- IRanges(start=c(2,5), end=c(3,7))
  range2 <- IRanges(start=1, end=3)
  for (compress in c(TRUE, FALSE)) {
    collection <- IRangesList(range1, range2, compress = compress)
    checkIdentical(restrict(collection, start=2, end=5),
                   IRangesList(IRanges(c(2, 5), c(3, 5)), IRanges(2, 3),
                               compress = compress))
    checkIdentical(restrict(collection, start=1, end=2),
                   IRangesList(IRanges(2, 2), IRanges(1, 2),
                               compress = compress))
    checkIdentical(restrict(collection, start=1, end=2, keep.all.ranges=TRUE),
                   IRangesList(IRanges(c(2, 3), c(2, 2)), IRanges(1, 2),
                               compress = compress))
  }
}

test_zoom_IntegerRanges <- function() {
  ir <- IRanges(c(1,5), c(3,10))
  checkIdentical(ir*1, ir)
  checkIdentical(ir*c(1,2), IRanges(c(1,6), c(3, 8)))
  checkIdentical(ir*-2, IRanges(c(-1,2), c(4, 13)))
  checkException(ir*NA_integer_, silent = TRUE)
  checkException(ir*numeric(), silent = TRUE)
  checkException(ir*c(1,2,1), silent = TRUE)
  checkException(ir[rep(1,3)]*c(1,2), silent = TRUE)
}

