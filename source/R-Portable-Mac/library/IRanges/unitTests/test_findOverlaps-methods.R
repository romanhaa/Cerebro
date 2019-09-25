###

test_findOverlaps_IntegerRanges <- function()
{
  ## .....
  ##    ....
  ##         ..
  ##  x
  ##  xx
  ##          xxx
  query <- IRanges(c(1, 4, 9), c(5, 7, 10))
  subject <- IRanges(c(2, 2, 10), c(2, 3, 12))

  result <- findOverlaps(query, subject, select = "first")
  checkIdentical(result, c(1L, NA, 3L))
  result <- findOverlaps(query, subject, select = "last")
  checkIdentical(result, c(2L, NA, 3L))
  result <- findOverlaps(query, subject, select = "arbitrary")
  checkIdentical(result, c(2L, NA, 3L))

  checkOverlap <- function(a, q, s, r, c, self=FALSE) {
    target <- Hits(q, s, r, c, sort.by.query=TRUE)
    if (self)
        target <- as(target, "SortedByQuerySelfHits")
    checkIdentical(t(a), t(target))
  }

  result <- findOverlaps(query, subject)
  checkOverlap(result, c(1, 1, 3), c(1, 2, 3), 3, 3)

  ## with 'maxgap'
  result <- findOverlaps(query, subject, maxgap = 0L)
  checkOverlap(result, c(1, 1, 2, 3), c(2, 1, 2, 3), 3, 3)

  ## with 'minoverlap'
  result <- findOverlaps(query, subject, minoverlap = 3L)
  checkOverlap(result, integer(0), integer(0), 3, 3)
  result <- findOverlaps(query, subject, minoverlap = 2L)
  checkOverlap(result, 1, 2, 3, 3)
  result <- findOverlaps(query, subject, minoverlap = 2L, select = "first")
  checkIdentical(result, c(2L, NA, NA))
  result <- findOverlaps(query, subject, minoverlap = 2L, select = "last")
  checkIdentical(result, c(2L, NA, NA))
  result <- findOverlaps(query, subject, minoverlap = 2L, select = "arbitrary")
  checkIdentical(result, c(2L, NA, NA))

  ## zero-width ranges
  query <- IRanges(9:14, 8:13)
  result <- findOverlaps(query, subject, minoverlap = 1L)
  checkOverlap(result, integer(0), integer(0), 6, 3)
  result <- findOverlaps(query, subject)
  checkOverlap(result, c(3, 4), c(3, 3), 6, 3)
  result <- findOverlaps(query, subject, maxgap = 0L)
  checkOverlap(result, 2:5, c(3, 3, 3, 3), 6, 3)
  result <- findOverlaps(query, subject, maxgap = 1L)
  checkOverlap(result, 1:6, c(3, 3, 3, 3, 3, 3), 6, 3)
  result <- findOverlaps(subject, query, minoverlap = 1L)
  checkOverlap(result, integer(0), integer(0), 3, 6)
  result <- findOverlaps(subject, query)
  checkOverlap(result, c(3, 3), c(3, 4), 3, 6)
  result <- findOverlaps(subject, query, maxgap = 0L)
  checkOverlap(result, c(3, 3, 3, 3), 2:5, 3, 6)
  result <- findOverlaps(subject, query, maxgap = 1L)
  checkOverlap(result, c(3, 3, 3, 3, 3, 3), 1:6, 3, 6)

  ## .....
  ##    ....
  ##         ..
  ##  xxxx
  ##  xxx
  query <- IRanges(c(1, 4, 9), c(5, 7, 10))
  subject <- IRanges(c(2, 2), c(5, 4))

  result <- findOverlaps(query, subject)
  checkOverlap(result, c(1, 1, 2, 2), c(1, 2, 1, 2), 3, 2)

  result <- findOverlaps(subject, query)
  checkOverlap(result, c(1, 1, 2, 2), c(1, 2, 1, 2), 2, 3)

  query <- IRanges(c(1, 4, 9, 11), c(5, 7, 10, 11))

  result <- findOverlaps(query)
  checkOverlap(result, c(1, 1, 2, 2, 3, 4), c(1, 2, 1, 2, 3, 4), 4, 4, TRUE)

  ## check case of identical subjects
  ## .....
  ##    .....
  ##         ..
  ##  xxxx
  ##  xxxx
  ##      xx
  ##      xxx
  ##      xx
  query <- IRanges(c(1, 4, 9), c(5, 7, 10))
  subject <- IRanges(c(2, 2, 6, 6, 6), c(5, 5, 7, 8, 7))
  result <- findOverlaps(query, subject)
  checkOverlap(result, c(1, 1, 2, 2, 2, 2, 2), c(1, 2, 1, 2, 3, 4, 5), 3, 5)

  subject <- IRanges(c(1, 6, 13), c(4, 9, 14)) # single points
  checkIdentical(findOverlaps(c(3L, 7L, 10L), subject, select = "first"),
                 c(1L, 2L, NA))
  checkIdentical(findOverlaps(c(3L, 7L, 10L), subject, select = "last"),
                 c(1L, 2L, NA))
  checkIdentical(findOverlaps(c(3L, 7L, 10L), subject, select = "arbitrary"),
                 c(1L, 2L, NA))
  checkIdentical(findOverlaps(IRanges(c(2,1),c(3,4)), subject),
                 Hits(1:2, c(1, 1), 2, 3, sort.by.query=TRUE))

  ## check other types of matching

  ## ..
  ##     ..
  ##   ....  
  ##    ......
  ## xxxx
  ##   xxxx
  ##     xxxxx
  ##      xxxx

  query <- IRanges(c(1, 5, 3, 4), width=c(2, 2, 4, 6))
  subject <- IRanges(c(1, 3, 5, 6), width=c(4, 4, 5, 4))

  ## 'start'
  result <- findOverlaps(query, subject, type = "start")
  checkOverlap(result, c(1, 2, 3), c(1, 3, 2), 4, 4)

  ## minoverlap > 1L
  result <- findOverlaps(query, subject, type = "start", minoverlap = 3L)
  checkOverlap(result, 3, 2, 4, 4)

  ## 'end'
  result <- findOverlaps(query, subject, type = "end")
  checkOverlap(result, c(2, 3, 4, 4), c(2, 2, 3, 4), 4, 4)
  result <- findOverlaps(subject, query, type = "end")
  checkOverlap(result, c(2, 2, 3, 4), c(2, 3, 4, 4), 4, 4)

  ## select = "first"
  result <- findOverlaps(query, subject, type = "end", select = "first")
  checkIdentical(result, c(NA, 2L, 2L, 3L))

  ## 'within'
  result <- findOverlaps(query, subject, type = "within")
  checkOverlap(result, c(1, 2, 2, 3), c(1, 2, 3, 2), 4, 4)

  ## 'equal'
  result <- findOverlaps(query, subject, type = "equal")
  checkOverlap(result, 3, 2, 4, 4)

  checkException(findOverlaps(query, NULL), silent = TRUE)
  checkException(findOverlaps(NULL, query), silent = TRUE)
}

test_subsetByOverlaps_IntegerRanges <- function() {
  x <- IRanges(9:12, 15)
  ranges <- IRanges(1, 10)
  checkIdentical(x[1:2], subsetByOverlaps(x, ranges))
  checkIdentical(x[3:4], subsetByOverlaps(x, ranges, invert=TRUE))
  checkIdentical(x[1:3], subsetByOverlaps(x, ranges, maxgap=0))
  checkIdentical(x[4], subsetByOverlaps(x, ranges, maxgap=0, invert=TRUE))

  x <- IRanges(c(1, 4, 9), c(5, 7, 10))
  ranges <- IRanges(c(6, 8, 10), c(7, 12, 14))
  checkIdentical(x[2:3], subsetByOverlaps(x, ranges))
  checkIdentical(x[1], subsetByOverlaps(x, ranges, invert=TRUE))
  checkIdentical(x, subsetByOverlaps(x, ranges, maxgap=0))
  checkIdentical(x[0], subsetByOverlaps(x, ranges, maxgap=0, invert=TRUE))
}

