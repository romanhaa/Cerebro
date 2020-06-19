checkMatching <- function(a, q, s, r, c) 
{
  mat <- cbind(queryHits = as.integer(q), subjectHits = as.integer(s))
  checkIdentical(as.matrix(a), mat)
  checkIdentical(c(queryLength(a), subjectLength(a)), as.integer(c(r, c)))
}

test_precede_follow_IntegerRanges <- function() 
{
  query <- IRanges(c(1, 3, 9), c(3, 7, 10))
  subject <- IRanges(c(3, 10, 2), c(3, 12, 5))
 
  checkIdentical(precede(query, subject), c(2L, 2L, NA))
  checkIdentical(precede(IRanges(), subject), integer())
  checkIdentical(precede(query, IRanges()), rep(NA_integer_, 3))
  checkIdentical(precede(query), c(3L, 3L, NA))
 
  checkIdentical(follow(query, subject), c(NA, NA, 3L))
  checkIdentical(follow(IRanges(), subject), integer())
  checkIdentical(follow(query, IRanges()), rep(NA_integer_, 3))
  checkIdentical(follow(query), c(NA, NA, 2L))

  checkMatching(precede(query, subject, select="all"),
                c(1, 2), c(2, 2), 3, 3)

  ## xxxx          
  ##  xxx
  ##         xx
  ##               xx
  ##               xxx
  ##      ..
  ##           ..
  ## ..        
  ##             ..
  ##                  ..

  subject <- IRanges(c(1, 2, 9, 15, 15), width=c(4, 3, 2, 2, 3))
  query <- IRanges(c(6, 11, 1, 13, 18), width=c(2, 2, 2, 2, 2))

  checkMatching(precede(query, subject, select="all"),
                c(1, 2, 2, 3, 4, 4), c(3, 4, 5, 3, 4, 5), 5, 5)
  checkMatching(precede(subject, query, select="all"),
                c(1, 2, 3, 4, 5), c(1, 1, 2, 5, 5), 5, 5)

  checkMatching(follow(query, subject, select="all"),
                c(1, 1, 2, 4, 5), c(1, 2, 3, 3, 5), 5, 5)
  checkMatching(follow(subject, query, select="all"),
                c(3, 4, 5), c(1, 4, 4), 5, 5)

  checkMatching(precede(query, select="all"),
                c(1, 2, 3, 4), c(2, 4, 1, 5), 5, 5)
  checkMatching(precede(subject, select="all"),
                c(1, 2, 3, 3), c(3, 3, 4, 5), 5, 5)

  checkMatching(follow(query, select="all"),
                c(1, 2, 4, 5), c(3, 1, 2, 4), 5, 5)
  checkMatching(follow(subject, select="all"),
                c(3, 3, 4, 5), c(1, 2, 3, 3), 5, 5)
}

test_nearest_IntegerRanges <- function() 
{
  query <- IRanges(c(1, 3, 9), c(2, 7, 10))
  subject <- IRanges(c(3, 5, 12), c(3, 6, 12))

  ## 2 possible results
  current <- nearest(query, subject)
  target1 <- c(1L, 1L, 3L)
  target2 <- c(1L, 2L, 3L)
  checkTrue(identical(target1, current) || identical(target2, current))
  checkIdentical(nearest(query), c(2L, 1L, 2L))
  checkIdentical(nearest(query, subject[c(2,3,1)]), c(3L, 3L, 2L))

  ## xxxx          
  ##  xxx
  ##         xx
  ##               xx
  ##               xxx
  ##      ..
  ##           ..
  ## ..        
  ##             ..
  ##                  ..

  subject <- IRanges(c(1, 2, 9, 15, 15), width=c(4, 3, 2, 2, 3))
  query <- IRanges(c(6, 11, 1, 13, 18), width=c(2, 2, 2, 2, 2))
 
  checkMatching(nearest(query, subject, select = "all"),
                c(1, 1, 1, 2, 3, 3, 4, 4, 5),
                c(1, 2, 3, 3, 1, 2, 4, 5, 5), 5, 5)
  checkMatching(nearest(subject, query, select = "all"),
                c(1, 2, 3, 4, 5, 5), c(3, 3, 2, 4, 4, 5), 5, 5)
 
  checkMatching(nearest(subject, select="all"),
                c(1, 2, 3, 3, 3, 3, 4, 5), c(2, 1, 1, 2, 4, 5, 5, 4), 5, 5)
  checkMatching(nearest(query, select="all"),
                c(1, 1, 2, 3, 4, 5), c(2, 3, 4, 1, 2, 4), 5, 5)
}

quiet <- suppressWarnings
test_distance_IntegerRanges <- function() 
{
  checkIdentical(quiet(distance(IRanges(), IRanges())), integer())

  ## adjacent, overlap, separated by 1
  query <- IRanges(c(1, 3, 9), c(2, 7, 10))
  subject <- IRanges(c(3, 5, 12), c(3, 6, 12))
  checkIdentical(quiet(distance(query, subject)), c(0L, 0L, 1L))
  ## recycling
  checkIdentical(quiet(distance(query[1:2], subject)), 
                 c(0L, 0L, 9L))
  ## zero-width
  target <- abs(-3:3)
  current <- sapply(-3:3, function(i) 
                 quiet(distance(shift(IRanges(4,3), i), IRanges(4,3))))
  checkIdentical(target, current)
  checkIdentical(quiet(distance(IRanges(4,3), IRanges(3,4))), 0L)
}

test_distanceToNearest_IntegerRanges <- function() 
{
  target <- Hits(sort.by.query=TRUE)
  current <- distanceToNearest(IRanges(), IRanges())
  checkIdentical(queryHits(current), queryHits(target))
  checkIdentical(subjectHits(current), subjectHits(target))
  checkIdentical(queryLength(current), queryLength(target))

  x <- IRanges(5, 10)
  subject <- IRanges(c(1, 1, 1), c(4, 5, 6))
  current <- distanceToNearest(x, subject, select="all")
  checkIdentical(subjectHits(current), 1:3)

  current <- distanceToNearest(x, rev(subject), select="all")
  checkIdentical(subjectHits(current), 1:3)

  current <- distanceToNearest(x, IRanges())
  checkIdentical(length(current), 0L)
  checkIdentical(queryLength(current), 1L)
  checkIdentical(subjectLength(current), 0L)

  x <- IRanges(c(2, 4, 12, 15), c(2, 3, 13, 14))
  subject <- IRanges(1, 10)
  current <- distanceToNearest(x, subject)
  checkIdentical(queryHits(current), 1:4)
  checkIdentical(mcols(current)$distance, c(0L, 0L, 1L, 4L))
}

