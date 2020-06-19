test_IRanges_names <- function() {
  range1 <- IRanges(start=c(1,2,3), end=c(5,2,8))
  
  checkIdentical(names(range1), NULL)
  nms <- c("a", NA, "b")
  names(range1) <- nms
  checkIdentical(names(range1), nms)
  checkTrue(validObject(nms))
  names(range1) <- NULL
  checkTrue(validObject(nms))
  checkIdentical(names(range1), NULL)
  names(range1) <- "a"
  checkTrue(validObject(range1))
  checkIdentical(names(range1), c("a", NA, NA))

  checkException(names(range1) <- c("a", "b", "c", "d"), silent = TRUE)
}

test_IntegerRanges_isDisjoint <- function() {
  ir1 <- IRanges(c(2,5,1), c(3,7,3))
  ir2 <- IRanges(c(2,9,5), c(3,9,6))
  ir3 <- IRanges(1, 5)
  checkIdentical(isDisjoint(ir1), FALSE)
  checkIdentical(isDisjoint(ir2), TRUE)
  checkIdentical(isDisjoint(ir3), TRUE)

  ## Handling of zero-width ranges
  current <- sapply(11:17,
                    function(i)
                      isDisjoint(IRanges(c(12, i), width=c(4, 0))))
  target <- rep(c(TRUE, FALSE, TRUE), c(2, 3, 2))
  checkIdentical(target, current)
}

test_IRanges_concatenate <- function() {
  range <- IRanges(start=c(1,2,3,1), end=c(5,2,8,3))
  srange <- split(range, start(range) == 1)
  checkIdentical(srange, IRangesList(`FALSE` = range[2:3], `TRUE` = range[c(1,4)]))
  checkIdentical(do.call(c, unname(as.list(srange))),
                 IRanges(c(2,3,1,1), c(2,8,5,3)))

  ir1 <- IRanges(1, 10)
  ir2 <- IRanges(c(1, 15), width=5)
  mcols(ir2) <- DataFrame(score=1:2)
  checkIdentical(mcols(c(ir1, ir2)),
                 DataFrame(score = c(NA, 1L, 2L)))
  
  ## Concatenating multiple IRanges object with varying mcols
  mcols(ir1) <- DataFrame(gc=0.78)
  ir12 <- c(ir1, ir2, ignore.mcols=TRUE)
  checkIdentical(mcols(ir12), NULL)
  target_mcols <- DataFrame(gc=c(0.78, NA, NA), score=c(NA, 1:2))
  mcols(ir12) <- target_mcols
  checkIdentical(c(ir1, ir2), ir12)
}

test_IRanges_annotation <- function() {
  range <- IRanges(c(1, 4), c(5, 7))
  mcols(range) <- DataFrame(a = 1:2)
  checkIdentical(mcols(range)[,1], 1:2)
  checkIdentical(mcols(range[2:1])[,1], 2:1)
  checkIdentical(mcols(c(range,range))[,1], rep(1:2,2))
}

