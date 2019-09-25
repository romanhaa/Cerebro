test_tile <- function() {
  ir <- IRanges()
  checkIdentical(tile(ir, n=3), IRangesList())
  checkIdentical(tile(ir, width=2), IRangesList())
  checkIdentical(tile(ir, n=0), IRangesList())
  ir <- IRanges(1, 4)
  checkIdentical(tile(ir, n=2), IRangesList(IRanges(c(1, 3), c(2, 4))))
  checkIdentical(tile(ir, n=2), tile(ir, width=2))
  ir <- IRanges(1, 5)
  checkIdentical(tile(ir, n=3), IRangesList(IRanges(c(1, 2, 4), c(1, 3, 5))))
  checkIdentical(tile(ir, n=3), tile(ir, width=2))
  ir <- IRanges(1, 4)
  checkIdentical(tile(ir, n=3), IRangesList(IRanges(1:3, c(1, 2, 4))))
  ir <- IRanges(1:3, width=5:3)
  checkIdentical(tile(ir, n=3),
                 IRangesList(IRanges(c(1, 2, 4), c(1, 3, 5)),
                             IRanges(c(2, 3, 4), c(2, 3, 5)),
                             IRanges(c(3, 4, 5), c(3, 4, 5))))
  checkIdentical(tile(ir, width=2),
                 IRangesList(IRanges(c(1, 2, 4), c(1, 3, 5)),
                             IRanges(c(2, 4), c(3, 5)),
                             IRanges(c(3, 4), c(3, 5))))
  checkIdentical(elementNROWS(tile(ir, width=4)), c(2L, 1L, 1L))
  checkException(tile(ir, n=4), silent=TRUE)
  checkException(tile(ir, width=-1), silent=TRUE)
  checkException(tile(ir, n=-1), silent=TRUE)
  ir <- setNames(IRanges(1:3, width = 10), letters[1:3])
  checkIdentical(names(ir), names(tile(ir, n = 2)))
  checkIdentical(names(ir), names(tile(ir, width = 3)))
}

test_slidingWindows <- function() {
    ir <- IRanges()
    checkIdentical(slidingWindows(ir, width=3), IRangesList())
    
    ir <- IRanges(1:3, width=5:3)
    checkIdentical(slidingWindows(ir, width=3, step=2),
                   IRangesList(IRanges(c(1, 3), c(3, 5)),
                               IRanges(c(2, 4), c(4, 5)),
                               IRanges(3, 5)))
    ir <- setNames(IRanges(1:3, width = 10), letters[1:3])
    checkIdentical(names(ir), names(slidingWindows(ir, width = 3)))
}
