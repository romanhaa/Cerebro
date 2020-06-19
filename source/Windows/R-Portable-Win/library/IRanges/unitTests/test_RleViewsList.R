test_RleViewsList <- function() {
    x1 <- rep(c(1L, 3L, NA, 7L, 9L), 1:5)
    x1Rle <- Rle(x1)
    x1Ranges <- IRanges(start = c(1, 3, 5, 7, 9), end = c(1, 13, 11, 10, 9))
    x2 <- rev(x1)
    x2Rle <- Rle(x2)
    x2Ranges <- IRanges(start = c(2, 4, 6, 8, 10), end = c(3, 9, 11, 13, 15))
    checkIdentical(RleViewsList(Views(x1Rle, x1Ranges), Views(x2Rle, x2Ranges)),
                   RleViewsList(rleList = RleList(x1Rle, x2Rle),
                                rangesList = IRangesList(x1Ranges, x2Ranges)))
    xRleViewsList <- RleViewsList(a = Views(x1Rle, x1Ranges), b = Views(x2Rle, x2Ranges))
    xList <-
      list(a = lapply(seq_len(length(xRleViewsList[[1]])),
                      function(i) window(x1, start = start(x1Ranges)[i],
                                         end = end(x1Ranges)[i])),
           b = lapply(seq_len(length(xRleViewsList[[2]])),
                      function(i) window(x2, start = start(x2Ranges)[i],
                                         end = end(x2Ranges)[i])))
    checkIdentical(c("a", "b"), names(viewApply(xRleViewsList, min)))
    checkIdentical(c("a", "b"), names(viewMins(xRleViewsList)))
    checkIdentical(c("a", "b"), names(viewMaxs(xRleViewsList)))
    checkIdentical(c("a", "b"), names(viewSums(xRleViewsList)))
    checkIdentical(c("a", "b"), names(viewMeans(xRleViewsList)))
    checkIdentical(c("a", "b"), names(viewWhichMins(xRleViewsList)))
    checkIdentical(c("a", "b"), names(viewWhichMaxs(xRleViewsList)))
    checkIdentical(c("a", "b"), names(viewRangeMins(xRleViewsList, na.rm = TRUE)))
    checkIdentical(c("a", "b"), names(viewRangeMaxs(xRleViewsList, na.rm = TRUE)))

    checkEqualsNumeric(unlist(lapply(xList, lapply, min)), unlist(viewMins(xRleViewsList)))
    checkEqualsNumeric(unlist(lapply(xList, lapply, min)), unlist(viewApply(xRleViewsList, min)))
    checkEqualsNumeric(unlist(lapply(xList, lapply, min, na.rm = TRUE)), unlist(viewMins(xRleViewsList, na.rm = TRUE)))
    checkEqualsNumeric(unlist(lapply(xList, lapply, min, na.rm = TRUE)), unlist(viewApply(xRleViewsList, min, na.rm = TRUE)))
    checkEqualsNumeric(unlist(lapply(xList, lapply, max)), unlist(viewMaxs(xRleViewsList)))
    checkEqualsNumeric(unlist(lapply(xList, lapply, max)), unlist(viewApply(xRleViewsList, max)))
    checkEqualsNumeric(unlist(lapply(xList, lapply, max, na.rm = TRUE)), unlist(viewMaxs(xRleViewsList, na.rm = TRUE)))
    checkEqualsNumeric(unlist(lapply(xList, lapply, max, na.rm = TRUE)), unlist(viewApply(xRleViewsList, max, na.rm = TRUE)))
    checkEqualsNumeric(unlist(lapply(xList, lapply, sum)), unlist(viewSums(xRleViewsList)))
    checkEqualsNumeric(unlist(lapply(xList, lapply, mean)), unlist(viewMeans(xRleViewsList)))
    checkEqualsNumeric(unlist(lapply(xList, lapply, sum)), unlist(viewApply(xRleViewsList, sum)))
    checkEqualsNumeric(unlist(lapply(xList, lapply, sum, na.rm = TRUE)), unlist(viewSums(xRleViewsList, na.rm = TRUE)))
    checkEqualsNumeric(unlist(lapply(xList, lapply, mean, na.rm = TRUE)), unlist(viewMeans(xRleViewsList, na.rm = TRUE)))
    checkEqualsNumeric(unlist(lapply(xList, lapply, sum, na.rm = TRUE)), unlist(viewApply(xRleViewsList, sum, na.rm = TRUE)))

    y1 <- rep(c(1.2, 3.4, NA, 7.8, 9.0), 1:5)
    y1Ranges <- IRanges(start = c(1, 3, 5, 7, 9), end = c(1, 13, 11, 10, 9))
    y1Rle <- Rle(y1)
    y2 <- rev(y1)
    y2Rle <- Rle(y2)
    y2Ranges <- IRanges(start = c(2, 4, 6, 8, 10), end = c(3, 9, 11, 13, 15))
    checkIdentical(RleViewsList(Views(y1Rle, y1Ranges), Views(y2Rle, y2Ranges)),
                   RleViewsList(rleList = RleList(y1Rle, y2Rle),
                                rangesList = IRangesList(y1Ranges, y2Ranges)))
    yRleViewsList <- RleViewsList(Views(y1Rle, y1Ranges), Views(y2Rle, y2Ranges))
    yList <-
      list(lapply(seq_len(length(yRleViewsList[[1]])),
                  function(i) window(y1, start = start(y1Ranges)[i],
                                     end = end(y1Ranges)[i])),
           lapply(seq_len(length(yRleViewsList[[2]])),
                  function(i) window(y2, start = start(y2Ranges)[i],
                                     end = end(y2Ranges)[i])))
    checkEqualsNumeric(unlist(lapply(yList, lapply, min)), unlist(viewMins(yRleViewsList)))
    checkEqualsNumeric(unlist(lapply(yList, lapply, min)), unlist(viewApply(yRleViewsList, min)))
    checkEqualsNumeric(unlist(lapply(yList, lapply, min, na.rm = TRUE)), unlist(viewMins(yRleViewsList, na.rm = TRUE)))
    checkEqualsNumeric(unlist(lapply(yList, lapply, min, na.rm = TRUE)), unlist(viewApply(yRleViewsList, min, na.rm = TRUE)))
    checkEqualsNumeric(unlist(lapply(yList, lapply, max)), unlist(viewMaxs(yRleViewsList)))
    checkEqualsNumeric(unlist(lapply(yList, lapply, max)), unlist(viewApply(yRleViewsList, max)))
    checkEqualsNumeric(unlist(lapply(yList, lapply, max, na.rm = TRUE)), unlist(viewMaxs(yRleViewsList, na.rm = TRUE)))
    checkEqualsNumeric(unlist(lapply(yList, lapply, max, na.rm = TRUE)), unlist(viewApply(yRleViewsList, max, na.rm = TRUE)))
    checkEqualsNumeric(unlist(lapply(yList, lapply, sum)), unlist(viewSums(yRleViewsList)))
    checkEqualsNumeric(unlist(lapply(yList, lapply, mean)), unlist(viewMeans(yRleViewsList)))
    checkEqualsNumeric(unlist(lapply(yList, lapply, sum)), unlist(viewApply(yRleViewsList, sum)))
    checkEqualsNumeric(unlist(lapply(yList, lapply, sum, na.rm = TRUE)), unlist(viewSums(yRleViewsList, na.rm = TRUE)))
    checkEqualsNumeric(unlist(lapply(yList, lapply, mean, na.rm = TRUE)), unlist(viewMeans(yRleViewsList, na.rm = TRUE)))
    checkEqualsNumeric(unlist(lapply(yList, lapply, sum, na.rm = TRUE)), unlist(viewApply(yRleViewsList, sum, na.rm = TRUE)))
}
