test_DataFrameList_construction <- function() {
    checkDFL2dfl <- function(DFL, dfl) {
        checkIdentical(lapply(as.list(DFL), as.data.frame), dfl)
    }
    data(airquality)
    data(swiss)
    checkDFL2dfl(DataFrameList(swiss, airquality), list(swiss, airquality))
}

test_SplitDataFrameList_construction <- function() {
    checkDFL2dfl <- function(DFL, dfl) {
        checkIdentical(lapply(as.list(DFL), as.data.frame), dfl)
    }
    striprownames <- function(x) {
        lapply(x, function(y) {
                   rownames(y) <- NULL
                   y
               })
    }
    data(airquality)
    data(swiss)
    aq <- DataFrame(airquality)
    sw <- DataFrame(swiss, row.names=rownames(swiss))

    aqsplit1 <- split(aq, aq[["Month"]])
    aqsplit2 <-
      SplitDataFrameList(lapply(split(airquality, airquality[["Month"]]),
                                as, "DataFrame"))
    checkIdentical(aqsplit1, aqsplit2)

    swsplit1 <- split(sw, sw[["Education"]])
    swsplit2 <-
      SplitDataFrameList(lapply(split(swiss, swiss[["Education"]]),
                                as, "DataFrame"))
    checkIdentical(swsplit1, swsplit2)

    for (compress in c(TRUE, FALSE)) {
        airqualitysplit <-
          striprownames(split(airquality, airquality[["Month"]]))
        aqsplit <-
          SplitDataFrameList(as.list(split(aq, aq[["Month"]])),
                             compress = compress)
        checkDFL2dfl(aqsplit, airqualitysplit)

        swisssplit <- split(swiss, swiss[["Education"]])
        swsplit <-
          SplitDataFrameList(as.list(split(sw, sw[["Education"]])),
                             compress = compress)
        checkDFL2dfl(swsplit, swisssplit)
    }
}

test_DataFrameList_subset <- function() {
    checkDFL2dfl <- function(DFL, dfl) {
        checkIdentical(lapply(as.list(DFL), as.data.frame), dfl)
    }
    data(airquality)
    data(swiss)

    DFL1 <- DataFrameList(swiss, airquality)
    dfl1 <- list(swiss, airquality)
    checkDFL2dfl(DFL1[], dfl1[])
    checkDFL2dfl(DFL1[1], dfl1[1])
    checkDFL2dfl(DFL1[2:1], dfl1[2:1])
    checkIdentical(as.data.frame(DFL1[[2]]), airquality)
    checkException(DFL1[[3]], silent = TRUE)

    DFL2 <- DataFrameList(s = swiss, a = airquality)
    dfl2 <- list(s = swiss, a = airquality)
    checkDFL2dfl(DFL2[], dfl2[])
    checkDFL2dfl(DFL2[1], dfl2[1])
    checkDFL2dfl(DFL2["a"], dfl2["a"])
    checkDFL2dfl(DFL2[c("a", "s")], dfl2[c("a", "s")])
    checkIdentical(as.data.frame(DFL2[["a"]]), airquality)
    checkIdentical(DFL2[["z"]], NULL)
}

test_SplitDataFrameList_subset <- function() {
    checkDFL2dfl <- function(DFL, dfl) {
        checkIdentical(lapply(as.list(DFL), as.data.frame), dfl)
    }
    data(swiss)
    sw <- DataFrame(swiss, row.names = rownames(swiss))

    for (compress in c(TRUE, FALSE)) {
        swsplit <-
          SplitDataFrameList(as.list(split(sw, sw[["Education"]])),
                             compress = compress)
        swisssplit <- split(swiss, swiss[["Education"]])

        checkDFL2dfl(swsplit[], swisssplit[])
        checkDFL2dfl(swsplit[1], swisssplit[1])
        checkDFL2dfl(swsplit[2:1], swisssplit[2:1])
        checkIdentical(as.data.frame(swsplit[[2]]), swisssplit[[2]])
        checkIdentical(swsplit[["A"]], NULL)
        checkException(swsplit[[30]], silent = TRUE)

        checkIdentical(as.list(swsplit[,1]),
                       split(swiss[[1]], swiss[["Education"]]))
        checkIdentical(as.list(swsplit[,"Examination"]),
                       split(swiss[["Examination"]], swiss[["Education"]]))
    }
}

test_SplitDataFrameList_as.data.frame <- function() {
    checkDFL2dfl <- function(DFL, dfl, compress) {
        target <- 
          data.frame(group = togroup(PartitioningByWidth(dfl)),
                     group_name = names(dfl)[togroup(PartitioningByWidth(dfl))],
                     do.call(rbind, dfl),
                     stringsAsFactors=FALSE, row.names=NULL)
        rownames(target) <- unlist(lapply(dfl, row.names), use.names = FALSE)
        checkIdentical(target, as.data.frame(DFL))
    }

    data(swiss)
    sw <- DataFrame(swiss, row.names = rownames(swiss))

    for (compress in c(TRUE, FALSE)) {
        swsplit <-
          SplitDataFrameList(as.list(split(sw, sw[["Education"]])),
                             compress = compress)
        swisssplit <- split(swiss, swiss[["Education"]])
        checkDFL2dfl(swsplit, swisssplit, compress)
    }
}

test_DataFrameList_replace <- function() {
    checkDFL2dfl <- function(DFL, dfl) {
        checkIdentical(lapply(as.list(DFL), as.data.frame), dfl)
    }
    data(airquality)
    data(swiss)

    DFL1 <- DataFrameList(swiss, airquality)
    dfl1 <- list(swiss, airquality)
    DFL1[] <- DFL1[1]
    dfl1[] <- dfl1[1]
    checkDFL2dfl(DFL1, dfl1)

    DFL1 <- DataFrameList(swiss, airquality)
    dfl1 <- list(swiss, airquality)
    DFL1[2] <- DFL1[1]
    dfl1[2] <- dfl1[1]
    checkDFL2dfl(DFL1, dfl1)

    DFL1 <- DataFrameList(swiss, airquality)
    dfl1 <- list(swiss, airquality)
    DFL1[[1]][[1]] <- DFL1[[1]][[1]] + 1L
    dfl1[[1]][[1]] <- dfl1[[1]][[1]] + 1L
    checkDFL2dfl(DFL1, dfl1)
}

test_SplitDataFrameList_replace <- function() {
    checkDFL2dfl <- function(DFL, dfl) {
        checkIdentical(lapply(as.list(DFL), as.data.frame), dfl)
    }
    striprownames <- function(x) {
        lapply(x, function(y) {
                   rownames(y) <- NULL
                   y
               })
    }
    data(airquality)
    data(swiss)
    swiss2 <- swiss
    rownames(swiss2) <- NULL
    sw2 <- DataFrame(swiss2)
    for (compress in c(TRUE, FALSE)) {
        swiss2split <- striprownames(split(swiss2, swiss2[["Education"]]))
        sw2split <-
          SplitDataFrameList(as.list(split(sw2, sw2[["Education"]])),
                             compress = compress)
        swiss2split[] <- swiss2split[1]
        sw2split[] <- sw2split[1]
        checkDFL2dfl(sw2split, swiss2split)

        swiss2split <- striprownames(split(swiss2, swiss2[["Education"]]))
        sw2split <-
          SplitDataFrameList(as.list(split(sw2, sw2[["Education"]])),
                             compress = compress)
        swiss2split[c(2, 4, 5)] <- swiss2split[1]
        sw2split[c(2, 4, 5)] <- sw2split[1]
        checkDFL2dfl(sw2split, swiss2split)

        swiss2split <- striprownames(split(swiss2, swiss2[["Education"]]))
        swiss2split <-
          lapply(swiss2split,
                 function(x) {x[["Examination"]] <- x[["Examination"]] + 1L; x})
        sw2split <-
          SplitDataFrameList(as.list(split(sw2, sw2[["Education"]])),
                             compress = compress)
        sw2split[,"Examination"] <- sw2split[,"Examination"] + 1L
        checkDFL2dfl(sw2split, swiss2split)

        swiss2split <- striprownames(split(swiss2, swiss2[["Education"]]))
        swiss2split <-
          lapply(swiss2split, function(x) {
                     x[["Examination"]][x[["Examination"]] > 22] <-
                       x[["Examination"]][x[["Examination"]] > 22] + 1L
                     x
                 })
        sw2split <-
          SplitDataFrameList(as.list(split(sw2, sw2[["Education"]])),
                             compress = compress)
        sw2split[sw2split[, "Examination"] > 22, "Examination"] <-
          sw2split[sw2split[, "Examination"] > 22,"Examination"] + 1L
        checkDFL2dfl(sw2split, swiss2split)
    }
}

test_DataFrameList_transform <- function() {
  DF <- DataFrame(state.division, state.region, state.area)
  DFL <- split(DF, DF$state.division) # NICER: split(DF, ~ state.devision)
  DFL <- transform(DFL, total.area=sum(state.area[state.region!="South"]),
                   fraction=ifelse2(total.area == 0, 0, state.area/total.area))

  ANS <- DataFrame(lapply(unlist(DFL, use.names=FALSE), unname))
  
  df <- as.data.frame(DF)
  df$total.area <-
    with(subset(df, state.region != "South"),
         sapply(split(state.area, state.division), sum))[df$state.division]
  df$fraction <- with(df, ifelse(total.area == 0, 0, state.area/total.area))
  df <- df[order(df$state.division),]
  rownames(df) <- NULL
  
  checkIdentical(ANS, DataFrame(df))
}
