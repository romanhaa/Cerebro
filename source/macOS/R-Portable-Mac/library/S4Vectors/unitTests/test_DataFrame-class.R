test_DataFrame_construction <- function() {
  score <- c(X=1L, Y=3L, Z=NA)
  counts <- c(10L, 2L, NA)

  ## na in rn
  checkException(DataFrame(score, row.names = c("a", NA, "b")), silent = TRUE)
  ## invalid rn length
  checkException(DataFrame(score, row.names = "a"), silent = TRUE)
 
  DF <- DataFrame() # no args
  checkTrue(validObject(DF))
  row.names <- c("one", "two", "three")
  DF <- DataFrame(row.names = row.names) # no args, but row.names
  checkTrue(validObject(DF))
  checkIdentical(rownames(DF), row.names)
 
  DF <- DataFrame(score) # single, unnamed arg
  checkTrue(validObject(DF))
  checkIdentical(DF[["score"]], score)
  DF <- DataFrame(score, row.names = row.names) # with row names
  checkTrue(validObject(DF))
  checkIdentical(rownames(DF), row.names)

  ## dups in rn
  row.names = c("a", "b", "a")
  DF <- DataFrame(score, row.names = row.names)
  checkTrue(validObject(DF))
  checkIdentical(rownames(DF), row.names)
  DF <- DataFrame(score=setNames(score, row.names))
  checkTrue(validObject(DF))
  checkIdentical(rownames(DF), row.names)
 
  DF <- DataFrame(vals = score) # named vector arg
  checkTrue(validObject(DF))
  checkIdentical(DF[["vals"]], score)
  DF <- DataFrame(counts, vals = score) # mixed named and unnamed
  checkTrue(validObject(DF))
  checkIdentical(DF[["vals"]], score)
  checkIdentical(DF[["counts"]], counts)
  DF <- DataFrame(score + score) # unnamed arg with invalid name expression
  checkTrue(validObject(DF))
  checkIdentical(DF[["score...score"]], score + score)

  mat <- cbind(score)
  DF <- DataFrame(mat) # single column matrix with column name
  checkTrue(validObject(DF))
  checkIdentical(DF[["score"]], unname(score))
  mat <- cbind(score, counts)
  DF <- DataFrame(mat) # two column matrix with col names
  checkTrue(validObject(DF))
  checkIdentical(DF[["score"]], unname(score))
  checkIdentical(DF[["counts"]], counts)
  colnames(mat) <- NULL
  DF <- DataFrame(mat) # two column matrix without col names
  checkTrue(validObject(DF))
  checkIdentical(DF[["V1"]], unname(score))

  sw <- DataFrame(swiss, row.names = rownames(swiss)) # a data.frame
  checkIdentical(as.data.frame(sw), swiss)
  rownames(swiss) <- NULL # strip row names to make them comparable
  sw <- DataFrame(swiss) # a data.frame
  checkIdentical(as.data.frame(sw), swiss)
  sw <- DataFrame(swiss[1:3,], score = unname(score))
  checkIdentical(as.data.frame(sw), data.frame(swiss[1:3,], score))
  sw <- DataFrame(score = score, swiss = swiss[1:3,]) # named data.frame/matrix
  checkIdentical(as.data.frame(sw),
                 data.frame(score = score, swiss = swiss[1:3,]))

  ## identity
  df <- DataFrame(A=I(list(1:3)))
  checkIdentical(as.data.frame(df), data.frame(A=I(list(1:3))))
  
  ## recycling
  DF <- DataFrame(1, score)
  checkIdentical(DF[[1]], rep(1, 3)) 
  checkIdentical(DF[[2]], score) 
}

test_DataFrame_coerce <- function() {
  ## need to introduce character() dim names
  checkTrue(validObject(as(matrix(0L, 0L, 0L), "DataFrame")))

  score <- c(X=1L, Y=3L, Z=NA)
  DF <- as(score, "DataFrame")
  checkTrue(validObject(DF))
  checkIdentical(DF[[1]], score)
}

test_DataFrame_subset <- function() {
  data(swiss)
  sw <- DataFrame(swiss)
  rn <- rownames(swiss)

  checkException(sw[list()], silent = TRUE) # non-atomic
  checkException(sw[NA], silent = TRUE) # column indices cannot be NA
  checkException(sw[100], silent = TRUE) # out of bounds col
  checkException(sw[,100], silent = TRUE)
  checkException(sw[1000,], silent = TRUE) # out of bounds row
  oldOpts <- options(warn=2)
  checkException(sw[1:3, drop=TRUE], silent = TRUE) # drop ignored
  checkException(sw[drop=TRUE], silent = TRUE)
  checkException(sw[foo = "bar"], silent = TRUE) # invalid argument
  options(oldOpts)
  checkException(sw[,"Fert"], silent = TRUE) # bad column name

  sw <- DataFrame(swiss)

  checkIdentical(sw[], sw) # identity subset
  checkIdentical(sw[,], sw)

  checkIdentical(sw[NULL], DataFrame(swiss[NULL])) # NULL subsetting
  checkIdentical(sw[,NULL], DataFrame(swiss[,NULL]))
  checkIdentical(as.data.frame(sw[NULL,]),
                 structure(data.frame(swiss[NULL,]), row.names = character()))

  rownames(sw) <- rn

  ## select columns
  checkIdentical(as.data.frame(sw[1:3]), swiss[1:3])
  checkIdentical(as.data.frame(sw[, 1:3]), swiss[1:3])
  ## select rows
  checkIdentical(as.data.frame(sw[1:3,]), swiss[1:3,])
  checkIdentical(as.data.frame(sw[1:3,]), swiss[1:3,])
  checkIdentical(as.data.frame(sw[sw[["Education"]] == 7,]),
                 swiss[swiss[["Education"]] == 7,])
  checkIdentical(as.data.frame(sw[Rle(sw[["Education"]] == 7),]),
                 swiss[swiss[["Education"]] == 7,])
  ## select rows and columns
  checkIdentical(as.data.frame(sw[4:5, 1:3]), swiss[4:5,1:3])

  checkIdentical(as.data.frame(sw[1]), swiss[1])  # a one-column data frame
  checkIdentical(sw[,"Fertility"], swiss[,"Fertility"])
  ## the same
  checkIdentical(as.data.frame(sw[, 1, drop = FALSE]), swiss[, 1, drop = FALSE])
  checkIdentical(sw[, 1], swiss[,1])      # a (unnamed) vector
  checkIdentical(sw[[1]], swiss[[1]])      # the same
  checkIdentical(sw[["Fertility"]], swiss[["Fertility"]])
  checkIdentical(sw[["Fert"]], swiss[["Fert"]]) # should return 'NULL'
  checkIdentical(sw[,c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)],
                 swiss[,c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)])

  checkIdentical(as.data.frame(sw[1,]), swiss[1,])       # a one-row data frame
  checkIdentical(sw[1,, drop=TRUE], swiss[1,, drop=TRUE]) # a list

  ## duplicate row, unique row names are created
  checkIdentical(as.data.frame(sw[c(1, 1:2),]), swiss[c(1,1:2),])

  ## NOTE: NA subsetting not yet supported for XVectors
  ##checkIdentical(as.data.frame(sw[c(1, NA, 1:2, NA),]), # mixin some NAs
  ##               swiss[c(1, NA, 1:2, NA),])

  checkIdentical(as.data.frame(sw["Courtelary",]), swiss["Courtelary",])
  subswiss <- swiss[1:5,1:4]
  subsw <- sw[1:5,1:4]
  checkIdentical(as.data.frame(subsw["C",]), subswiss["C",]) # partially matches
  ## NOTE: NA subsetting not yet supported for XVectors
  ##checkIdentical(as.data.frame(subsw["foo",]), # bad row name
  ##               subswiss["foo",]) 
  ##checkIdentical(as.data.frame(sw[match("C", row.names(sw)), ]), 
  ##               swiss[match("C", row.names(sw)), ]) # no exact match
}

test_DataFrame_dimnames_replace <- function() {
  data(swiss)
  cn <- paste("X", seq_len(ncol(swiss)), sep = ".")
  sw <- DataFrame(swiss)
  colnames(sw) <- cn
  checkIdentical(colnames(sw), cn)
  cn <- as.character(seq_len(ncol(swiss)))
  colnames(sw) <- cn
  colnames(swiss) <- cn
  checkIdentical(colnames(sw), colnames(swiss))
  colnames(sw) <- cn[1]
  colnames(swiss) <- cn[1]
  checkIdentical(colnames(sw), colnames(swiss))
  rn <- seq(nrow(sw))
  rownames(sw) <- rn
  checkIdentical(rownames(sw), as.character(rn))
  checkException(rownames(sw) <- rn[1], silent = TRUE)
  rownames(sw) <- rep(rn[1], nrow(sw))
  checkIdentical(rownames(sw), as.character(rep(rn[1], nrow(sw))))
  rn[1] <- NA
  checkException(rownames(sw) <- rn, silent = TRUE)
}

test_DataFrame_replace <- function() {
  score <- c(1L, 3L, NA)
  counts <- c(10L, 2L, NA)

  DF <- DataFrame(score) # single, unnamed arg

  DF[["counts"]] <- counts
  checkIdentical(DF[["counts"]], counts)
  DF[[3]] <- score
  checkIdentical(DF[[3]], score)
  DF[[3]] <- NULL # deletion
  DF[["counts"]] <- NULL
  DF$counts <- counts
  checkIdentical(DF$counts, counts)

  checkException(DF[[13]] <- counts, silent = TRUE) # index must be < length+1
  checkException(DF[["tooshort"]] <- counts[1:2], silent = TRUE)

  sw <- DataFrame(swiss, row.names = rownames(swiss)) # a data.frame
  sw1 <- sw; swiss1 <- swiss
  sw1[] <- 1L; swiss1[] <- 1L
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[] <- 1; swiss1[] <- 1
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1["Education"] <- 1; swiss1["Education"] <- 1
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[,"Education"] <- 1; swiss1[,"Education"] <- 1
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1["Courtelary",] <- 1; swiss1["Courtelary",] <- 1
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[1:3] <- 1; swiss1[1:3] <- 1
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[,1:3] <- 1; swiss1[,1:3] <- 1
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[2:4,1:3] <- 1; swiss1[2:4,1:3] <- 1
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[2:4,-c(2,4,5)] <- 1; swiss1[2:4,-c(2,4,5)] <- 1
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[,1:3] <- sw1[,2:4]; swiss1[,1:3] <- swiss1[,2:4]
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[2:4,] <- sw1[1:3,]; swiss1[2:4,] <- swiss1[1:3,]
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[2:4,1:3] <- sw1[1:3,2:4]; swiss1[2:4,1:3] <- swiss1[1:3,2:4]
  checkIdentical(as.data.frame(sw1), swiss1)

  sw1 <- sw; swiss1 <- swiss
  sw1["NewCity",] <- NA; swiss1["NewCity",] <- NA
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[nrow(sw1)+(1:2),] <- NA; swiss1[nrow(swiss1)+(1:2),] <- NA
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1["NewCol"] <- seq(nrow(sw1)); swiss1["NewCol"] <- seq(nrow(sw1))
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[ncol(sw1)+1L] <- seq(nrow(sw1)); swiss1[ncol(swiss1)+1L] <- seq(nrow(sw1))
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[,"NewCol"] <- seq(nrow(sw1)); swiss1[,"NewCol"] <- seq(nrow(sw1))
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1["NewCity","NewCol"] <- 0
  swiss1["NewCity","NewCol"] <- 0
  checkIdentical(as.data.frame(sw1), swiss1)

  sw1 <- sw; swiss1 <- swiss
  sw1["NewCity",] <- DataFrame(NA); swiss1["NewCity",] <- data.frame(NA)
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[nrow(sw1)+(1:2),] <- DataFrame(NA)
  swiss1[nrow(swiss1)+(1:2),] <- data.frame(NA)
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1["NewCol"] <- DataFrame(seq(nrow(sw1)))
  swiss1["NewCol"] <- data.frame(seq(nrow(sw1)))
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[ncol(sw1)+1L] <- DataFrame(seq(nrow(sw1)))
  swiss1[ncol(swiss1)+1L] <- data.frame(seq(nrow(sw1)))
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1[,"NewCol"] <- DataFrame(seq(nrow(sw1)))
  swiss1[,"NewCol"] <- data.frame(seq(nrow(sw1)))
  checkIdentical(as.data.frame(sw1), swiss1)
  sw1 <- sw; swiss1 <- swiss
  sw1["NewCity","NewCol"] <- DataFrame(0)
  swiss1["NewCity","NewCol"] <- data.frame(0)
  checkIdentical(as.data.frame(sw1), swiss1)

  sw1 <- sw
  mcols(sw1) <- DataFrame(id = seq_len(ncol(sw1)))
  sw1["NewCol"] <- DataFrame(seq(nrow(sw1)))
  checkIdentical(mcols(sw1, use.names=TRUE),
                 DataFrame(id = c(seq_len(ncol(sw1)-1), NA),
                           row.names = colnames(sw1)))
}

## combining
test_DataFrame_combine <- function() {
  data(swiss)
  rn <- rownames(swiss)
  sw <- DataFrame(swiss, row.names=rn)
  swisssplit <- split(swiss, swiss$Education)
 
  ## rbind
  checkIdentical(rbind(DataFrame(), DataFrame()), DataFrame())
  score <- c(X=1L, Y=3L, Z=NA)
  DF <- DataFrame(score)
  checkIdentical(rbind(DF, DF)[[1]], c(score, score))
  zr <- sw[FALSE,]
  checkIdentical(rbind(DataFrame(), zr, zr[,1:2]), zr)
  checkIdentical(as.data.frame(rbind(DataFrame(), zr, sw)), swiss)
  target <- do.call(rbind, swisssplit)
  current <- do.call(rbind, lapply(swisssplit, DataFrame))
  rownames(target) <- rownames(current) <- NULL
  checkIdentical(target, as.data.frame(current))
  DF <- DataFrame(A=I(list(1:3)))
  df <- as.data.frame(DF)
  checkIdentical(as.data.frame(rbind(DF, DF)), rbind(df, df))
  
  ## combining factors
  df1 <- data.frame(species = c("Mouse", "Chicken"), n = c(5, 6))
  DF1 <- DataFrame(df1)
  df2 <- data.frame(species = c("Human", "Chimp"), n = c(1, 2))
  DF2 <- DataFrame(df2)
  df12 <- rbind(df1, df2)
  rownames(df12) <- NULL
  checkIdentical(as.data.frame(rbind(DF1, DF2)), df12)
 
  checkIdentical(rownames(rbind(sw, DataFrame(swiss))),
                 c(rownames(swiss), rownames(swiss)))
  checkIdentical(rownames(do.call(rbind, lapply(swisssplit, DataFrame))),
                 unlist(lapply(swisssplit, rownames), use.names=FALSE))

  checkException(rbind(sw[,1:2], sw), silent = TRUE)
  other <- sw
  colnames(other)[1] <- "foo"
  checkException(rbind(other, sw), silent = TRUE)
}

test_DataFrame_looping <- function() {
  data(iris)
  actual <- by(iris, iris$Species, nrow)
  ## a bit tricky because of the 'call' attribute
  attr(actual, "call")[[1]] <- as.name("by")
  iris <- DataFrame(iris, row.names=rownames(iris))
  checkIdentical(actual, by(iris, iris$Species, nrow))
}

test_DataFrame_annotation <- function() {
  df <- DataFrame(x = c(1L, 3L, NA), y = c(10L, 2L, NA))
  mcols(df) <- DataFrame(a = 1:2)
  checkIdentical(mcols(df)[,1], 1:2)
  checkIdentical(mcols(df[2:1])[,1], 2:1)
  checkIdentical(mcols(cbind(df,df))[,1], rep(1:2,2))
  df$z <- 1:3
  checkIdentical(mcols(df, use.names=TRUE),
                 DataFrame(a = c(1L, 2L, NA), row.names = c("x", "y", "z")))
}

## '[<-' setter
test_DataFrame_Setter <- function() {

  .SingleBracket <- function(df0, df1, idx) {
    target <- df0
    for (i in seq_len(length(df0))[idx])
        target[[i]] <- df1[[i]]
    df <- df0
    df[idx] <- df1[idx]
    stopifnot(identical(target, df))
    df <- DataFrame(df0)
    df[idx] <- DataFrame(df1)[idx]
    if (!identical(DataFrame(target), df))
        FALSE
    else
        TRUE
  }

  df0 <- data.frame(x=11:12, y=21:22, z=31:32)
  df1 <- data.frame(matrix(LETTERS[1:6], ncol=3))

  checkTrue(.SingleBracket(df0, df1, c(FALSE, FALSE, TRUE)))
  checkTrue(.SingleBracket(df0, df1, c(TRUE, FALSE, TRUE)))
  checkTrue(.SingleBracket(df0, df1, c(TRUE, TRUE, TRUE)))
  checkTrue(.SingleBracket(df0, df1, TRUE))
 
  target <- df0
  target[] <- df1[]
  df <- DataFrame(df0)
  df[] <- DataFrame(df1)[]
  checkIdentical(DataFrame(target), df)

  for (i in c('a', 'c', 'e')) {
      DF <- DataFrame(A=1:5, row.names=letters[1:5])
      df <- data.frame(A=1:5, row.names=letters[1:5])
      DF[i, 'B'] <- df[i, 'B'] <- 1
      checkIdentical(as.data.frame(DF), df)
  }
}

test_DataFrame_droplevels <- function() {
  df <- DataFrame(state.name, state.region, state.region.rle=Rle(state.region))
  df2 <- head(df)
  checkIdentical(lapply(droplevels(df2), levels),
                 list(state.name=NULL,
                      state.region=c("South", "West"),
                      state.region.rle=c("South", "West")))
}

test_DataFrame_transform <- function() {
  DF <- DataFrame(state.name, state.region, state.area)
  df <- as.data.frame(DF)

  checkIdentical(transform(DF), DF)

  TF <- transform(DF, log.area = log(state.area), ratio = log.area / state.area)
  tf <- transform(transform(df, log.area = log(state.area)),
                  ratio = log.area / state.area)
  checkIdentical(tf, as.data.frame(TF))
}
