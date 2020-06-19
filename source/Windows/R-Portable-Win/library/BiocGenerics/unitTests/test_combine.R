###

checkDataFramesEqual <- function(obj1, obj2)
{
    checkTrue(identical(row.names(obj1), row.names(obj2)))
    checkTrue(identical(colnames(obj1), colnames(obj2)))
    checkTrue(all(sapply(colnames(obj1), function(nm) identical(obj1[[nm]], obj2[[nm]]))))
}

test_combine_df <- function()
{
    ## no warnings
    x <- data.frame(x=1:5,y=letters[1:5], row.names=letters[1:5])
    y <- data.frame(z=3:7,y=letters[c(3:5,1:2)], row.names=letters[3:7])
    z <- combine(x,y)
    checkDataFramesEqual(x, z[1:5, colnames(x)])
    checkDataFramesEqual(y, z[3:7, colnames(y)])

    ## an error -- content mismatch
    x <- data.frame(x=1:3, y=letters[1:3], row.names=letters[1:3])
    y <- data.frame(z=2:4, y=letters[1:3], row.names=letters[2:4])
    checkException(suppressWarnings(combine(x,y)), silent=TRUE)

    ## a warning -- level coercion
    oldw <- options("warn")
    options(warn=2)
    on.exit(options(oldw))
    x <- data.frame(x=1:2, y=letters[1:2], row.names=letters[1:2])
    y <- data.frame(z=2:3, y=letters[2:3], row.names=letters[2:3])
    checkException(combine(x,y), silent=TRUE)
    options(oldw)
    checkDataFramesEqual(suppressWarnings(combine(x,y)),
                         data.frame(x=c(1:2, NA),
                                    y=letters[1:3],
                                    z=c(NA, 2:3),
                                    row.names=letters[1:3]))
}

test_combine_df_preserveNumericRows <- function()
{
    dfA <- data.frame(label=rep("x", 2), row.names=1:2)
    dfB <- data.frame(label=rep("x", 3), row.names=3:5)
    dfAB <- combine(dfA, dfB)
    ## preserve integer row names if possible
    checkEquals(1:5, attr(dfAB, "row.names"))

    ## silently coerce row.names to character
    dfC <- data.frame(label=rep("x", 2), row.names=as.character(3:4))
    dfAC <- combine(dfA, dfC)
    checkEquals(as.character(1:4), attr(dfAC, "row.names"))
}

test_combine_df_NoRow <- function()
{
    x <- data.frame(x=1,y=letters[1])[FALSE,]
    y <- data.frame(z=1,y=letters[1])[FALSE,]
    z <- combine(x,x)
    checkTrue(identical(dim(z), as.integer(c(0,2))))
    x <- data.frame(x=1,y=letters[1])[FALSE,]
    y <- data.frame(z=1,y=letters[1])
    z <- combine(x,y)
    checkTrue(identical(dim(z), as.integer(c(1,3))))
    checkTrue(is.na(z$x))
    z <- combine(y,x)
    checkTrue(identical(dim(z), as.integer(c(1,3))))
    checkTrue(is.na(z$x))
}

test_combine_df_OneRow <- function()
{
    x <- data.frame(x=1,y=letters[1], row.names=letters[1])
    y <- data.frame(z=3,y=letters[1], row.names=letters[2])
    z <- combine(x,y)
    checkTrue(identical(dim(z), as.integer(c(2,3))))
    checkTrue(z$x[[1]]==1)
    checkTrue(all(is.na(z$x[[2]]), is.na(z$z[[1]])))
    z <- combine(x,data.frame())
    checkTrue(identical(dim(z), as.integer(c(1,2))))
    checkTrue(all(z[,1:2]==x[,1:2]))
    z <- combine(data.frame(),x)
    checkTrue(identical(dim(z), as.integer(c(1,2))))
    checkTrue(all(z[,1:2]==x[,1:2]))
}

test_combine_df_NoCol <- function()
{
    ## row.names
    obj1 <- data.frame(numeric(20), row.names=letters[1:20])[,FALSE]
    obj <- combine(obj1, obj1)
    checkTrue(identical(obj, obj1))
    ## no row.names -- fails because row.names not recoverable from data.frame?
    obj1 <- data.frame(numeric(20))[,FALSE]
    obj <- combine(obj1, obj1)
    checkTrue(all(dim(obj)==dim(obj1)))
}

test_combine_df_NoCommonCols <- function()
{
    x <- data.frame(x=1:5, row.names=letters[1:5])
    y <- data.frame(y=3:7, row.names=letters[3:7])
    z <- combine(x,y)
    checkTrue(all(dim(z)==as.integer(c(7,2))))
    checkTrue(all(z[1:5,"x"]==x[,"x"]))
    checkTrue(all(z[3:7,"y"]==y[,"y"]))
    checkTrue(all(which(is.na(z))==6:9))
}

test_combine_df_Empty <- function()
{
    z <- combine(data.frame(), data.frame())
    checkTrue(identical(dim(z), as.integer(c(0,0))))
    x <- data.frame(x=1,y=letters[1], row.names=letters[1])
    z <- combine(x,data.frame())
    checkTrue(identical(dim(z), as.integer(c(1,2))))
    checkTrue(identical(z["a",1:2], x["a",1:2]))
    z <- combine(data.frame(), x)
    checkTrue(identical(dim(z), as.integer(c(1,2))))
    checkTrue(identical(z["a",1:2], x["a",1:2]))
}

test_combine_df_AsIs <- function()
{
    x <- data.frame(x=I(1:5),y=I(letters[1:5]), row.names=letters[1:5])
    y <- data.frame(z=I(3:7),y=I(letters[3:7]), row.names=letters[3:7])
    z <- combine(x,y)
    checkTrue(all(sapply(z, class)=="AsIs"))
}

test_combine_df_ColNamesSuffix <- function()
{
    obj1 <- data.frame(a=1:5, a.x=letters[1:5])
    obj2 <- data.frame(a=1:5, a.y=LETTERS[1:5], b=5:1)
    obj <- combine(obj1, obj2)
    checkDataFramesEqual(obj,
                         data.frame(a=1:5, a.x=letters[1:5], a.y=LETTERS[1:5], b=5:1))
}

test_combine_3df <- function()
{
    ## data.frame's are tricky, because c(df, list(...)) unlists df
    x <- data.frame(x=1:5,
                    y=factor(letters[1:5], levels=letters[1:8]),
                    row.names=letters[1:5])
    y <- data.frame(z=3:7,
                    y=factor(letters[3:7], levels=letters[1:8]),
                    row.names=letters[3:7])
    w <- data.frame(w=4:8,
                    y=factor(letters[4:8], levels=letters[1:8]),
                    row.names=letters[4:8])
    res <- combine(w, x, y)

    e <- data.frame(w=c(4:8, rep(NA, 3)),
                    y=c(letters[c(4:8, 1:3)]),
                    x=c(4:5, rep(NA, 3), 1:3),
                    z=as.integer(c(4:7, rep(NA, 3), 3)),
                    row.names=letters[c(4:8, 1:3)])
    checkIdentical(e, res)
}

test_combine_df_POSIXct <- function()
{
    ## class(x) can have length > 1 as in Sys.time()
    t0 <- Sys.time()
    df1 <- data.frame(i = 1:3, t = rep(t0, 3), row.names=letters[1:3])
    df2 <- data.frame(i = 1:3, t = c(t0, t0 + 500, t0 + 1000),
                      row.names=c("a", "d", "e"))
    e <- data.frame(i = c(1L, 2L, 3L, 2L, 3L),
                    t = c(t0, t0, t0, t0 + 500, t0 + 1000),
                    row.names=c("a", "b", "c", "d", "e"))
    res <- combine(df1, df2)
    checkIdentical(e, res)
}


test_combine_df_WithNamedArgs <- function() {
    x <- data.frame(x=1:5,
                    y=factor(letters[1:5], levels=letters[1:8]),
                    row.names=letters[1:5])
    y <- data.frame(z=3:7,
                    y=factor(letters[3:7], levels=letters[1:8]),
                    row.names=letters[3:7])
    w <- data.frame(w=4:8,
                    y=factor(letters[4:8], levels=letters[1:8]),
                    row.names=letters[4:8])
    checkIdentical(combine(w, y, x), combine(w, x, y=y))
    checkIdentical(combine(w, y, x), combine(w, y=y, x))
    checkIdentical(combine(x, y, w), combine(w, y=y, x=x))
    checkIdentical(combine(x, y, w), combine(y=y, x=x, w))
}

test_combine_mat <- function()
{
    ## dimnames
    m <- matrix(1:20, nrow=5, dimnames=list(LETTERS[1:5], letters[1:4]))
    checkEquals(m, combine(m, m))
    checkEquals(m, combine(m[1:3,], m[4:5,]))
    checkEquals(m, combine(m[,1:3], m[,4, drop=FALSE]))
    ## overlap
    checkEquals(m, combine(m[1:3,], m[3:5,]))
    checkEquals(m, combine(m[,1:3], m[,3:4]))
    checkEquals(matrix(c(1:3, NA, NA, 6:8, NA, NA, 11:15, NA, NA, 18:20),
                       nrow=5,
                       dimnames=list(LETTERS[1:5], letters[1:4])),
                combine(m[1:3,1:3], m[3:5, 3:4]))
    ## row reordering
    checkEquals(m[c(1,3,5,2,4),], combine(m[c(1,3,5),], m[c(2,4),]))
    ## Exceptions
    checkException(combine(m, matrix(0, nrow=5, ncol=4)),
                   silent=TRUE)         # types differ
    checkException(combine(m, matrix(0L, nrow=5, ncol=4)),
                   silent=TRUE)         # attributes differ
    m1 <- matrix(1:20, nrow=5)
    checkException(combine(m, m1), silent=TRUE) # dimnames required
}

test_combine_mat_DifferentModes <- function()
{
    m <- matrix(1:20, nrow=5, dimnames=list(LETTERS[1:5], letters[1:4]))
    n <- matrix(as.numeric(1:20),
                nrow=5, dimnames=list(LETTERS[1:5], letters[1:4]))
    res <- combine(m, n)                # modes coerced to same
    checkEquals("numeric", mode(res))
    n <- matrix(as.character(1:20),
                nrow=5, dimnames=list(LETTERS[1:5], letters[1:4]))
    checkException(combine(m, n))       # modes differ
}

