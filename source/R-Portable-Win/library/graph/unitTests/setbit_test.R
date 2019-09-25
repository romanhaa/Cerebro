.testbit <- graph:::testbit
setbitv <- graph:::setbitv
setbit <- graph:::setbit
bitToMat <- graph:::bitToMat
makebits <- graph:::makebits

test_setbitv <- function()
{
    len <- 5L * 8L
    xx <- makebits(len)
    for (i in seq_len(len)) {
        checkEquals(FALSE, .testbit(xx, i))
        xx <- setbitv(xx, i, FALSE)
        checkEquals(FALSE, .testbit(xx, i))
    }

    xx <- makebits(len)
    for (i in seq_len(len)) {
        checkEquals(FALSE, .testbit(xx, i))
        xx <- setbitv(xx, i, TRUE)
        checkEquals(TRUE, .testbit(xx, i), msg=paste("iter", i))
    }
    for (i in seq_len(len)) checkEquals(TRUE, .testbit(xx, i))
    for (i in seq_len(len)) xx <- setbitv(xx, i, FALSE)
    for (i in seq_len(len)) checkEquals(FALSE, .testbit(xx, i))
}

test_setbitv_vectorized <- function()
{
    len <- 3L * 8L
    xx <- makebits(len)
    xxall <- setbitv(xx, 1:10, rep(TRUE, 10))
    checkTrue(all(.testbit(xxall, 1:10)))
    checkTrue(!any(.testbit(xxall, 11:len)))
}

test_setbit_basics <- function()
{
    xx <- makebits(40L)
    tf <- function(n) rawToBits(setbit(xx, n)[1])
    got <- as.logical(do.call(rbind, lapply(1:8, tf)))
    dim(got) <- c(8, 8)
    want <- matrix(FALSE, nrow=8, ncol=8)
    diag(want) <- TRUE
    checkEquals(want, got)

    tf <- function(n) rawToBits(setbit(xx, n)[2])
    got <- as.logical(do.call(rbind, lapply(9:16, tf)))
    dim(got) <- c(8, 8)
    want <- matrix(FALSE, nrow=8, ncol=8)
    diag(want) <- TRUE
    checkEquals(want, got)

    x2 <- setbit(xx, 38)
    for (i in 1:4) checkTrue(!as.logical(x2[i]))
    checkTrue(as.logical(x2[5]))
}

test_testbit <- function()
{
    xx <- makebits(40)
    for (i in 1:(5 * 8)) {
        checkTrue(!.testbit(xx, i))
    }
    checkTrue(!any(.testbit(xx, 1:40)))

    xx <- as.raw(rep(255L, 5))
    for (i in 1:(5 * 8)) {
        checkTrue(.testbit(xx, i), i)
    }
    checkTrue(all(.testbit(xx, 1:40)))

    xx <- setbit(as.raw(5), 23)
    checkTrue(.testbit(xx, 23))
    checkEquals(c(TRUE, TRUE), .testbit(xx, c(23, 23)))
    checkEquals(c(FALSE, TRUE), .testbit(xx, c(21, 23)))
    checkEquals(c(TRUE, FALSE), .testbit(xx, c(23, 24)))
}

rand_bitarray_matrix <- function(nrow, nset)
{
    bv <- makebits(nrow^2, bitdim=c(nrow, nrow))
    idx <- sample(1:(nrow^2), nset)
    setbitv(bv, idx, rep(1L, length(idx)))
}

test_bitToMat <- function()
{
    make_mat <- function(nrow)
    {
        matrix(sample(c(0L, 1L), nrow^2, replace = TRUE), nrow = nrow)
    }

    do_test <- function(nrow)
    {
        m <- make_mat(nrow)
        bv <- makebits(nrow^2, bitdim=c(nrow, nrow))
        bv <- setbitv(bv, which(m == 1L), rep(1L, sum(m)))
        checkEquals(m, bitToMat(bv))
    }

    sizes <- c(1, 2, 5, 13, 26)
    reps <- c(2, 3, 25, 25, 25)
    for (i in seq_along(sizes)) {
        size <- sizes[i]
        for (j in seq_len(reps[i]))
            do_test(size)
    }
}

test_bitarray_transpose <- function()
{
    nreps <- 25L
    sizes <- as.integer(c(1, 2, 4, 5, 6, 7, 23, 24))
    for (size in sizes) {
        for (i in seq_len(nreps)) {
            v0 <- rand_bitarray_matrix(size, ceiling((size^2) %/% 2))
            m0 <- bitToMat(v0)
            want <- t(m0)
            vt <- .Call(graph:::graph_bitarray_transpose, v0)
            mt <- bitToMat(vt)
            checkEquals(want, mt)
        }
    }
}
