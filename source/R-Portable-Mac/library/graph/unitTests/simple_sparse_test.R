.indexToCoord <- graph:::.indexToCoord
.coordToIndex <- graph:::.coordToIndex
set.seed(0x389d)

check_coord_txfm <- function(n)
{
    n <- as.integer(n)
    seqn <- seq_len(n)
    seqnn <- seq_len(n * n)
    coord <- cbind(rep(seqn, n), rep(seqn, each = n))
    index <- .coordToIndex(coord[ , 1L], coord[, 2L], n)
    checkEquals(seqnn, index)

    got_coord <- .indexToCoord(seqnn, n)
    checkEquals(coord, got_coord)

    ## now check order preservation
    p1 <- sample(seqnn)
    p1coord <- coord[p1, ]
    if (n == 1L) {
        dim(p1coord) <- c(1L, 2L)
    }
    checkEquals(p1, .coordToIndex(p1coord[,1L], p1coord[,2L], n))
    checkEquals(p1coord, .indexToCoord(p1, n))
}

test_coordinate_txfm <- function()
{
    for (i in 1:25) {
        check_coord_txfm(i)
    }
}


