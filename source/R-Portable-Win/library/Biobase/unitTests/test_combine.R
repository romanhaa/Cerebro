testCombineOne <- function() {
    data(sample.ExpressionSet)
    checkIdentical(sample.ExpressionSet,
                   combine(sample.ExpressionSet))
}

testCombineTwo <- function() {
    data(sample.ExpressionSet)
    x <- sample.ExpressionSet
    checkTrue(all.equal(x, combine(x[,1:5],x[,seq(6, ncol(x))])))
}

testCombineThree <- function() {
    data(sample.ExpressionSet)
    x <- sample.ExpressionSet
    y <- combine(x[,1:5],x[, 6:10], x[,seq(11, ncol(x))])
    checkTrue(all.equal(x, y))
}

