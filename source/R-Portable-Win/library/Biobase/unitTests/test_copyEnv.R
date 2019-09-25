basic_test <- function() {
    env1 <- new.env(hash=TRUE, parent=emptyenv())
    env1$a <- 1:3
    env1$b <- "foo"
    env1[[".hidden"]] <- 3.1
    z <- copyEnv(env1)
    checkEquals(env1$a, z$a)
    checkEquals(env1$b, z$b)
    checkTrue(!exists(".hidden", z))
    z <- copyEnv(env1, all.names=TRUE)
    checkEquals(env1[[".hidden"]], z[[".hidden"]])
    env1$a[1] <- 10L
    checkEquals(1L, z$a[1])
}
