test_unwritableDirectories <- function()
{
    if (.Platform$OS.type != "unix")
        return()
    .unwritableDirectories <- BiocInstaller:::.unwritableDirectories

    dir.create(f <- tempfile(), mode="400")
    dir.create(g <- tempfile())
    on.exit(Sys.chmod(f, mode="777"))
    res <- suppressWarnings({
        .unwritableDirectories(c(f, g))
    })
    checkIdentical(f, res)
}
