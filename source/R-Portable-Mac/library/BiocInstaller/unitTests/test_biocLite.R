test_reposType <- function()
{
    .rRepos <- BiocInstaller:::.rRepos
    .githubRepos <- BiocInstaller:::.githubRepos

    r <- "foo"
    http <- c("http://foo.bar/baz", "https://foo.bar/baz")
    github <- c("foo/bar", "foo/bar@baz")
    all <- c(r, http, github)

    checkIdentical(c(r, http), .rRepos(all))
    checkIdentical(github, .githubRepos(all))
}
