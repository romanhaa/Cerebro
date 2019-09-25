test_getOrgPkg_load_only <- function() {
    ## check that map between chip and org package works with loaded
    ## but not attached org package
    if ("package:org.Hs.eg.db" %in% search()) {
        detach("package:org.Hs.eg.db")
        on.exit(attachNamespace("org.Hs.eg.db"))
    }
    pkg <- "hgu95av2.db"
    env <- loadNamespace(pkg)
    db <- get(pkg, env)
    keys <- head(AnnotationDbi::keys(db))
    df <- AnnotationDbi::select(db, keys, "SYMBOL")
    checkIdentical(c(6L, 2L), dim(df))
}
