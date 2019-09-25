testNew <- function() {
    ## create instances of non-virtual objects with simple call to "new"
    nms <- ls(getNamespace("Biobase"),all=TRUE)
    classes <- gsub(".__C__", "", nms[grep(".__C__", nms)])
    isVirtual <- sapply(classes, function(nm) getClass(nm)@virtual)
    res <- lapply(classes[!isVirtual],
                  function(x) suppressWarnings(new(x)))
}

test_MIAME_construction <- function() {
    checkTrue(validObject(new("MIAME")))
    checkTrue(validObject(MIAME()))
    checkIdentical(new("MIAME", name="mytest"), MIAME(name="mytest"))
}
