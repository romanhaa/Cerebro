## these assume we're in Biobase/inst/unitTests

createCurrentInstances <- function(instanceDir = "VersionedClass_data") {
    require("Biobase")
    classes <- getClasses("package:Biobase")
    isVirtual <- sapply(classes, function(nm) getClass(nm)@virtual)
    nonvirtualClasses <- classes[!isVirtual]

    instances <- sub(".Rda", "",
                     list.files(path=file.path(instanceDir, "devel"), pattern=".*.Rda"))

    need <- nonvirtualClasses[!nonvirtualClasses %in% instances]
    if (length(need)!=0) {
        cat("need:", need, "\n")
        lapply(need, function(cls) {
            cat("creating", cls, "\n")
            assign(cls,  new(cls))
            save(list=cls,
                 file=file.path(instanceDir, "devel", paste(cls, ".Rda", sep="")))
        })
    } else cat("no instances need creating\n")
}

createComponentClasses <- function(ExpressionSet, vers="devel", instanceDir = "VersionedClass_data") {
    MIAME <- experimentData(ExpressionSet)
    AnnotatedDataFrame <- phenoData(ExpressionSet)
    cat("creating MIAME\n")
    save(MIAME, file=file.path(instanceDir, vers, "MIAME.Rda"))
    cat("creating AnnotatedDataFrame\n")
    save(AnnotatedDataFrame, file=file.path(instanceDir, vers, "AnnotatedDataFrame.Rda"))
}
