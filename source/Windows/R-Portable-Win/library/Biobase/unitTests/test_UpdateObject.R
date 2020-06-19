testUpdateObjectToDefaults <- function() {
    x <- 1:10
    checkTrue(identical(1:10, updateObjectTo(x, 10:1)))
    x <- as.numeric(1:10)
    checkTrue(identical(as.integer(1:10), updateObjectTo(x, integer())))
    checkTrue(!identical(as.numeric(1:10), updateObjectTo(x, integer())))
}

testUpdateObjectToSetClass <- function() {
    setClass("A",
             representation(x="numeric"),
             prototype=prototype(x=1:10),
             where=.GlobalEnv)
    a <- new("A")
    a1 <- new("A",x=10:1)
    checkTrue(identical(a, updateObjectTo(a, a1)))

    setClass("B",
             representation(x="numeric"),
             where=.GlobalEnv)
    b <- new("B")
    checkException(updateObjectTo(a, b), silent=TRUE)

    setAs("A", "B", function(from) {
        b <- new("B")
        b@x <- from@x
        b
    }, where=.GlobalEnv)
    obj <- updateObjectTo(a,b)
    checkTrue(class(obj)=="B")
    checkIdentical(obj@x, a@x)
    removeMethod("coerce", c("A","B"), where=.GlobalEnv)
    removeClass("B", where=.GlobalEnv)
    removeClass("A", where=.GlobalEnv)
}

testUpdateExpressionSet <- function() {
    opts <- options()
    options(warn=-1)
    obj <- new("ExpressionSet")
    checkTrue(all.equal(obj, updateObject(obj)))
    checkTrue(!identical(new("ExpressionSet"), updateObject(obj))) # different environments
    obj <- new("ExpressionSet", storage.mode="list")
    checkTrue(identical(obj, updateObject(obj)))
    checkTrue(identical(new("ExpressionSet", storage.mode="list"), updateObject(obj))) # same class -- list

    data(sample.ExpressionSet)
    classVersion(sample.ExpressionSet)["eSet"] <- "1.0.0"
    checkException(validObject(sample.ExpressionSet), silent=TRUE)

    suppressMessages(obj <- updateObject(sample.ExpressionSet))
    checkTrue(isVersioned(obj))
    checkTrue(all(isCurrent(obj)))
    checkTrue(validObject(obj))
    checkTrue(identical(lapply(ls(assayData(obj), all=TRUE), function(x) x),
                        lapply(ls(assayData(sample.ExpressionSet),all=TRUE), function(x) x)))
    checkTrue(identical(annotation(obj), annotation(sample.ExpressionSet)))

    suppressMessages(obj1a <- updateObjectTo(sample.ExpressionSet, new("ExpressionSet")))
    ## next better written as(sample.ExpressionSet, "MultiSet")
    suppressMessages(obj1b <- updateObjectTo(sample.ExpressionSet, new("MultiSet")))
    obj2 <- updateObject(obj)           # stop after eSet
    options(opts)
}

testUpdateESetMisc <- function() {
    opts <- options()
    options(warn=-1)

    idx <- c("phenoData", "experimentData", "featureData")
    fun <- function(nm)
        isS4(eval(parse(text=paste(nm,"(obj)", sep=""))))

    load(system.file("unitTests", "VersionedClass_data", "devel",
                     "sample.exprSet.rda", package="Biobase"))
    suppressMessages(obj <- as(sample.exprSet, "ExpressionSet"))
    checkTrue(validObject(obj, complete=TRUE))
    checkTrue(all(sapply(idx, fun)))

    load(system.file("unitTests", "VersionedClass_data", "devel",
                     "sample.eSet.rda", package="Biobase"))
    obj <- as(sample.eSet, "MultiSet")
    checkTrue(validObject(obj, complete=TRUE))
    checkTrue(all(sapply(idx, fun)))

    load(system.file("unitTests", "VersionedClass_data", "devel", "eset.rda",
                     package="Biobase"))
    obj <- as(eset, "ExpressionSet")
    checkTrue(validObject(obj, complete=TRUE))
    checkTrue(all(sapply(idx, fun)))

    options(opts)
}

testUpdateMiscPreviousInstances <- function() {
    opts <- options("warn")
    options(warn=-1)
    on.exit(options(opts))

    rda <- dir(system.file("unitTests", "VersionedClass_data",
                           package="Biobase"), full.names=TRUE,
               recursive=TRUE, pattern="^([^(ExpressionSet)]).*\\.Rda")

    ok <- sapply(rda, function(nm) {
        env <- new.env(parent=emptyenv())
        load(nm, env)
        tryCatch({
            eapply(env, function(elt) {
                suppressMessages(obj <- updateObject(elt))
                checkTrue(isS4(obj))
                checkTrue(validObject(obj, complete=TRUE))
            })
            TRUE
        }, error=function(...) FALSE)
    })
    checkTrue(all(ok),
              msg=sprintf("failed: '%s'", paste(rda[!ok], collapse="' '")))
}

testUpdatePreviousExpressionSet <- function() {
    opts <- options("warn")
    options(warn=-1)
    on.exit(options(opts))

    rda <- dir(system.file("unitTests", "VersionedClass_data",
                           package="Biobase"), full.names=TRUE,
               recursive=TRUE, pattern="^ExpressionSet.*\\.Rda")

    ok <- sapply(rda, function(nm) {
        env <- new.env(parent=emptyenv())
        load(nm, env)
        tryCatch({
            eapply(env, function(elt) {
                suppressMessages(obj <- updateObject(elt))
                checkTrue(validObject(obj, complete=TRUE))
                ## S4
                idx <- c("phenoData", "experimentData", "featureData")
                ok <- sapply(idx, function(nm) {
                    isS4(eval(parse(text=paste(nm,"(obj)", sep=""))))
                })

                checkTrue(all(ok))
                ## content
                checkIdentical(exprs(obj),
                               slot(elt, "assayData")[["exprs"]])
                checkIdentical(pData(phenoData(obj)),
                               slot(slot(elt, "phenoData"), "data"))
                checkIdentical(varMetadata(phenoData(obj)),
                               slot(slot(elt, "phenoData"), "varMetadata"))
                nms <- names(getSlots("MIAME"))
                nms <- nms[!nms %in% ".__classVersion__"]
                lapply(nms, function(nm)
                       checkIdentical(slot(experimentData(obj), nm),
                                      slot(slot(elt, "experimentData"),
                                           nm)))
            })
            TRUE
        }, error=function(...) FALSE)
    })

    checkTrue(all(ok),
              msg=sprintf("failed: '%s'", paste(rda[!ok], collapse="' '")))
}
