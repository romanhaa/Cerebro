###

test_updateObject_list <- function()
{
    setClass("A",
             representation(x="numeric"), prototype(x=1:10),
             where=.GlobalEnv)
    a <- new("A")
    l <- list(a,a)
    checkTrue(identical(l, updateObject(l)))

    setMethod("updateObject", "A",
              function(object, ..., verbose=FALSE) {
                  if (verbose) message("updateObject object = 'A'")
                  object@x <- -object@x
                  object
              },
              where=.GlobalEnv)

    obj <- updateObject(l)
    checkTrue(identical(lapply(l, function(elt) { elt@x <- -elt@x; elt }),
                        obj))
    removeMethod("updateObject", "A", where=.GlobalEnv)
    removeClass("A", where=.GlobalEnv)
}

test_updateObject_env <- function()
{
    opts <- options()
    options(warn=-1)
    e <- new.env()
    e$x=1
    e$.x=1
    obj <- updateObject(e)
    checkTrue(identical(e,obj))         # modifies environment

    lockEnvironment(e)
    obj <- updateObject(e)              # copies environment
    checkTrue(identical(lapply(ls(e, all=TRUE), function(x) x),
                        lapply(ls(obj, all=TRUE), function(x) x)))
    checkTrue(!identical(e, obj))       # different environments

    e <- new.env()
    e$x=1
    e$.x=1
    lockBinding("x", e)
    checkException(updateObject(e), silent=TRUE)

    lockEnvironment(e)
    obj <- updateObject(e)
    checkTrue(TRUE==bindingIsLocked("x", obj)) # R bug, 14 May, 2006, fixed
    checkTrue(FALSE==bindingIsLocked(".x", obj))
    options(opts)
}

test_updateObject_defaults <- function()
{
    x <- 1:10
    checkTrue(identical(x, updateObject(x)))
}

test_updateObject_S4 <- function()
{
    setClass("A",
             representation=representation(
               x="numeric"),
             prototype=list(x=1:5),
             where=.GlobalEnv)
    .__a__ <- new("A")
    setClass("A",
             representation=representation(
               x="numeric",
               y="character"),
             where=.GlobalEnv)
    checkException(validObject(.__a__), silent=TRUE)      # now out-of-date
    .__a__@x <- 1:5
    a <- updateObject(.__a__)
    checkTrue(validObject(a))
    checkIdentical(1:5, a@x)
    removeClass("A", where=.GlobalEnv)
}

test_updateObject_setClass <- function()
{
    setClass("A",
             representation(x="numeric"),
             prototype=prototype(x=1:10),
             where=.GlobalEnv)
    a <- new("A")
    checkTrue(identical(a,updateObject(a)))
    removeClass("A", where=.GlobalEnv)
}

test_updateObject_refClass <- function()
{
    cls <- ".__test_updateObject_refClassA"
    .A <- setRefClass(cls, fields=list(x="numeric", y="numeric"),
                      where=.GlobalEnv)

    a <- .A()
    checkTrue(all.equal(a, updateObject(a)))

    a <- .A(x=1:5, y=5:1)
    checkTrue(all.equal(a, updateObject(a)))

    .A <- setRefClass(cls, fields=list(x="numeric", y="numeric", z="numeric"),
                      where=.GlobalEnv)
    checkTrue(all.equal(.A(x=1:5, y=5:1, z=numeric()), updateObject(a)))

    .A <- setRefClass(cls, fields=list(x="numeric"))
    warn <- FALSE
    value <- withCallingHandlers(updateObject(a), warning=function(w) {
        txt <- "dropping fields(s) 'y' from object = '.__test_updateObject_refClassA'"
        warn <<- identical(txt, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    checkTrue(warn)
    checkTrue(all.equal(.A(x=1:5), value))
    
    removeClass(cls, where=.GlobalEnv)
}
