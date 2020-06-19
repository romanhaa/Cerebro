## Use of unsafeSetSlot is currently 'safe' inside a repelacement
## method, because under specific conditions the object is duplicated
## (referenced exactly once) on the way in.  This test tries to catch
## any changes to R that eliminate this duplication, and hence
## invalidates the unsafe slot assignments used in Biobase

testUnsafeSetSlotInSlotAssignment1 <- function() {
    ## calls like set(x) <- value 
    setClass("A", representation=representation(x="numeric"),
             where=.GlobalEnv)
    on.exit(removeClass("A", .GlobalEnv))
    setGeneric("set<-",
               function(object, value) standardGeneric("set<-"),
               where=.GlobalEnv)
    setReplaceMethod("set",
                     signature=signature(object="A"),
                     function(object, value) {
                         Biobase:::unsafeSetSlot(object, "x", value)
                     },
                     where=.GlobalEnv)
    on.exit(removeMethod("set<-", c(object="A"), .GlobalEnv))
    on.exit(removeGeneric("set<-", .GlobalEnv))

    b <- a <- new("A", x=1:10)
    set(a) <- 10:1
    checkIdentical(a@x, 10:1)
    checkIdentical(b@x, 1:10)
}

testUnsafeSetSlotInSlotAssignment2 <- function() {
    ## calls like x <- "set<-"(x, value)
    setClass("A", representation=representation(x="numeric"),
             where=.GlobalEnv)
    on.exit(removeClass("A", .GlobalEnv))
    setGeneric("set<-",
               function(object, value) standardGeneric("set<-"),
               where=.GlobalEnv)
    setReplaceMethod("set",
                     signature=signature(object="A"),
                     function(object, value) {
                         object@x <- value
                         Biobase:::unsafeSetSlot(object, "x", value)
                     },
                     where=.GlobalEnv)
    on.exit(removeMethod("set<-", c(object="A"), .GlobalEnv))
    on.exit(removeGeneric("set<-", .GlobalEnv))

    b <- a <- new("A", x=1:10)
    a <- "set<-"(a, 10:1)
    checkIdentical(a@x, 10:1)
    checkIdentical(b@x, 1:10)
}
