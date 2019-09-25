basicProps <- list(weight=1, color="blue", friends=c("bob", "alice"))


testCreation <- function() {
    aset <- new("attrData", defaults=basicProps)
    checkEquals(TRUE, is(aset, "attrData"))
}


testDefaultAttributesGetting <- function() {
    aset <- new("attrData", defaults=basicProps)

    ## Get the entire list
    checkEquals(basicProps, attrDefaults(aset))

    ## Get a single attribute
    checkEquals(basicProps$weight, attrDefaults(aset, attr="weight"))
    checkEquals(basicProps$friends, attrDefaults(aset, attr="friends"))

    ## It is an error to ask for an undefined attr
    checkException(attrDefaults(aset, attr="NOSUCHATTRIBUTE"), silent=TRUE)
    ## You can only ask for one attr at a time
    checkException(attrDefaults(aset, attr=c("weight", "friends")), silent=TRUE)
}


testDefaultAttributesSetting <- function() {
    aset <- new("attrData", defaults=basicProps)

    ## Edit default value, type changes are allowed
    attrDefaults(aset, attr="weight") <- 100
    checkEquals(100, attrDefaults(aset, attr="weight"))

    ## Add a new attribute
    attrDefaults(aset, attr="size") <- c(1, 2)
    checkEquals(c(1, 2), attrDefaults(aset, attr="size"))
        
    ## This is sort of dangerous, but for now, we'll allow it.
    ## I would prefer the attributes to be write-once and there
    ## forever.  Or at least a special interface to remove unwanted...
    newProps <- list(dir="home", size=100)
    attrDefaults(aset) <- newProps
    checkEquals(newProps, attrDefaults(aset))
}


testItemGettingAndSettingSimple <- function() {
    aset <- new("attrData", defaults=basicProps)

    ## access to defaults
    checkEquals(1, attrDataItem(aset, x="k1", attr="weight")[[1]])
    expect <- as.list(rep(1, 3))
    names(expect) <- letters[1:3]
    checkEquals(expect, attrDataItem(aset, x=letters[1:3], attr="weight"))

    ## mixed custom/defaults
    attrDataItem(aset, x="k1", attr="weight") <- 900
    checkEquals(900, attrDataItem(aset, x="k1", attr="weight")[[1]])

    ## Retrieve entire attribute list 
    expect <- basicProps
    expect[["weight"]] <- 900
    checkEquals(expect, attrDataItem(aset, x=c("k1", "newone"))[[1]])
    checkEquals(basicProps, attrDataItem(aset, x=c("k1", "newone"))[[2]])

    ## error on unknown attrs
    checkException(attrDataItem(aset, "k1", "UNKNOWN"), silent=TRUE)
    checkException(attrDataItem(aset, "k1", "UNKNOWN") <- "BAD", silent=TRUE)
}


testItemGettingAndSettingVectorized <- function() {
    aset <- new("attrData", defaults=basicProps)
    keys <- c("k1", "k2", "k3")
    
    ## Set multiple with same value 1
    attrDataItem(aset, x=keys, attr="weight") <- 222
    expectPerElem <- basicProps
    expectPerElem[["weight"]] <- 222
    for(k in keys)
      checkEquals(expectPerElem, attrDataItem(aset, k)[[1]])

    ## Set multiple with same value 2
    complexVal <- list(a=as.list(1:3), b="ccc", c=1:5)
    attrDataItem(aset, x=keys, attr="weight") <- list(complexVal)
    expectPerElem <- basicProps
    expectPerElem[["weight"]] <- complexVal
    checkEquals(expectPerElem, attrDataItem(aset, "k1")[[1]])
    checkEquals(expectPerElem, attrDataItem(aset, "k2")[[1]])
    checkEquals(expectPerElem, attrDataItem(aset, "k3")[[1]])

    ## Set multiple with same value 3
    complexVal <- list(a=as.list(1:3), b="ccc", c=1:5, d="extra")
    attrDataItem(aset, x=keys, attr="weight") <- list(complexVal)
    expectPerElem <- basicProps
    expectPerElem[["weight"]] <- complexVal
    for(k in keys)
      checkEquals(expectPerElem, attrDataItem(aset, k)[[1]])

    ## Set multiple with distinct values 1
    wVect <- c(10, 20, 30)
    attrDataItem(aset, x=keys, attr="weight") <- wVect
    for (i in 1:length(wVect))
      checkEquals(wVect[i], attrDataItem(aset, keys[i], "weight")[[1]])

    ## Set multiple with distinct values 2
    wVect <- list(list(a=1), list(a=2), list(a=3))
    attrDataItem(aset, x=keys, attr="weight") <- wVect
    for (i in 1:length(wVect))
      checkEquals(wVect[[i]], attrDataItem(aset, keys[i], "weight")[[1]])
}


testItemRemovalSimple <- function() {
    aset <- new("attrData", defaults=basicProps)

    attrDataItem(aset, x="k1", attr="weight") <- 900
    checkEquals(900, attrDataItem(aset, x="k1", attr="weight")[[1]])

    removeAttrDataItem(aset, x="k1") <- NULL
    checkEquals(1, attrDataItem(aset, x="k1", attr="weight")[[1]])
}


testNames <- function() {
    aset <- new("attrData", defaults=basicProps)

    attrDataItem(aset, x=letters, attr="weight") <- 1:26
    checkEquals(letters, names(aset))

    names(aset) <- LETTERS
    checkEquals(LETTERS, names(aset))

    ## must have right length
    checkException(names(aset) <- letters[1:4], silent=TRUE)
    ## can't have NA
    bad <- letters
    bad[11] <- NA
    checkException(names(aset) <- bad, silent=TRUE)
    ## can't have duplicates
    bad[11] <- "b"
    checkException(names(aset) <- bad, silent=TRUE)
}
