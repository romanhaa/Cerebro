nonvirtualClasses <- NULL
instanceDir <-
    system.file("unitTests", "VersionedClass_data", package="Biobase")

.setUp <- function() {
    nms <- ls(getNamespace("Biobase"),all=TRUE)
    classes <- gsub(".__C__", "", nms[grep(".__C__", nms)])
    classes <- classes[!classes %in% "phenoData"]
    isVirtual <- sapply(classes, function(nm) getClass(nm)@virtual)
    nonvirtualClasses <<- classes[!isVirtual]
}

testNewVersioned <- function() {
  new("Versioned")
  new("VersionedBiobase")
  ## Biobase:::.unversionedObj

  v <- new("Versioned", versions=list(x="1.0.0"))
  checkTrue(classVersion(v)["x"]=="1.0.0")

  ## use class definition defined in RUnit.R
##   a <- new("A")
##   checkTrue(all(classVersion(a) == classVersion(new("VersionedBiobase", versions=list(A="1.0.1")))))
##   checkTrue(all(a@x==10:1))
##   a <- new("A", x=1:10)
##   checkTrue(all(classVersion(a) == classVersion(new("VersionedBiobase", versions=list(A="1.0.1")))))
##   checkTrue(all(a@x==1:10))
##   a <- new("A", x=1:10, versions=list(x="1.0.1"))
##   checkTrue(all(a@x==1:10))
##   checkTrue(all(classVersion(a) == classVersion(new("VersionedBiobase", versions=list(A="1.0.1",x="1.0.1")))))
}

testIsVersioned <- function() {
  checkTrue(isVersioned(new("VersionedBiobase")))
  checkTrue(isVersioned("VersionedBiobase"))
  checkTrue(!isVersioned(1:10))
}

testClassVersion <- function() {
  classVersion(new("Versioned"))
  classVersion(new("VersionedBiobase"))
  checkTrue(is(classVersion(new("VersionedBiobase")),"Versions"))

  checkTrue(all(classVersion("VersionedBiobase") == classVersion(new("VersionedBiobase"))))
}

testClassVersionReplace <- function() {
  ref <- obj <- new("VersionedBiobase")
  classVersion(obj)["x"] <- "1.0.0"
  checkTrue(all(classVersion(obj)[names(classVersion(ref))] == classVersion(ref)))
  checkTrue(classVersion(obj)["x"] == "1.0.0")
                              
  y <- new("Versions", y="1.0.1")
  classVersion(obj)[names(y)] <- y
  checkTrue(all(classVersion(obj)[names(classVersion(ref))] == classVersion(ref)))
  checkTrue(classVersion(obj)["x"] == "1.0.0")
  checkTrue(classVersion(obj)["y"] == "1.0.1")
  checkTrue(classVersion(obj)["y"] != "1.0.0")

  obj <- ref
  classVersion(obj) <- y
  checkTrue(classVersion(obj)["y"] == "1.0.1")
}

testClassVersionSubset <- function() {
    obj <- new("Versions", x="1.0.0")
    checkTrue(obj[1]=="1.0.0")
    checkTrue(obj["x"]=="1.0.0")
    checkException(obj["y"], silent=TRUE)
}

testClassVersionCompare <- function() {
    obj <- new("Versions", x="1.0.0", y="2.0.1")
    checkTrue(all(obj == c(x="1.0.0", y="2.0.1")))
    checkTrue(all(obj == c(y="2.0.1", x="1.0.0")),
              msg="versions in different order")
    checkTrue(!any(obj == c(y="1.0.0", x="2.0.1")),
              msg="incorrectly named elements")
    checkTrue(!all(obj == c(y="1.0.0", x="1.0.0")),
              msg="one element incorrect elements")
    checkException(all(obj == c(x="1.0.0", z="2.0.1")),
                   msg="different version elements",
                   silent=TRUE)
    checkException(obj == c(x="1.0.0"),
                   msg="different version lengths",
                   silent=TRUE)

    ## as above, but comparing version objects
    checkTrue(all(obj == new("Versions", x="1.0.0", y="2.0.1")))
    checkTrue(all(obj == new("Versions", y="2.0.1", x="1.0.0")))
    checkTrue(!any(obj == new("Versions", x="2.0.1", y="1.0.0")))
    checkTrue(!all(obj == new("Versions", x="1.0.0", y="1.0.0")))
    checkException(all(obj == new("Versions", x="1.0.0", z="2.0.1")),
                   msg="different version elements",
                   silent=TRUE)
    checkException(obj == new("Versions", x="1.0.0"),
                   msg="different version lengths",
                   silent=TRUE)
}

testIsCurrent <- function() {
  checkTrue(is.na(isCurrent(1:10)))
  checkTrue(all(isCurrent(new("VersionedBiobase"))))
}

testDevelInstanceArchived <- function() {
    ## archived devel instance?
    instances <- sub(".Rda", "",
                     list.files(path=file.path(instanceDir, "devel"),
                                pattern=".*.Rda"))
    checkTrue(all(nonvirtualClasses %in% instances),
              msg=paste(nonvirtualClasses[!nonvirtualClasses %in% instances],
                collapse=" "))
}

testDevelInstanceIsCurrent <- function() {
    ## overall class is current
    instanceEnv <- new.env(parent=emptyenv())
    lapply(nonvirtualClasses, function(cls) {
        cls <- paste(cls, ".Rda", sep="")
        load(file.path(instanceDir, "devel", cls), instanceEnv)
    })
    instances <-
      sub(".Rda", "",
          list.files(path=file.path(instanceDir, "devel"), pattern=".*.Rda"))
    instances <- instances[!instances %in% "Versioned"]
    current <- sapply(instances, function(obj) {
        vers <- isCurrent(get(obj, env=instanceEnv))[c("S4", obj)]
        all(vers[!is.na(vers)])
    })
    currentv <- current[!is.na(current)]
    checkTrue(all(currentv), msg=paste(names(currentv)[!currentv], collapse=" "))
}
