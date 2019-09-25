fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
obo <- getOBOCollection(fl)

test_obo_OBOCollection <- function() {
    obj <- OBOCollection()
    checkTrue(validObject(obj))
    checkEquals(c(0,3), dim(GSEABase:::.kv(obj)))
    checkEquals(c(0,1), dim(GSEABase:::.stanza(obj)))
    checkEquals(c(0,1), dim(GSEABase:::.obo_subset(obj)))

    obj <- OBOCollection(letters[1:5])
    checkTrue(validObject(obj))
    checkEquals(c(5,3), dim(GSEABase:::.kv(obj)))
    checkEquals(c(5,1), dim(GSEABase:::.stanza(obj)))
    checkEquals(c(0,1), dim(GSEABase:::.obo_subset(obj)))
    ## value should be NA_character_
    checkTrue(all(is.na(GSEABase:::.kv(obj)$value)))
    
    obj <- OBOCollection(c("GO:0008967", "GO:0015119"))
    checkTrue(validObject(obj))
    checkEquals(c(2,3), dim(GSEABase:::.kv(obj)))
    checkEquals(c(2,1), dim(GSEABase:::.stanza(obj)))
    checkEquals(c(0,1), dim(GSEABase:::.obo_subset(obj)))
    ## terms should be filled with definitions
    checkTrue(!any(is.na(GSEABase:::.kv(obj)$value)))
}


test_getOBOCollection <- function() {
    obj <- getOBOCollection(fl, evidenceCode="TAS")
    checkEquals(106, length(ids(obj)))
    checkEquals("TAS", evidenceCode(obj))
}

test_obo_subsets <- function() {
    s <- c("Generic GO slim", "GOA and proteome slim",
           "Plant GO slim", "Yeast GO slim", "Prokaryotic GO subset")
    names(s) <- c("goslim_generic", "goslim_goa",
                  "goslim_plant", "goslim_yeast", "gosubset_prok")

    checkIdentical(s, subsets(obo))
    checkIdentical(s, subsets(obo, "named"))
    checkIdentical(names(s), subsets(obo, "key"))
    checkIdentical(as.vector(s), subsets(obo, "value"))
    checkIdentical(paste(names(s), " (", s, ")", sep=""),
                   subsets(obo, "full"))
}

test_obo_subset <- function() {
    checkIdentical(obo, obo[])

    obo0 <- obo[evidenceCode="TAS"]
    checkIdentical("TAS", evidenceCode(obo0))
    checkIdentical(ids(obo), ids(obo0))

    obo1 <- obo["goslim_yeast"]
    checkIdentical("goslim_yeast", subsets(obo1, "key"))
    checkEquals(46, length(ids(obo1)))

    obo2 <- obo[c("goslim_yeast", "gosubset_prok")]
    checkIdentical(c("goslim_yeast", "gosubset_prok"),
                   subsets(obo2, "key"))
    checkEquals(90, length(ids(obo2)))

    obo3 <- obo[c("goslim_yeast", "gosubset_prok"),
                evidenceCode="TAS"]
    checkIdentical(ids(obo2), ids(obo3))
    checkIdentical("TAS", evidenceCode(obo3))

    checkException(obo[ids="GO:0000003"],
                   msg="ids must match those implied by subsets",
                   silent=TRUE)
    checkTrue(validObject(obo[ids=ids(obo)]))
}

test_obo_as_graphNEL_empty_obj <- function() {
    validObject(as(OBOCollection(), "graphNEL"))
    validObject(as(graphNEL(), "OBOCollection"))
}

test_obo_as_graphNEL <- function() {
    obo0 <- obo["goslim_goa"]
    ## some loss of information here
    g <- as(obo0, "graphNEL")
    ## round trip should be identical
    x <- as(g, "OBOCollection")
    checkEquals(41, length(ids(x)))
    g1 <- as(x, "graphNEL")
    checkEquals(41, length(nodes(g1)))
    checkEquals(30, length(unlist(edges(g1))))
    checkIdentical(g, g1)
    checkIdentical(x, as(g1, "OBOCollection"))
}

test_goSlim <- function() {
    id <- GOCollection("GO:0005618")
    res <- goSlim(id, obo, "CC")
    checkIdentical(26L, nrow(res))
    checkIdentical(4L, sum(res$Count))

    id <- GOCollection(rep("GO:0005618", 2))
    res <- goSlim(id, obo, "CC")
    checkIdentical(26L, nrow(res))
    checkIdentical(8L, sum(res$Count))
}
