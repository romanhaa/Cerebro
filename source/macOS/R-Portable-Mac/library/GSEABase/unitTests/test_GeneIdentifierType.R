constructors <- local({
    idTypes <- names(slot(getClass("GeneIdentifierType"), "subclasses"))
    idTypes[!idTypes %in% c("AnnotationIdentifier", "KEGGFrameIdentifier")]
})

test_GeneIdentifierType_Constructors <- function() {
    ## do they exist and return the correct class?
    frameIds <- paste0("GOAll", "FrameIdentifier")
    for (i in seq_along(constructors)) {
        res <- do.call(constructors[[i]], list())
        checkTrue(validObject(res))
        checkTrue(is(res, constructors[[i]]))
        if (!constructors[[i]] %in% frameIds) {
            res <- do.call(constructors[[i]], list("foo"))
            checkTrue(validObject(res))
            checkTrue(is(res, constructors[[i]]))
        }
    }

    ## Required slot for AnnotationIdentifier
    res <- AnnotationIdentifier(annotation="123")
    checkTrue(validObject(res))
    checkTrue(is(res, "AnnotationIdentifier"))
    checkTrue("hgu95av2"==annotation(AnnotationIdentifier("hgu95av2")))
}

test_GeneIdentifierType_geneIdType <- function() {
    data(sample.ExpressionSet)

    gs <- GeneSet(sample.ExpressionSet[100:110],
                  setName="123", setIdentifier="456")
    suppressWarnings(geneIdType(gs) <- "EntrezIdentifier")
    checkTrue(validObject(gs, complete=TRUE))
    checkTrue(is(geneIdType(gs), "EntrezIdentifier"))

    gs <- GeneSet(sample.ExpressionSet[100:110],
                  setName="123", setIdentifier="456")
    suppressWarnings(geneIdType(gs) <- EntrezIdentifier())
    checkTrue(validObject(gs, complete=TRUE))
    checkTrue(is(geneIdType(gs), "EntrezIdentifier"))

    ## duplicate gene names exception
    gs <- GeneSet(sample.ExpressionSet[100:200],
                  setName="123", setIdentifier="456")
    opt <- options(warn=2)
    on.exit(options(opt))
    checkException(geneIdType(gs, verbose=TRUE) <- EntrezIdentifier(),
                   silent=TRUE)
}

test_GeneIdentifierType_mapIdentifiers_toAnnotation <- function() {
    map <- getAnnMap("SYMBOL", "hgu95av2")

    gss <- getBroadSets(system.file("extdata", "Broad.xml",
                                    package="GSEABase"))
    suppressWarnings({
        res <- mapIdentifiers(gss[[1]], AnnotationIdentifier("hgu95av2"))
    })
    checkTrue(validObject(res))
    gids <- mget(geneIds(gss[[1]]), revmap(map), ifnotfound=NA)
    gids <- gids[!is.na(gids)]
    checkEquals(sort(unique(unlist(gids))), sort(geneIds(res)))
    checkIdentical(res,
                   mapIdentifiers(res, AnnotationIdentifier()))
    checkIdentical(res,
                   mapIdentifiers(res, AnnotationIdentifier("hgu95av2")))
    checkException(mapIdentifiers(res, AnnotationIdentifier("dfs")),
                   silent=TRUE)
}

test_GeneIdentifierType_mapIdentifiers_toAnnotation_via_Dbi <- function()  {
    map <- getAnnMap("SYMBOL", "hgu95av2")
    gss <- getBroadSets(system.file("extdata", "Broad.xml",
                                    package="GSEABase"))
    suppressMessages(suppressWarnings({
        res <- mapIdentifiers(gss[[1]], AnnotationIdentifier("hgu95av2.db"))
    }))
    checkTrue(validObject(res))
    gids <- mget(geneIds(gss[[1]]), revmap(map), ifnotfound=NA)
    gids <- gids[!is.na(gids)]
    checkEquals(sort(unique(unlist(gids))), sort(geneIds(res)))
}

test_GeneIdentifierType_mapIdentifiers_from_to_Annotation <- function() {
    gs <- GeneSet("EGF")
    geneIdType(gs) <- AnnotationIdentifier("org.Hs.eg.db")
    gs1 <- gs
    geneIdType(gs1) <- AnnotationIdentifier("org.Hs.eg.db")
    checkIdentical(gs, gs1)
}

test_GeneIdentifierType_mapIdentifiers_from_Annotation <- function() {
    data(sample.ExpressionSet)
    gs <- GeneSet(sample.ExpressionSet[100:110])
    ## AnnotationIdentifier --> AnnotationIdentifier
    checkIdentical(gs,
                   mapIdentifiers(gs, AnnotationIdentifier()))
    checkIdentical(gs,
                   mapIdentifiers(gs, AnnotationIdentifier("hgu95av2")))
    checkException(mapIdentifiers(gs, AnnotationIdentifier("fdsf")),
                   silent=TRUE)
    res <- mapIdentifiers(gs, EntrezIdentifier())
    checkTrue(validObject(res))
    checkIdentical(geneIdType(res), EntrezIdentifier("hgu95av2"))

    ## AnnotationIdentifier --> Other Annotation-based Identifier
    res <- mapIdentifiers(gs, EntrezIdentifier("hgu95av2"))
    checkTrue(validObject(res))
    checkIdentical(geneIdType(res), EntrezIdentifier("hgu95av2"))
    checkException(mapIdentifiers(gs, EntrezIdentifier("sdfs")),
                   silent=TRUE)
}

test_GeneIdentifierType_mapIdentifiers_AnnDbBimap <- function() {
    library(org.Hs.eg.db)
    gs <- mapIdentifiers(GeneSet("4214"), SymbolIdentifier(),
                         org.Hs.egSYMBOL)
    checkIdentical("MAP3K1", geneIds(gs))
    checkIdentical("Symbol", geneIdType(geneIdType(gs)))
}

test_GeneIdentifierType_mapIdentifiers_verbose_warnings <-
    function()
{
    ## duplicate gene names exception
    gs <- GeneSet(sample.ExpressionSet[100:200],
                  setName="123", setIdentifier="456")
    opt <- options(warn=2)
    on.exit(options(opt))
    checkException(mapIdentifiers(gs,  EntrezIdentifier(), verbose=TRUE),
                   silent=TRUE)
}


test_GeneIdentifierType_show <- function() {
    f <- function(x) capture.output(show(x))
    checkEquals("geneIdType: Annotation", f(AnnotationIdentifier()))
    checkEquals("geneIdType: Annotation (hgu95av2)",
                f(AnnotationIdentifier("hgu95av2")))
    checkEquals("geneIdType: EntrezId", f(EntrezIdentifier()))
    checkEquals("geneIdType: EntrezId (hgu95av2)",
                f(EntrezIdentifier("hgu95av2")))
    checkEquals("geneIdType: Genename", f(GenenameIdentifier()))
    checkEquals("geneIdType: Genename (hgu95av2)",
                f(GenenameIdentifier("hgu95av2")))
}

test_GeneIdentifierType_mapIdentifiers_isNullMap <- function() {
    f <- GSEABase:::.mapIdentifiers_isNullMap
    ai <- AnnotationIdentifier()
    aia <- AnnotationIdentifier("hgu95av2")
    checkTrue(!f(ai, ai, FALSE))        # not mappable at all
    checkTrue(!f(ai, aia, FALSE))       # 
    checkTrue(f(aia, aia, FALSE))
}

test_GeneIdentifierType_mapIdentifiers_isMappable <- function() {
    f <- GSEABase:::.mapIdentifiers_isMappable
    ai <- AnnotationIdentifier()
    aia <- AnnotationIdentifier("hgu95av2")
    si <- SymbolIdentifier()

    ## both with annotation -- ok only when annotation same
    checkTrue(f(aia, aia))
    checkException(f(aia, AnnotationIdentifier("hgu133plus2")),
                   silent=TRUE)

    ## one with annotation -- ok
    checkTrue(f(ai, aia))
    checkTrue(f(aia, si))
    checkTrue(f(si, aia))

    ## neither with annotation -- err
    checkException(f(ai, ai), silent=TRUE)
}

test_GeneIdentifierType_mapIdentifiers_normalize <- function() {
    f <- GSEABase:::.mapIdentifiers_normalize
    ai <- AnnotationIdentifier()
    aia <- AnnotationIdentifier("hgu95av2")
    si <- SymbolIdentifier()
    sia <- SymbolIdentifier("hgu95av2")

    checkIdentical(list(aia, aia), f(ai, aia))
    checkIdentical(list(aia, aia), f(aia, ai))
    checkIdentical(list(aia, aia), f(aia, aia))

    checkIdentical(list(sia, aia), f(si, aia))
    checkIdentical(list(sia, aia), f(sia, ai))
    checkIdentical(list(sia, aia), f(sia, aia))

    checkException(f(ai, ai), silent=TRUE)
}

test_GeneIdentifierType_mapIdentifiers_selectMaps <- function() {
    f <- GSEABase:::.mapIdentifiers_selectMaps
    ## all fully specified and mappable at this point
    ## single maps
    chk1 <- function(map, id1, id2) {
        res <- f(id1, id2)
        checkTrue(length(res)==1)
        all.equal(map, res[[1]])
    }
    ## two maps
    chk2 <- function(map1, map2, id1, id2) {
        res <- f(id1, id2)
        checkTrue(length(res)==2)
        all.equal(map1, res[[1]])
        all.equal(map2, res[[2]])
    }

    pkg <- "hgu95av2"
    ## need to look at EntrezId too, as these could also be dispatched
    ## on for chip-based packages
    ai <- AnnotationIdentifier(pkg)
    ei <- EntrezIdentifier(pkg)
    si <- SymbolIdentifier(pkg)
    gi <- GenenameIdentifier(pkg)
    chk1(hgu95av2ENTREZID, ai, ei)
    chk1(hgu95av2SYMBOL, ai, si)
    chk1(revmap(hgu95av2ENTREZID), ei, ai)
    chk1(revmap(hgu95av2SYMBOL), si, ai)
    chk2(revmap(hgu95av2ENTREZID), hgu95av2SYMBOL, ei, si)
    chk2(revmap(hgu95av2SYMBOL), hgu95av2ENTREZID, si, ei)
    chk2(revmap(hgu95av2GENENAME), hgu95av2SYMBOL, gi, si)
    chk2(revmap(hgu95av2SYMBOL), hgu95av2GENENAME, si, gi)

    pkg <- "org.Hs.eg.db"
    ei <- EntrezIdentifier(pkg)
    si <- SymbolIdentifier(pkg)
    gi <- GenenameIdentifier(pkg)

    chk1(org.Hs.egSYMBOL, ei, si)
    chk1(revmap(org.Hs.egSYMBOL), si, ei)
    chk2(revmap(org.Hs.egSYMBOL), org.Hs.egGENENAME, si, gi)
}

test_GeneIdentifierType_mapIdentifiers_map <- function() {
    f <- GSEABase:::.mapIdentifiers_map
    ## ids 300:310 of sample.ExpressionSet; not 1:1 maps below; These
    ## are hand-validated

    aids <- c("31539_r_at", "31540_at", "31541_at", "31542_at",
              "31543_at", "31544_at", "31545_at", "31546_at",
              "31547_at", "31548_at", "31549_at")

    eids <- c("3604", "5275", "2312", "58503", "2299", "6222", "6141",
              "141", "4142")

    sids <- c("TNFRSF9", "SERPINB13", "FLG", "OPRPN", "FOXI1",
              "RPS18", "RPL18", "ADPRH", "MAS1")

    gids <- c("TNF receptor superfamily member 9",
              "serpin family B member 13",
              "filaggrin",  "opiorphin prepropeptide",
              "forkhead box I1", "ribosomal protein S18",
              "ribosomal protein L18",
              "ADP-ribosylarginine hydrolase",
              "MAS1 proto-oncogene, G protein-coupled receptor")

    pkg <- "hgu95av2"
    ai <- AnnotationIdentifier(pkg)
    ei <- EntrezIdentifier(pkg)
    si <- SymbolIdentifier(pkg)

    checkEquals(eids, f(aids, ai, ei))
    checkEquals(sids, f(aids, ai, si))
    ## the following is convenient, not a general property
    checkEquals(sids, f(f(aids, ai, ei), ei, si))

    pkg <- "org.Hs.eg.db"
    sids <- c("TNFRSF9", "SERPINB13", "FLG", "OPRPN", "FOXI1",
              "RPS18", "RPL18", "ADPRH", "MAS1")
    ei <- EntrezIdentifier(pkg)
    si <- SymbolIdentifier(pkg)
    gi <- GenenameIdentifier(pkg)

    checkEquals(sids, f(eids, ei, si))
    checkEquals(eids, f(sids, si, ei))
    checkEquals(gids, f(sids, si, gi))
}

test_GeneIdentifierType_mapIdentifiers_revMap <- function() {
    f <- GSEABase:::.mapIdentifiers_revMap
    pkg <- "hgu95av2"
    si <- SymbolIdentifier(pkg)
    ai <- AnnotationIdentifier(pkg)

    ## no annotation for IMAGINARY
    sids <- c("IMAGINARY", "XPO1", "LBR")
    checkEquals(list("288_s_at" = "LBR", "37729_at" = "XPO1"),
                f(sids, si, ai))
    ## Multiple annotations for MAP2
    sids <- c("MAP2", "XPO1", "LBR")
    checkEquals(list("183_at" = "MAP2", "1972_s_at" = "MAP2",
                     "219_i_at" = "MAP2", "220_r_at" = "MAP2",
                     "288_s_at" = "LBR", "35422_at" = "MAP2",
                     "37729_at" = "XPO1"),
                f(sids, si, ai))
    ## "35422_at"  "220_r_at" map to same sym
    aids <- c("35422_at", "220_r_at", "41229_at")
    checkEquals(list(MAP2 = c("35422_at", "220_r_at"), NFIB = "41229_at"),
                f(aids, ai, si))

    ## two-step maps
    gi <- GenenameIdentifier(pkg)
    sids <- c("IMAGINARY", "XPO1", "LBR")
    checkEquals(list("exportin 1" = "XPO1",
                     "lamin B receptor" = "LBR"),
                f(sids, si, gi))
    sids <- c("MAP2", "XPO1", "LBR")
    checkEquals(list("exportin 1" = "XPO1",
                     "lamin B receptor" = "LBR",
                     "microtubule associated protein 2" = "MAP2"),
                f(sids, si, gi))
}

## edge cases

test_GeneIdentifierType_mapIdentifiers_nullAmbiguity <- function() {
    ## Original bug: 
    ##     1: Ambiguous method selection for "mapIdentifiers", target "GeneSet#AnnotationIdentifier#NullIdentifier" (the first of the signatures shown will be used)
    ##     GeneSet#AnnotationIdentifier#GeneIdentifierType
    ##     GeneSet#GeneIdentifierType#NullIdentifier
    opts <- options(warn=2)
    on.exit(options(opts))
    gs <- GeneSet(setName="123", setIdentifier="345")
    geneIdType(gs) <- AnnotationIdentifier("xyz")
    checkTrue(validObject(gs))
}

test_GeneIdentifierType_AnnoOrEntrez <- function() {
    checkIdentical(AnnotationIdentifier("hgu95av2.db"),
                   AnnoOrEntrezIdentifier("hgu95av2.db"))
    checkIdentical(AnnotationIdentifier("hgu95av2"),
                   AnnoOrEntrezIdentifier("hgu95av2"))
    checkIdentical(EntrezIdentifier("org.Hs.eg.db"),
                   AnnoOrEntrezIdentifier("org.Hs.eg.db"))
}
