.colorBroadSets <- function() {
    gss <- getBroadSets(system.file("extdata", "Broad.xml",
                                    package="GSEABase"))
    lapply(gss, function(gs) {
        gcs <- GeneColorSet(gs, phenotype="<undefined>")
        df <- coloring(gcs)
        df[,"geneColor"] <- as.factor(rep(c("P", "M"), length=nrow(df)))
        df[,"phenotypeColor"] <- as.factor(rep(LETTERS[1:3], length=nrow(df)))
        coloring(gcs) <- df
        gcs
    })
}

.check_subset_Ok <- function(expected, gs) {
    checkTrue(validObject(gs))
    len <- length(geneIds(expected))
    checkIdentical(len, length(geneIds(gs)))
    checkIdentical(len, length(geneColor(gs)))
    checkIdentical(len, length(phenotypeColor(gs)))
    checkTrue(all(geneIds(expected)==geneIds(gs)))
    checkTrue(all(geneColor(expected)==geneColor(gs)))
    checkTrue(all(phenotypeColor(expected)==phenotypeColor(gs)))
}

test_GCS_ConstructorNoColorSetArgs <- function() {
    checkException(GeneColorSet(setIdentifier="123",
                                setName="Set name"),
                   silent=TRUE)
}

test_GCS_ConstructorDefaultArgs <- function() {
    gs <- GeneColorSet(setIdentifier="123",
                       setName="Set name",
                       phenotype="A phenotype")
    checkTrue(validObject(gs, complete=TRUE))
    checkTrue(length(geneIds(gs))==0)
    checkTrue(phenotype(gs)=="A phenotype")
    checkTrue(length(geneColor(gs))==0)
    checkTrue(length(phenotypeColor(gs))==0)
}

test_GCS_ConstructorAllColorArgs <- function() {
    ## appropriate default colors
    gs <- GeneColorSet(setIdentifier="123",
                       setName="Set name",
                       phenotype="A phenotype",
                       geneIds=LETTERS[1:24])
    checkTrue(validObject(gs, complete=TRUE))
    checkIdentical(geneIds(gs), LETTERS[1:24])
    checkTrue(phenotype(gs)=="A phenotype")
    checkTrue(length(geneColor(gs))==24)
    checkTrue(length(phenotypeColor(gs))==24)

    ## correct color contents
    gfactor <- factor(rep(c("high", "low"), 12))
    pfactor <- factor(rep(c("big", "medium", "small"), 8))
    gs <- GeneColorSet(setIdentifier="123",
                       setName="Set name",
                       phenotype="A phenotype",
                       geneIds=LETTERS[1:24],
                       geneColor=gfactor,
                       phenotypeColor=pfactor)
    checkTrue(validObject(gs, complete=TRUE))
    checkTrue(phenotype(gs)=="A phenotype")
    checkIdentical(geneColor(gs), gfactor)
    checkIdentical(phenotypeColor(gs), pfactor)
}

test_GCS_show <- function() {
    gs <- GeneColorSet(setIdentifier="123",
                       setName="Set name",
                       phenotype="A phenotype",
                       geneIds=LETTERS[1:24])
    con <- textConnection("tmp", open="w", local=TRUE)
    sink(con)
    on.exit(sink())
    show(gs)
}

test_GCS_colorizeReplace <- function() {
    gcs <- .colorBroadSets()[[1]]
    df <- coloring(gcs)
    gc <- as.factor(rep(c("Up", "Down"), length=nrow(df)))
    pc <- as.factor(rep(c("Bigger", "Smaller", "Same"),
                        length=nrow(df)))
    df$geneColor <- gc
    df$phenotypeColor <- pc
    coloring(gcs) <- df

    checkIdentical(geneColor(gcs), gc)
    checkIdentical(phenotypeColor(gcs), pc)
}

test_GCS_colorizeReplaceRetainGeneOrder <- function() {
    gcs <- .colorBroadSets()[[1]]
    ogeneIds <- geneIds(gcs)
    coloring(gcs) <- coloring(gcs)[sample(ogeneIds, length(ogeneIds)),]
    checkIdentical(geneIds(gcs), ogeneIds)
}

test_GCS_intersect <- function() {
    gcss <- .colorBroadSets()

    res <- GSEABase::intersect(gcss[[1]], gcss[[2]])
    checkTrue(validObject(res, complete=TRUE))
    checkTrue(length(geneIds(res))==0)
    checkTrue(length(urls(res))==4)
    checkIdentical(levels(geneColor(gcss[[1]])),
                   levels(geneColor(res)))
    checkIdentical(levels(phenotypeColor(gcss[[1]])),
                   levels(phenotypeColor(res)))

    res <- GSEABase::intersect(gcss[[1]], gcss[[1]])
    checkTrue(validObject(res, complete=TRUE))
    checkIdentical(geneIds(gcss[[1]]), geneIds(res))
    checkIdentical(urls(gcss[[1]]), urls(res))
    checkIdentical(levels(geneColor(gcss[[1]])),
                   levels(geneColor(res)))
    checkIdentical(levels(phenotypeColor(gcss[[1]])),
                   levels(phenotypeColor(res)))
}

test_GCS_intersectDifferentColors <- function() {
    gcs1 <- .colorBroadSets()[[1]]
    gcs2 <- gcs1
    phenotype(gcs2) <- paste(phenotype(gcs2), "A")
    geneColor(gcs2) <- factor(rep(c("Q", "R"), length=length(geneIds(gcs2))))
    ## warning about synthetic phenotype
    oldOpts <- options(warn=2)
    on.exit(options(oldOpts))
    checkException(res <- GSEABase::intersect(gcs1, gcs2), silent=TRUE)
    options(oldOpts)
    ## 
    suppressWarnings(res <- GSEABase::intersect(gcs1, gcs2))
    checkTrue(phenotype(res) != phenotype(gcs1))
    checkIdentical(2L, length(levels(geneColor(res))))
    checkIdentical(3L, length(levels(phenotypeColor(res))))
    checkTrue(!any(levels(geneColor(res)) == levels(geneColor(gcs2))))
}

test_GCS_union <- function() {
    gcss <- .colorBroadSets()

    res <- union(gcss[[1]], gcss[[2]])
    checkTrue(validObject(res, complete=TRUE))
    checkIdentical(length(geneIds(res)),
                   sum(sapply(gcss,
                              function(x) length(geneIds(x)))))
    checkTrue(all(geneIds(res) %in% c(geneIds(gcss[[1]]),
                                    geneIds(gcss[[2]]))))
    checkIdentical(levels(geneColor(gcss[[1]])),
                   levels(geneColor(res)))
    checkIdentical(levels(phenotypeColor(gcss[[1]])),
                   levels(phenotypeColor(res)))
    checkTrue(all(urls(res) %in% unlist(sapply(gcss, urls))))

    res <- union(gcss[[1]], gcss[[1]])
    checkTrue(validObject(res))
    checkIdentical(geneIds(res), geneIds(gcss[[1]]))
    checkIdentical(geneColor(res), geneColor(gcss[[1]]))
    checkIdentical(phenotypeColor(res), phenotypeColor(gcss[[1]]))
    checkIdentical(urls(gcss[[1]]), urls(res))
}

test_GCS_setdiff <- function() {
    gcss <- .colorBroadSets()

    res <- GSEABase::setdiff(gcss[[1]], gcss[[2]])
    checkTrue(validObject(res, complete=TRUE))
    checkIdentical(geneIds(gcss[[1]]), geneIds(res))
    checkIdentical(geneColor(gcss[[1]]), geneColor(res))
    checkIdentical(phenotypeColor(gcss[[1]]), phenotypeColor(res))
    
    res <- GSEABase::setdiff(gcss[[2]], gcss[[1]])
    checkTrue(validObject(res, complete=TRUE))
    checkIdentical(geneIds(gcss[[2]]), geneIds(res))
    checkIdentical(geneColor(gcss[[2]]), geneColor(res))
    checkIdentical(phenotypeColor(gcss[[2]]), phenotypeColor(res))

    res <- GSEABase::setdiff(gcss[[1]], gcss[[1]])
    checkTrue(validObject(res, complete=TRUE))
    checkTrue(length(geneIds(res))==0)
    checkIdentical(levels(geneColor(gcss[[1]])), levels(geneColor(res)))
    checkIdentical(levels(phenotypeColor(gcss[[1]])), levels(phenotypeColor(res)))
}

test_GCS_LogicalNonOverlapping <- function() {
    gcss <- .colorBroadSets()
    gs1 <- gcss[[1]]
    gs2 <- gcss[[2]]

    gs12 <- gs1 & gs2
    checkTrue(length(geneIds(gs12))==0)

    gs12 <- gs1 | gs2
    checkTrue(length(geneIds(gs12))==length(c(geneIds(gs1), geneIds(gs2))))
    checkTrue(all(geneIds(gs1) %in% geneIds(gs12)))
    checkTrue(all(geneIds(gs2) %in% geneIds(gs12)))
    checkIdentical(geneIds(GSEABase::setdiff(gs12, gs1)), geneIds(gs2))
    checkIdentical(geneIds(GSEABase::setdiff(gs12, gs2)), geneIds(gs1))
}

test_GCS_LogicalOverlapping <- function() {
    gcss <- .colorBroadSets()
    gs1 <- gcss[[1]]

    idx1 <- sample(seq_along(geneIds(gs1)), 20)
    idx2 <- sample(seq_along(geneIds(gcss[[2]])), 20)

    gs2 <-
        GeneColorSet(type=geneIdType(gs1),
                     setName="123", setIdentifier="456",
                     geneIds=c(geneIds(gs1)[idx1], geneIds(gcss[[2]])[idx2]),
                     phenotype=phenotype(gs1),
                     geneColor=factor(c(
                       as.character(geneColor(gs1))[idx1],
                       as.character(geneColor(gcss[[2]]))[idx2])),
                     phenotypeColor=factor(c(
                       as.character(phenotypeColor(gs1))[idx1],
                       as.character(phenotypeColor(gcss[[2]]))[idx2])))
    checkTrue(all(geneIds(gs1 | gs2) %in%
                  geneIds(GSEABase::setdiff(gs1, gs2) |
                        (gs1 & gs2) |
                        GSEABase::setdiff(gs2, gs1))))
}

test_GCS_subset <- function() {
    gcs <- .colorBroadSets()[[1]]

    expected <- new(class(gcs), gcs,
                    geneIds=geneIds(gcs)[4:1],
                    geneColor=factor(as.character(
                      geneColor(gcs)[4:1])),
                    phenotypeColor=factor(as.character(
                      phenotypeColor(gcs)[4:1])))

    .check_subset_Ok(expected, gcs[4:1])
    .check_subset_Ok(expected, gcs[ geneIds(gcs)[4:1] ])

    max <- length(geneIds(gcs))
    checkException(gcs[(max-1):(max+1)], silent=TRUE)
    checkException(gcs["adfas"], silent=TRUE)
}

test_GCS_subset2 <- function() {
    gcs <- .colorBroadSets()[[1]]

    expected <- c(geneId=geneIds(gcs)[[2]],
                  geneColor=as.character(geneColor(gcs)[[2]]),
                  phenotypeColor=as.character(phenotypeColor(gcs)[[2]]))

    checkIdentical(expected, gcs[[2]])
    checkIdentical(expected, gcs[[ geneIds(gcs)[[2]] ]])

    max <- length(geneIds(gcs))
    checkException(gcs[["sdfsf"]], silent=TRUE)
    checkException(gcs[[ max+1 ]], silent=TRUE)

    checkIdentical(expected, do.call("$",list(gcs, geneIds(gcs)[[2]])))
}

test_GCS_mapIdentifiers <- function() {
    f <- GSEABase:::.mapIdentifiers_map
    ai <- AnnotationIdentifier("hgu95av2")
    si <- SymbolIdentifier("hgu95av2")
    gi <- GenenameIdentifier("hgu95av2")
    ## 1:1 annotations
    gcs <- GeneColorSet(geneIds=c("XPO1", "LBR"),
                        geneIdType=si,
                        phenotype="pheno data",
                        geneColor = factor(c("increase","increase")),
                        phenotypeColor = factor(c("B","A")))
    gcs1 <- mapIdentifiers(gcs, ai)
    checkTrue(validObject(gcs1))
    checkEquals(f(geneIds(gcs), si, ai), geneIds(gcs1))
    checkEquals(geneColor(gcs), geneColor(gcs1))
    checkEquals(phenotypeColor(gcs), phenotypeColor(gcs1))
    ## 0 annotations for IMAGINARY
    gcs <- GeneColorSet(geneIds=c("IMAGINARY", "XPO1", "LBR"),
                        geneIdType=si,
                        phenotype="pheno data",
                        geneColor = factor(c("increase","increase","increase"),
                          levels=c("increase", "decrease")),
                        phenotypeColor = factor(c("A","B","A")))
    gcs1 <- mapIdentifiers(gcs, ai)
    checkTrue(validObject(gcs1))
    checkEquals(f(geneIds(gcs), si, ai), geneIds(gcs1))
    checkEquals(geneColor(gcs)[2:3], geneColor(gcs1))
    checkEquals(phenotypeColor(gcs)[2:3], phenotypeColor(gcs1))
    ## 5 annotations for MAP2
    gcs <- initialize(gcs, geneIds=c("MAP2", "XPO1", "LBR"))
    gcs1 <- mapIdentifiers(gcs, ai)
    checkTrue(validObject(gcs1))
    checkEquals(f(geneIds(gcs), si, ai), geneIds(gcs1))
    checkEquals(phenotypeColor(gcs)[c(rep(1,5), 2:3)], phenotypeColor(gcs1))
    checkEquals(geneColor(gcs)[c(rep(1,5), 2:3)], geneColor(gcs1))
    ## 2-step map; ends up being 1:1
    gcs1 <- mapIdentifiers(gcs, gi)
    checkTrue(validObject(gcs1))
    checkEquals(f(geneIds(gcs), si, gi), geneIds(gcs1))
    checkEquals(phenotypeColor(gcs), phenotypeColor(gcs1))
    checkEquals(geneColor(gcs), geneColor(gcs1))
}

test_GCS_mapIdentifiers_exceptions <- function() {
    ai <- AnnotationIdentifier("hgu95av2")
    si <- SymbolIdentifier("hgu95av2")
    gcs <- GeneColorSet(geneIds=c("183_at", "1972_s_at", "219_i_at",
                          "220_r_at", "35422_at", "37729_at",
                          "288_s_at"),
                        geneIdType=ai,
                        phenotype="phenotype",
                        geneColor=factor(c("decrease", "increase",
                          "increase", "increase", "increase",
                          "increase", "increase")),
                        phenotypeColor=factor(rep("A", 7)))
    checkException(mapIdentifiers(gcs, si), silent=TRUE)
}

test_GCS_mapIdentifiers_NullIdentifier <- function() {
    si <- SymbolIdentifier("hgu95av2")
    ni <- NullIdentifier()
    gcs <- GeneColorSet(geneIds=c("XPO1", "LBR"),
                        geneIdType=si,
                        phenotype="pheno data",
                        geneColor = factor(c("increase","increase")),
                        phenotypeColor = factor(c("B","A")))
    gcs1 <- mapIdentifiers(gcs, ni)
    checkTrue(validObject(gcs1))
    checkEquals(geneIds(gcs), geneIds(gcs1))
    checkEquals(geneColor(gcs), geneColor(gcs1))
    checkEquals(phenotypeColor(gcs), phenotypeColor(gcs))
    
    gcs <- mapIdentifiers(gcs1, si)
    checkTrue(validObject(gcs1))
    checkEquals(geneIds(gcs), geneIds(gcs1))
    checkEquals(geneColor(gcs), geneColor(gcs1))
    checkEquals(phenotypeColor(gcs), phenotypeColor(gcs))
}
