testNew <- function() {
    ## default constructor
    obj <- new("NChannelSet")
    checkTrue(length(channelNames(obj))==0)
}

testNew2 <- function() {
    ## assayData from arguments not matching slots
    obj <- new("NChannelSet",
               R=matrix(0,10,5), G=matrix(0,10,5))
    checkTrue(all(c("R", "G") %in% channelNames(obj)))
    checkTrue(all(dim(obj)==c(10,5)))
}

testNew3 <- function() {
    ## explicit assayData, infered phenoData
    obj <- new("NChannelSet",
               assayData=assayDataNew(
                 R=matrix(0,10,5),
                 G=matrix(1,10,5)))
    checkTrue(all(c("R", "G") %in% channelNames(obj)))
    checkTrue(all(c(10,5)==dim(obj)))

    obj <- new("NChannelSet",
               R=matrix(0,10,5, dimnames=list(NULL, LETTERS[1:5])),
               G=matrix(0,10,5, dimnames=list(NULL, letters[1:5])))
    checkTrue(all(c("R", "G") %in% channelNames(obj)))
    checkTrue(all(c(10,5)==dim(obj)))

    obj <- new("NChannelSet",
               R=matrix(0,10,5, dimnames=list(NULL, LETTERS[1:5])),
               G=matrix(0,10,5))
    checkTrue(all(c("R", "G") %in% channelNames(obj)))
    checkTrue(all(c(10,5)==dim(obj)))

    obj <- new("NChannelSet",
               R=matrix(0,10,5),
               G=matrix(0,10,5, dimnames=list(NULL, LETTERS[1:5])))
    checkTrue(all(c("R", "G") %in% channelNames(obj)))
    checkTrue(all(c(10,5)==dim(obj)))
}

testNew4 <- function() {
    ## no assayData, zero-row phenoData
    data <- data.frame(x=numeric(),
                       y=numeric(),
                       z=numeric())
    varMetadata <- data.frame(labelDescription=character(ncol(data)),
                              channel=factor(rep("_ALL_", ncol(data))))
    phenoData <- new("AnnotatedDataFrame",
                     data=data, varMetadata=varMetadata)
    checkTrue(validObject(new("NChannelSet", phenoData=phenoData)))
}

testNew5 <- function() {
    ## explicit assayData, phenoData
    R <- matrix(0, 10, 3, dimnames=list(NULL, LETTERS[1:3]))
    G <- matrix(1, 10, 3, dimnames=list(NULL, LETTERS[1:3]))
    assayData <- assayDataNew(R=R, G=G)

    data <- data.frame(x=numeric(ncol(R)),
                       y=numeric(ncol(R)),
                       z=numeric(ncol(R)))
    varMetadata <- data.frame(labelDescription=character(ncol(data)),
                              channel=factor(rep("_ALL_", ncol(data)),
                                levels=c("R", "G", "_ALL_")))
    phenoData <- new("AnnotatedDataFrame",
                     data=data, varMetadata=varMetadata)

    obj <- new("NChannelSet", assayData=assayData, phenoData=phenoData)
    checkTrue(nrow(varMetadata)==ncol(data))
    checkTrue(all(varMetadata(obj)[["channel"]] == "_ALL_"))
}

testNew5a <- function() {
    ## explicit assayData, phenoData; no names
    R <- matrix(0, 10, 5)
    G <- matrix(1, 10, 5)
    assayData <- assayDataNew(R=R, G=G)

    data <- data.frame(x=numeric(ncol(R)))
    varMetadata <- data.frame(labelDescription=character(ncol(data)),
                              channel=factor(
                                rep("__ALL__", ncol(data)),
                                levels=c("R", "G", "_ALL_")))
    phenoData <- new("AnnotatedDataFrame",
                     data=data, varMetadata=varMetadata)
    obj <- new("NChannelSet",
               assayData=assayData, phenoData=phenoData)
}
                             
testNew6 <- function() {
    ## silently add 'channel' to varMetadata
    R <- matrix(0, 10, 5)
    G <- matrix(1, 10, 5)
    assayData <- assayDataNew(R=R, G=G)
    data <- data.frame(x=numeric(ncol(R)))
    phenoData <- new("AnnotatedDataFrame", data=data)
    obj <- new("NChannelSet",
               assayData=assayData, phenoData=phenoData)
    checkTrue(validObject(obj))

    varMetadata <- data.frame(labelDescription=character(ncol(data)))
    phenoData <- new("AnnotatedDataFrame", data=data)
    obj <- new("NChannelSet",
               assayData=assayData, phenoData=phenoData)
    checkTrue(validObject(obj))
}

testAssayDataGets <- function() {
    assayData <- assayDataNew(R = matrix(0, 10, 5), G = matrix(1,10,5))
    obj <- NChannelSet(assayData = assayData)
    exp <- -1 * assayData(obj)[["G"]]
    assayDataElement(obj, "G") <- exp
    checkIdentical(exp, assayData(obj)[["G"]])

    exp <- assayData(obj)[["R"]]
    assayDataElement(obj, "G") <- NULL
    checkIdentical("R", channelNames(obj))
    checkIdentical("R", assayDataElementNames(obj))
    checkIdentical(exp, assayData(obj)[["R"]])
}

testDifferentSampleNames <- function() {
    ## channels have different identifiers
    assayData <- assayDataNew(R = matrix(0,10,5,
                                dimnames=list(NULL, LETTERS[1:5])),
                              G = matrix(1,10,5,
                                dimnames=list(NULL, letters[1:5])))
    obj <- new("NChannelSet", assayData = assayData)
    checkTrue(sampleNames(obj)[["R"]] == LETTERS[1:5] &&
              sampleNames(obj)[["G"]] == letters[1:5])
}

testSampleNamesUpdate <- function() {
    assayData <- assayDataNew(R = matrix(0,10,5),
                              G = matrix(1,10,5))
    obj <- new("NChannelSet", assayData = assayData)

    sampleNames(obj) <- list(R=LETTERS[1:5])
    checkTrue(validObject(obj))
    checkTrue(all(sampleNames(obj)[["R"]] == LETTERS[1:5]))

    sampleNames(obj) <- list(R=LETTERS[5:1], G=letters[5:1])
    checkTrue(validObject(obj))
    checkTrue(all(sampleNames(obj)[["R"]] == LETTERS[5:1]) &&
              all(sampleNames(obj)[["G"]] == letters[5:1]))

    ## drop to character when smaple names identic
    sampleNames(obj) <- LETTERS[1:5]
    checkIdentical(LETTERS[1:5], sampleNames(obj))
    checkIdentical(list(G=LETTERS[1:5], R=LETTERS[1:5]),
                   eapply(assayData(obj), colnames))

    checkException(sampleNames(obj) <- list(LETTERS[1:5]),
                   silent=TRUE)         # unnamed repleacement
    checkException(sampleNames(obj) <- list(X=LETTERS[1:5]) ,
                   silent=TRUE)         # misnamed repleacement
}

testSelectChannels <- function() {
    obj <- new("NChannelSet",
               assayData=assayDataNew(
                 R=matrix(0,10,5),
                 G=matrix(1,10,5)))
    objR <- selectChannels(obj, "R")
    checkTrue(channelNames(objR) == "R")
    checkIdentical(dim(obj), dim(objR))

    objRG <- selectChannels(obj, c("R", "G"))
    checkTrue(all(c("R", "G") %in% channelNames(objRG)))
    checkIdentical(dim(obj), dim(objRG))

    checkException(selectChannels(obj, c("G", "G")), silent=TRUE)
}

testSelectChannelsCommonPhenoData <- function() {
    obj <- new("NChannelSet",
               assayData=assayDataNew(
                 R=matrix(0,10,5),
                 G=matrix(0,10,5)),
               phenoData=new("AnnotatedDataFrame",
                 data=data.frame(
                   r=letters[1:5],
                   g=LETTERS[1:5],
                   both=1:5),
                 varMetadata=data.frame(
                   labelDescription=c(
                     "r data", "g data", "both data"),
                   channel=factor(
                     c("R", "G", "_ALL_")))))

    objR <- selectChannels(obj, "R")
    checkIdentical("R", channelNames(objR))
    checkTrue(all(c("r", "both") %in% names(pData(objR)) &
                  !"g" %in% names(pData(objR))))
    checkIdentical(c("R", "_ALL_"), levels(varMetadata(objR)[["channel"]]))

    objG <- selectChannels(obj, "G")
    checkIdentical("G", channelNames(objG))
    checkTrue(all(c("g", "both") %in% names(pData(objG)) &
                  !"r" %in% names(pData(objG))))
    checkIdentical(c("G", "_ALL_"), levels(varMetadata(objG)[["channel"]]))

    objRG <- selectChannels(obj, c("R", "G"))
    checkTrue(all(c("R", "G") %in% channelNames(objRG)))
    checkTrue(all(c("r", "g", "both") %in% names(pData(objRG))))
    checkTrue(all(c("R", "G") %in% levels(varMetadata(objRG)[["channel"]])))
}

testChannel <- function() {
    obj <- new("NChannelSet",
               assayData=assayDataNew(
                 R=matrix(1,10,5),
                 G=matrix(2,10,5)),
               phenoData=new("AnnotatedDataFrame",
                 data=data.frame(
                   r=letters[1:5],
                   g=LETTERS[1:5],
                   both=1:5),
                 varMetadata=data.frame(
                   labelDescription=c(
                     "r data", "g data", "both data"),
                   channel=factor(
                     c("R", "G", "_ALL_")))))
    robj <- channel(obj, "R")
    checkTrue(validObject(robj))
    checkTrue(all(varLabels(robj) %in% c("r", "both")))
    checkIdentical(assayDataElementNames(robj), "exprs")
    gobj <- channel(obj, "G")
    checkTrue(validObject(gobj))
    checkTrue(all(varLabels(gobj) %in% c("g", "both")))
    checkIdentical(assayDataElementNames(gobj), "exprs")
}

testChannelNames_replace <- function()
{
    obj <- NChannelSet(R=matrix(-1, 5, 5), G=matrix(1, 5, 5))
    checkIdentical(c("G", "R"), channelNames(obj))
    channelNames(obj) <- c(Gn="G", Rd="R")   ## rename
    checkIdentical(c("Gn", "Rd"), channelNames(obj))
    channelNames(obj) <- c("Rd", "Gn")       ## reorder
    checkIdentical(c("Rd", "Gn"), channelNames(obj))
    checkTrue(all(assayData(obj)[["Gn"]] == 1))

    exp <- "'value' elements must include all channelNames()"
    obs <- tryCatch(channelNames(obj) <- "X", error=conditionMessage)
    checkIdentical(exp, obs)
    obs <- tryCatch(channelNames(obj) <- c(X="Gn"), error=conditionMessage)
    checkIdentical(exp, obs)

    exp <- "duplicated channelNames are not allowed"
    obs <- tryCatch(channelNames(obj) <- c(X="Gn", X="Rd"),
                    error=conditionMessage)
    checkIdentical(exp, obs)
   
}
