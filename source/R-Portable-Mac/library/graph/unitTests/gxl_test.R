simpleWithAttributes <- system.file("GXL/attributesExample.gxl", package="graph")

if (nchar(simpleWithAttributes) == 0) stop("bad gxl path")

testGxlNodes <- function() {
    con <- file(simpleWithAttributes)
    tryCatch({
        g <- fromGXL(con)
        eNodes <- c("p", "q", "v", "w")
        checkEquals(eNodes, nodes(g))
    }, finally=close(con))
}


testGxlEdges <- function() {
    con <- file(simpleWithAttributes)
    tryCatch({
        g <- fromGXL(con)
        eEdges <- list(p=c("v", "q"),
                       q="w",
                       v=character(0),
                       w=character(0))
        checkEquals(eEdges, edges(g))
    }, finally=close(con))
}


testGxlNodeAttrs <- function() {
    con <- file(simpleWithAttributes)
    tryCatch({
        g <- fromGXL(con)
        checkEquals(316, nodeData(g, "w", "line")[[1]])
        checkEquals(225, nodeData(g, "v", "line")[[1]])
        
        checkEquals("main.c", nodeData(g, "p", "file")[[1]])
        checkEquals(555, nodeData(g, "p", "code")[[1]])
        checkEquals(1.234, nodeData(g, "p", "rate")[[1]])
        checkEquals(TRUE, nodeData(g, "p", "pass")[[1]])
        checkEquals(FALSE, nodeData(g, "p", "fail")[[1]])
    checkTrue(is.na(nodeData(g, "p", "line")[[1]]))
    }, finally=close(con))
}


testNodeEdgeOrderDoesNotMatter <- function() {
    gxlFile <- system.file("GXL/outOfOrderExample.gxl", package="graph")
    con <- file(gxlFile)
    tryCatch({
        g <- fromGXL(con)
        checkEquals(c("p", "v", "q", "w"), nodes(g))
    }, finally=close(con))
}
    

testUndirectedReading <- function() {
    gxlFileGz <- system.file("GXL/graphExample-02.gxl.gz", package="graph")
    con <- gzfile(gxlFileGz, open="rb")
    tryCatch({
        g <- fromGXL(con)
        checkEquals(10, numEdges(g))
        checkEquals(9, numNodes(g))
    }, finally=close(con))
}
