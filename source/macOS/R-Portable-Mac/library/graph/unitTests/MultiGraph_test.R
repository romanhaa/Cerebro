set.seed(0x12a9b)

## alias subsetEdgeSets; remove once it is exported


make_directed_MultiGraph <- function(use.factors = TRUE)
{
    ft1 <- data.frame(from=c("a", "a", "a", "b", "b"),
                        to=c("b", "c", "d", "a", "d"),
                      weight=c(1, 3.1, 5.4, 1, 2.2),
                      stringsAsFactors = use.factors)

    ft2 <- data.frame(from=c("a", "a", "a", "x", "x", "c"),
                        to=c("b", "c", "x", "y", "c", "a"),
                      weight=c(3.4, 2.6, 1, 1, 1, 7.9),
                      stringsAsFactors = use.factors)

    esets <- list(e1=ft1, e2=ft2)

    g <- MultiGraph(esets)

    list(esets=esets, g=g)
}

make_mixed_MultiGraph <- function(use.factors = TRUE)
{
    ft1 <- data.frame(from=c("a", "a", "a", "b", "b"),
                      to=c("b", "c", "d", "a", "d"),
                      weight=c(1, 3.1, 5.4, 1, 2.2),
                      stringsAsFactors = use.factors)

    ft2 <- data.frame(from=c("a", "a", "a", "x", "x", "c"),
                      to=c("b", "c", "x", "y", "c", "a"),
                      weight=c(3.4, 2.6, 1, 1, 1, 7.9),
                      stringsAsFactors = use.factors)

    ft3 <- data.frame(from=c("a", "a", "x", "x", "x"),
                      to  =c("b", "c", "a", "y", "c"),
                      weight=c(1:5),
                      stringsAsFactors = use.factors)

    esets <- list(e1=ft1, e2=ft2, e3=ft3, e4=ft2[FALSE, ],
                  e5=ft3[FALSE, ])

    g <- MultiGraph(esets, directed = c(TRUE, TRUE, FALSE, TRUE, FALSE))
    list(esets=esets, g=g)
}

randMultiGraph <- function(numNodes, numEdges)
{
    ftlist <- lapply(numEdges, function(ne) {
        graph:::randFromTo(numNodes, ne)
    })
    nn <- ftlist[[1L]]$nodes
    edgeSets <- lapply(ftlist, function(x) x[["ft"]])
    names(edgeSets) <- paste("e", seq_len(length(edgeSets)), sep="")
    MultiGraph(edgeSets, nodes = nn)
}

randFromTo2 <- function(numNodes, numEdges,
                            weightFun = function(N) rep(1L, N),
                            directed = TRUE)
{
    if (numNodes > 2^15) stop("too many nodes: ", numNodes)
    maxEdges <- numNodes * numNodes
    nodeNames <- sprintf("%010d", seq_len(numNodes))
    x <- c(rep(1L, numEdges), rep(0L, maxEdges - numEdges))
    idx <- which(sample(x) == 1L)
    to_i <- ((idx - 1L) %/% numNodes) + 1L
    from_i <- ((idx - 1L) %% numNodes) + 1L
    from <- nodeNames[from_i]
    to <- nodeNames[to_i]
    w <- weightFun(length(from))
    if (!directed) {
        tmp <- graph:::.mg_undirectEdges(from, to, w)
        from <- tmp$from
        to <- tmp$to
        w <- tmp$weight
        df <- data.frame(from = from, to = to, weight = w,
                         stringsAsFactors = FALSE)
        df <- df[!duplicated(df), ]
    } else {
        df <- data.frame(from = from, to = to, weight = w,
                         stringsAsFactors = FALSE)
    }
    list(nodes = nodeNames, ft = df)
}

sort_esets <- function(esets)
{
    ## sorting is based on column major ordering
    lapply(esets, function(ft) {
        ft <- ft[order(ft$to, ft$from), ]
    })
}

make_unique_ft <- function(ftdata)
{
    ## ftdata is a list with components $nodes and $ft
    ## $ft is a data.frame with columns 'from', 'to', and 'weight'
    ft <- ftdata[["ft"]]
    el <- paste(ft[["from"]], ft[["to"]], sep = "\t")
    dups <- duplicated(el)
    ftdata[["ft"]] <- ft[!dups, ]
    ftdata
}

test_basic_accessors <- function()
{
    basic <- make_mixed_MultiGraph()
    esets <- basic$esets
    g <- basic$g

    checkEquals(6L, numNodes(g))
    checkEquals(c("a", "b", "c", "d", "x", "y"), nodes(g))
    checkEquals(c(e1=5L, e2=6L, e3=5L, e4=0L, e5=0L), numEdges(g))
    checkEquals(structure(c(TRUE, TRUE, FALSE, TRUE, FALSE),
                          .Names = paste("e", 1:5, sep="")),
                isDirected(g))
}

test_no_edge_sets <- function()
{
    g1 <- MultiGraph(list(), nodes = letters)
    g2 <- MultiGraph(list(), nodes = letters, directed = FALSE)
    for (g in list(g1, g2)) {
        checkEquals(26L, numNodes(g))
        checkEquals(letters, nodes(g))
        checkEquals(list(), numEdges(g))
        checkEquals(list(), eweights(g))
        checkEquals(list(), eweights(g, "="))
        checkEquals(list(), isDirected(g))
    }
    tcon = textConnection(NULL, "w")
    sink(tcon)
    show(g1)
    sink()
    checkEquals("MultiGraph with 26 nodes and 0 edge sets",
                textConnectionValue(tcon))
    close(tcon)
}

test_create_empty_edgeSets <- function()
{
    df1 <- data.frame(from=c("a", "b"),
                       to=c("b", "c"), weight=c(1, 1))
    esets <- list(e1 = df1, empty1 = df1[FALSE, ])
    g <- MultiGraph(esets)
    checkEquals(c(e1=2L, empty1=0L), numEdges(g))

    dg <- MultiGraph(esets, directed = FALSE)
    checkEquals(c(e1=2L, empty1=0L), numEdges(dg))
}

test_edgeSets_arg_checking <- function()
{
    ## data.frame's in edgeSets list must have names:
    ## from, to, weights
    df0 <- data.frame(fr=c("a", "b"),
                      to=c("b", "c"), weights=c(1, 1))
    checkException(MultiGraph(list(e1=df0)))

    ## edgeSets must be named list or empty list
    checkException(MultiGraph(NULL))
    checkException(MultiGraph(list(data.frame(from=c("a", "b"),
                                              to=c("b", "c"),
                                              weights=c(1, 1)))))
}

test_no_nodes <- function()
{
    mg <- MultiGraph(list())
    checkEquals(0L, length(nodes(mg)))
    checkEquals(list(), numEdges(mg))
}

test_create_infer_nodes <- function()
{
    basic <- make_directed_MultiGraph()
    esets <- basic$esets
    g <- basic$g
    checkEquals(6L, numNodes(g))
    checkEquals(c(e1=5L, e2=6L), numEdges(g))
}

## test constructor with nodes arg given (singleton case)
test_bad_nodes_in_create <- function()
{
    basic <- make_directed_MultiGraph()
    esets <- basic$esets
    esets[[1]]$to <- as.character(esets[[1]]$to)
    bad_values <- c(NA, "a\n", "", "z|a", "a\t")
    for (v in bad_values) {
        tmp <- esets
        tmp[[1]]$to[3] <- v
        checkException(MultiGraph(tmp))
    }
}

test_dup_edges_is_an_error <- function()
{
    ## directed case
    ft1 <- data.frame(from=c("a", "a", "a"),
                        to=c("b", "c", "b"),
                      weight=c(1, 3.1, 5.4))
    checkException(MultiGraph(list(e1=ft1)))

    ## undirected case
    ft2 <- data.frame(from=c("a", "a", "b"),
                        to=c("b", "c", "a"),
                      weight=c(1, 3.1, 5.4))
    ## ok if directed
    junk <- MultiGraph(list(e1=ft2))
    checkException(MultiGraph(list(e1=ft2), directed = FALSE))
}

test_edgeWeights_create <- function()
{
    basic <- make_mixed_MultiGraph()
    esets <- basic$esets
    g <- basic$g

    esets <- sort_esets(esets)
    got <- eweights(g)
    checkIdentical(list(e1 = esets[[1L]][, "weight"]), got[1])
    checkIdentical(list(e2 = esets[[2L]][, "weight"]), got[2])
    ## undirected case normalizes edges by sorting, always putting the
    ## node that sorts first as from.
    checkIdentical(list(e3 = c(1L, 2L, 3L, 5L, 4L)), got[3])
}

test_edgeWeights_edge_names <- function()
{
    basic <- make_mixed_MultiGraph()
    esets <- basic$esets
    g <- basic$g

    wv <- eweights(g, names.sep = "=>")
    esets <- sort_esets(esets)

    want <- paste(esets[[1]]$from, esets[[1]]$to, sep = "=>")
    checkEquals(want, names(wv[[1]]))

    want <- c("a|b"=1L, "a|c"=2L, "a|x"=3L, "c|x"=5L, "x|y"=4L)
    checkIdentical(list(e3 = want), eweights(g, "|")[3])
}

test_supports_self_loops <- function()
{
    esets <- list(e1 = data.frame(from = c("a", "a"), to = c("a", "b"),
                  weight = c(1, 2)))
    g <- MultiGraph(esets)
    checkEquals(c(e1 = 2), numEdges(g))
}

test_isDirected <- function()
{
    g <- make_mixed_MultiGraph()$g
    checkEquals(c(e1=TRUE, e2=TRUE, e3=FALSE, e4=TRUE, e5=FALSE),
                isDirected(g))
}

test_ugraph_via_isDirected <- function()
{
    g <- make_mixed_MultiGraph()$g
    ## verify precondition
    want <- c(TRUE, TRUE, FALSE, TRUE, FALSE)
    names(want) <- paste("e", 1:5, sep="")
    checkEquals(want, isDirected(g))

    ug <- ugraph(g)
    want[] <- FALSE
    checkEquals(want, isDirected(ug))
}

test_ugraph_for_undirected_edge_sets <- function()
{
    df1 <- data.frame(from=c("x", "a", "b"),
                        to=c("a", "b", "x"),
                      weight=c(1, 2, 3))
    g <- MultiGraph(list(e1=df1), directed=FALSE)
    ug <- ugraph(g)
    checkEquals(nodes(g), nodes(ug))
    checkEquals(numEdges(g), numEdges(ug))
    checkEquals(isDirected(g), isDirected(ug))
    ## verify attributes have been dropped
    checkEquals(rep(1L, 3), eweights(ug)[[1]])
}

test_ugraph_for_directed_edge_sets <- function()
{
    df1 <- data.frame(from=c("x", "a", "b", "x", "b", "c"),
                        to=c("a", "x", "x", "b", "a", "x"),
                      weight=1:6)
    g <- MultiGraph(list(e1=df1), directed=TRUE)
    checkEquals(6, numEdges(g)[[1]])
    ug <- ugraph(g)
    checkEquals(nodes(g), nodes(ug))
    checkTrue(!isDirected(ug)[1])
    checkEquals(4, numEdges(ug)[[1]])
    checkEquals(rep(1L, 4), eweights(ug)[[1]])
    checkEquals(c("a=b", "a=x", "b=x", "c=x"),
                names(eweights(ug, "=")[[1]]))
}

mg_equals <- function(g1, g2)
{
    checkEquals(nodes(g1), nodes(g2))
    checkEquals(isDirected(g1), isDirected(g2))
    checkEquals(numEdges(g1), numEdges(g2))
    checkEquals(eweights(g1, "==>"), eweights(g2, "==>"))
}

test_edgeSetIntersect0_trivial <- function()
{
    ## Verify 0 and 1 edge set cases for directed/undirected
    df <- data.frame(from="a", to="b", weight=1L)
    mgs <- list(
                ## empty edge sets
                MultiGraph(list(), nodes = letters),
                MultiGraph(list(), nodes = letters, directed = FALSE),
                ## single edge set
                MultiGraph(list(e1=df)),
                MultiGraph(list(e1=df), directed = FALSE))
    for (g in mgs) {
        mg_equals(g, edgeSetIntersect0(g))
    }
    ## Verify empty intersection for disjoint graphs
    df1 <- data.frame(from="a", to="b", weight=1L)
    df2 <- data.frame(from="x", to="y", weight=1L)
    g <- MultiGraph(list(e1=df1, e2=df2))
    gu <- MultiGraph(list(e1=df1, e2=df2), directed = FALSE)
    want <- MultiGraph(list(), nodes = c("a", "b", "x", "y"))
    mg_equals(want, edgeSetIntersect0(g))
    mg_equals(want, edgeSetIntersect0(gu))
}
# 
test_edgeSetIntersect0_directed_1 <- function()
{
    ## non-trivial directed intersect
    g <- make_directed_MultiGraph()$g
    gi <- edgeSetIntersect0(g)
    ## TODO: do we want the minimal node set or not?
    ## checkEquals(c("a", "b", "c"), nodes(gi))
    checkEquals(nodes(g), nodes(gi))    # original nodes
    checkEquals(c(e1_e2=2L), numEdges(gi)[1L])
    checkEquals("e1_e2", names(numEdges(gi)))
    w <- c(as.numeric(NA), as.numeric(NA))
    names(w) <- c("a=>b", "a=>c")
    checkEquals(list(e1_e2=w), eweights(gi, "=>"))
}

test_edgeSetIntersect0_random <- function()
{
    make_data <- function(nsets, nn, ne, ns,
                          type=c("directed", "undirected")) {
        ## nsets: number of edge sets
        ## nn: number of nodes
        ## ne: number of edges
        ## ns: number of shared edges
        directed <- switch(match.arg(type),
                           directed=TRUE,
                           undirected=FALSE,
                           mixed=sample(c(TRUE, FALSE), nsets, replace=TRUE))
        grouped <- randFromTo2(nn, (ne * nsets) + ns, directed = all(directed))$ft
        ## for the undirected case, we will end up with fewer edges so
        ## need to adjust.
        ne <- (nrow(grouped) - ns) %/% nsets
        shared <- grouped[1:ns, ]
        starts <- seq(ns, nrow(grouped) - ne, by = ne) + 1L
        esets <- vector("list", nsets)
        names(esets) <- paste("e", 1:nsets, sep = "")
        for (i in seq_along(esets)) {
            z <- grouped[seq(starts[[i]], starts[[i]] + ne - 1L), ]
            z <- rbind(shared, z)
            esets[[i]] <- z
        }
        list(shared=shared,
             g=MultiGraph(esets, directed = directed),
             esets = esets)
    }
    do_test <- function(d)
    {
        gi <- edgeSetIntersect0(d$g)
        checkEquals(nrow(d$shared), numEdges(gi)[[1]])
        all_directed <- all(isDirected(d$g))
        checkEquals(all_directed, isDirected(gi)[[1]])
        checkEquals(nodes(d$g), nodes(gi))
    }

    for (t in c("directed", "undirected")) {
        for (i in 1:10) {
            do_test(make_data(2, 10, 10, 3, type = t))
            do_test(make_data(3, 10, 10, 3, type = t))
            do_test(make_data(3, 10, 10, 1, type = t))
            do_test(make_data(3, 11, 20, 6, type = t))
        }
    }
}

test_subSetEdgeSets_single <- function(){
    g <- make_directed_MultiGraph()$g
    gi <- subsetEdgeSets(g, "e1")
    checkEquals(nodes(g), nodes(gi))
    checkEquals(c(e1 = 5L), numEdges(gi)[1L])
    checkEquals("e1", names(numEdges(gi)))
    w <- c(1.0, 1.0, 3.1, 5.4, 2.2)
    names(w) <- c("b=>a", "a=>b", "a=>c", "a=>d", "b=>d")
    checkEquals(list(e1=w), eweights(gi, "=>"))
}

test_subSetEdgeSets_multiple <- function() {
    g <- make_mixed_MultiGraph()$g
    gi <- subsetEdgeSets(g, c("e1","e3"))
    checkEquals(nodes(g), nodes(gi))
    checkEquals(c(e1 = 5L, e3 = 5L), numEdges(gi))
    checkEquals(c("e1", "e3"), names(numEdges(gi)))
    w1 <- c( 1.0, 1.0, 3.1, 5.4, 2.2)
    w2 <- c(1, 2, 3, 5, 4)
    names(w1) <- c("b=>a", "a=>b", "a=>c", "a=>d", "b=>d")
    names(w2) <- c("a=>b", "a=>c", "a=>x", "c=>x", "x=>y")
    checkEquals(list(e1 = w1, e3 = w2), eweights(gi, "=>"))
}

test_subSetEdgeSets_no_duplicate_edgeSets <- function() {
    g <- make_directed_MultiGraph()$g
    checkException(subsetEdgeSets(g, c("e1", "e1")))
}

test_subSetEdgeSets_no_such_edgeSet <- function() {
    g <- make_directed_MultiGraph()$g
    checkException(subsetEdgeSets(g, "notAnEdgeSet"))
}

test_subSetEdgeSets_empty_edgeSet <- function() {
    g <- make_directed_MultiGraph()$g
    gi <- subsetEdgeSets(g, character(0))
    checkEquals(nodes(g), nodes(gi))
    checkEquals("list", class(numEdges(gi)))
    checkEquals(0, length(numEdges(gi)))
    checkEquals(character(0), names(numEdges(gi)))
}

test_extractFromTo_Directed <- function(use.factors=TRUE){
    ft1 <- data.frame(from=c("a", "a", "a", "b", "b"),
                      to=c("b", "c", "d", "a", "d"),
                      weight=c(1, 3.1, 5.4, 1, 2.2),
                      stringsAsFactors = use.factors)

    ft2 <- data.frame(from=c("a", "a", "a", "x", "x", "c"),
                      to=c("b", "c", "x", "y", "c", "a"),
                      weight=c(3.4, 2.6, 1, 1, 1, 7.9),
                      stringsAsFactors = use.factors)
    esets <- list(e1=ft1, e2=ft2)
    g <- MultiGraph(esets)
    res <- extractFromTo(g)
    ft1 <- ft1[do.call(order ,ft1["to"]),]
    rownames(ft1) <- 1:5
    ft2 <- ft2[do.call(order ,ft2["to"]),]
    rownames(ft2) <- 1:6
    checkEquals(list(e1 = ft1, e2 = ft2), res)
}

test_extractFromTo_UnDirected <- function(use.factors=TRUE){
    ft1 <- data.frame(from=c("a", "a", "a", "b", "b"),
                      to=c("b", "c", "d", "c", "d"),
                      weight=c(1, 3.1, 5.4, 1, 2.2),
                      stringsAsFactors = use.factors)

    ft2 <- data.frame(from=c("a", "a", "a", "x", "x", "c"),
                      to=c("b", "a", "x", "y", "c", "a"),
                      weight=c(3.4, 2.6, 1, 1, 1, 7.9),
                      stringsAsFactors = use.factors)
    esets <- list(e1=ft1, e2=ft2)

    g <- MultiGraph(esets,directed=c(FALSE,FALSE))
    res <- extractFromTo(g)

    ft1 <- ft1[do.call(order ,ft1["to"]),]
    rownames(ft1) <- 1:5

    ft2["from"] <- factor(c("a", "a", "a", "a", "c", "x"))
    ft2["to"] <- factor(c("a", "b", "c", "x", "x", "y"))
    ft2["weight"] <- c(2.6, 3.4, 7.9, 1, 1, 1)
    checkEquals(list(e1 = ft1, e2 = ft2), res)
}

test_degree_Mixed <- function(use.factors=TRUE){
    ft1 <- data.frame(from=c("a", "a", "a", "b", "b"),
                      to=c("b", "c", "d", "c", "d"),
                      weight=c(1, 3.1, 5.4, 1, 2.2),
                      stringsAsFactors = use.factors)

    ft2 <- data.frame(from=c("a", "a", "a", "x", "x", "c"),
                      to=c("b", "a", "x", "y", "c", "a"),
                      weight=c(3.4, 2.6, 1, 1, 1, 7.9),
                      stringsAsFactors = use.factors)
    esets <- list(e1=ft1, e2=ft2)
    g <- MultiGraph(esets,directed=c(FALSE,TRUE))
    deg <- degree(g)

    e1Degree  <- as.numeric(c(3, 3, 2, 2, 0, 0))
    attributes(e1Degree) <- list(names=c("a","b","c","d","x","y"))

    inDegree <- as.numeric(c(2, 1, 1, 0, 1, 1))
    attributes(inDegree) <- list(names=c("a","b","c","d","x","y"))
    outDegree <- as.numeric(c(3, 0, 1, 0, 2, 0))
    attributes(outDegree) <- list(names=c("a","b","c","d","x","y"))
    res <- list(e1 = e1Degree, e2 = list(inDegree = inDegree,
                               outDegree=outDegree))
    checkEquals(res, deg)
}

checkSubGraph <- function(g, subG) {
    nds <- nodes(g)
    subNodes <- nodes(subG)
    origFromTo <- extractFromTo(g)
    subFromTo <- extractFromTo(subG)
    sapply(names(origFromTo), function(x){
                indx <- (origFromTo[[x]]$from %in% subNodes) &
                (origFromTo[[x]]$to %in% subNodes)
                origdf = origFromTo[[x]]
                want <- origdf[(origdf$from %in% subNodes) & (origdf$to %in% subNodes),]
                subdf <- subFromTo[[x]]
                checkEquals(as.character(want$from), as.character(subdf$from))
                checkEquals(as.character(want$to), as.character(subdf$to))
                checkEquals(g@edge_sets[[x]]@weights[indx], subG@edge_sets[[x]]@weights)
            })
}

test_basic_subGraph <- function() {
    g <- make_mixed_MultiGraph()$g
    nds <- nodes(g)[1:3]
    sg <- subGraph(nds, g)
    checkSubGraph(g,sg)
}

test_large_subGraph <- function() {
    df1 <- graph:::randFromTo(1000L, 10001L, directed = TRUE,
                              weightFun = seq_len)
    df2 <- graph:::randFromTo(1000L, 10001L, directed = FALSE,
                              weightFun = seq_len)
    g <- MultiGraph(list(e1= df1$ft, e2 = df2$ft))
    nds <- sample( graph:::nodes(g), 100)
    subG <- subGraph(nds, g)
    checkSubGraph(g, subG)
}

test_basic_mgToGraphAM <- function() {
    g <- make_mixed_MultiGraph()$g
    res <- extractGraphAM(g)
    checkGraphAMObj(res, g)
}

test_large_mgToGraphAM  <- function() {
    df1 <- graph:::randFromTo(800L, 90L, directed = TRUE,
                              weightFun = seq_len)
    df2 <- graph:::randFromTo(800L, 60L, directed = FALSE,
                              weightFun = seq_len)
    g <- MultiGraph(list(e1= df1$ft, e2 = df2$ft))
    res <- extractGraphAM(g)
    checkGraphAMObj(res,g)
}

checkGraphAMObj <- function(am, mg){
    nds <- nodes(mg)
    dr <- isDirected(mg)
    checkEquals(names(mg@edge_sets),names(am))
    sapply(names(am), function(x){
                mat <- as(am[[x]], "matrix")
                checkEquals(colnames(mat),rownames(mat))
                checkEquals(colnames(mat), nds)
                wtMg <- graph:::edgeSetToMatrix(nds, mg@edge_sets[[x]], dr[[x]])
                checkEquals(mat, wtMg)
    })

}

test_mixed_MultiGraph_Intersect <- function(use.factors=TRUE) {

    ft1 <- data.frame(from=c("a", "a", "a", "b", "b"),
                      to=c("b", "c", "d", "a", "d"),
                      weight=c(1, 3.1, 5.4, 1, 2.2),
                      stringsAsFactors = use.factors)

    ft2 <- data.frame(from=c("a", "a", "a", "x", "x"),
                      to=c("b", "c", "x", "y", "c"),
                      weight=c(3.4, 2.6, 1, 1, 1),
                      stringsAsFactors = use.factors)

    ft3 <- data.frame(from=c("a", "a", "x", "x", "x"),
                      to  =c("b", "c", "a", "y", "c"),
                      weight=c(1:5),
                      stringsAsFactors = use.factors)

    esets <- list(e1=ft1, e2=ft2, e3=ft3, e4=ft2[FALSE, ],
                  e5=ft3[FALSE, ])

    g1 <- MultiGraph(esets, directed = c(TRUE, FALSE, TRUE, TRUE, FALSE))

    ft1 <- data.frame(from=c("a", "b"),
                      to=c("d", "d"),
                      weight=c(5.4, 2.2),
                      stringsAsFactors = use.factors)

    ft2 <- data.frame(from=c("a", "a", "a"),
                      to=c("b", "c", "x"),
                      weight=c(3.4, 2.6, 1),
                      stringsAsFactors = use.factors)

    esets <- list(e1=ft1, e2=ft2)

    g2 <- MultiGraph(esets, directed = c(TRUE, FALSE))
    res <- graphIntersect(g1, g2)
    checkEquals(nodes(res), c("a", "b", "c", "d", "x"))
    checkEquals(isDirected(res),
                structure(c(TRUE, FALSE), names = c("e1", "e2")))
    df <- extractFromTo(res)
    checkEquals(names(df), c("e1", "e2"))
    df1 <- data.frame(from = c("a", "b"), to = c("d", "d"), weight = c(5.4, 2.2))
    checkEquals(df$e1, df1)   
    df2 <- data.frame(from = c("a", "a", "a"), to = c("b", "c", "x"), weight = c(3.4, 2.6, 1))
    checkEquals(df$e2, df2)   
}

test_mixed_MultiGraph_Union <- function(use.factors=TRUE) {

    ft1 <- data.frame(from=c("a", "a", "a", "b", "b"),
                      to  =c("b", "c", "d", "a", "d"),
                      weight=c(1, 3.1, 5.4, 1, 2.2),
                      stringsAsFactors = use.factors)
    ft2 <- data.frame(from=c("a", "a"),
                      to=c("b", "c"),
                      weight=c(3.4, 2.6),
                      stringsAsFactors = use.factors)
    ft3 <- data.frame(from=c("a", "a"),
                      to  =c("d", "b"),
                      weight=c(1,2),
                      stringsAsFactors = use.factors)
    esets <- list(e1=ft1, e2=ft2, e3=ft3, e4=ft2[FALSE, ],
                  e5=ft3[FALSE, ])
    g1 <- MultiGraph(esets, directed = c(TRUE, FALSE, TRUE, TRUE, FALSE))

    ft1 <- data.frame(from=c("a", "a", "b"),
                      to=c("b", "x", "z"),
                      weight=c(6, 5, 2),
                      stringsAsFactors = use.factors)
    ft2 <- data.frame(from=c("a", "a", "a"),
                      to=c("a", "x", "y"),
                      weight=c(1, 2, 3),
                      stringsAsFactors = use.factors)
    esets <- list(e1=ft1, e2=ft2)
    g2 <- MultiGraph(esets, directed = c(TRUE, FALSE))
    
    res <- graphUnion(g1, g2)
    checkEquals(nodes(res), c("a", "b", "c", "d", "x", "y", "z"))
    checkEquals(names(res@edge_sets), c("e1", "e2", "e3", "e4", "e5"))
    checkEquals(isDirected(res), structure(c(TRUE, FALSE, TRUE, TRUE, FALSE), 
                                           names = c("e1", "e2", "e3", "e4", "e5")))
    
    df <- extractFromTo(res)
    checkEquals(names(df), c("e1", "e2", "e3", "e4", "e5"))
    df1 <- data.frame(from = c("b", "a", "a", "a", "b", "a", "b"), 
                      to = c("a", "b", "c", "d", "d", "x", "z"),
                      weight = c(1, NA, 3.1, 5.4, 2.2, 5.0, 2.0))
    checkEquals(df$e1, df1)   

    df2 <- data.frame(from = c("a", "a", "a", "a", "a"), 
                      to = c("a", "b", "c", "x", "y"),
                      weight = c(1, 3.4, 2.6, 2, 3 ))
    checkEquals(df$e2, df2)  

    df3 <- data.frame(from = c("a", "a"), 
                      to = c("b", "d"),
                      weight = c(2, 1))
    checkEquals(df$e3, df3)

    df4 <- data.frame(from = factor(), to = factor(), weight = numeric())
    checkEquals(df$e4, df4)
    checkEquals(df$e5, df4)
}
# 
test_MultiGraph_To_graphBAM <- function(use.factors=TRUE) {

    ft1 <- data.frame(from=c("a", "a", "a", "b", "b"),
                      to=c("b", "c", "d", "a", "d"),
                      weight=c(1, 3.1, 5.4, 1, 2.2),
                      stringsAsFactors = use.factors)
    
    ft2 <- data.frame(from=c("a", "a", "a", "x", "x"),
                      to=c("b", "c", "x", "y", "c"),
                      weight=c(3.4, 2.6, 1, 1, 1),
                      stringsAsFactors = use.factors)

    esets <- list(e1 = ft1, e2 = ft2, e3 = ft2[FALSE, ])
    g1 <- MultiGraph(esets, directed = c(TRUE, FALSE, TRUE))
    res <- extractGraphBAM(g1) 
    checkEquals(names(res), c("e1", "e2", "e3"))
  
    bam1 <- graphBAM(ft1, nodes=nodes(g1), edgemode = "directed")
    checkEquals(bam1, res$e1)

    bam2 <- graphBAM(ft2, nodes=nodes(g1), edgemode = "undirected")
    checkEquals(bam2, res$e2)

    bam3 <- graphBAM( ft2[FALSE, ], nodes=nodes(g1), edgemode = "directed")
    checkEquals(bam3, res$e3)
  
    res <- extractGraphBAM(g1, "e1") 
    checkEquals(bam1, res$e1)
    res <- extractGraphBAM(g1, c("e2", "e3"))
    target <- structure(list(bam2, bam3), names = c("e2", "e3"))
    checkEquals(target, res)
}  

test_MultiGraph_nodeAttributes <- function() {

    mg <- make_directed_MultiGraph()$g
    nds <- nodes(mg)
    checkException(nodeData(mg, n = c("a"), attr = "color"))

    checkException( nodeData(mg, n = "z", attr = "color") <- "red")
    nodeDataDefaults(mg, attr = "color") <- "violet"
    nodeDataDefaults(mg, attr = "class") <- "unknown"
    nodeData(mg, n = c("a", "x") , attr = "color") <- c("red", "green")
    nodeData(mg, attr = "class") <- "high"
    
    current  <- nodeData(mg, attr = "color")
    target <- structure(list("red", "violet", "violet", "violet", "green",
                             "violet") , names = nds)
    checkEquals(target, current)
    
    current <- nodeData(mg, attr = "class")
    target <- structure( as.list(rep("high",6)), names = nds)
    checkEquals(target, current)
    
    sg <- subGraph(c("a", "x"), mg)
    current <- nodeData(sg, attr = "color")
    target <- structure( list("red","green"), names = c("a", "x"))
    checkEquals(target, current)

    nodeDataDefaults(mg, attr="st") <- "unknown"
    nodeData(mg, n = c("b", "d"), attr = "st") <- mg@edge_sets$e1
    current <- nodeData(mg, attr = "st")
    target <- structure(list("unknown",  mg@edge_sets$e1, "unknown", 
                mg@edge_sets$e1, "unknown", "unknown"),
                   names = c( "a", "b", "c", "d", "x", "y"))
    checkEquals(target, current)

}


test_MultiGraph_edgeAttributes_directed  <- function() {
    ## Check directed edge
    mg <- make_mixed_MultiGraph()$g

    mgEdgeDataDefaults(mg, "e1", attr = "color") <- "violet"
    mgEdgeData(mg, "e1", from = c("a"), to = c("b"), attr = "color") <- "red"
    current  <- mgEdgeData(mg, "e1", attr = "color")
     
    nms <- paste( c("b", "a", "a", "a", "b"), c("a","b", "c", "d", "d"),sep ="|")
    target <- structure( list( "violet", "red", "violet", "violet", "violet"), names = nms)
    checkEquals(target, current)

    mgEdgeData(mg, "e1", to = "d",  attr = "color") <- "green"
    current <- mgEdgeData(mg, "e1", attr = "color") 
    target <- structure( list("violet", "red", "violet", "green", "green"), names = nms)
    checkEquals(target, current)
   
    mgEdgeData(mg, "e1", from = "b", attr = "color") <- c("pink")
    current <- mgEdgeData(mg, "e1", attr = "color") 
    target <- structure( list("pink", "red", "violet", "green" ,"pink"), names = nms)
    checkEquals(target, current)
   
    checkException(mgEdgeData(mg, "e1", attr = "class"))
    checkException(mgEdgeData(mg, "e9", attr = "color"))
 
}

test_MultiGraph_edgeAttributes_directed_S4  <- function() {
    ## Check directed edge
    mg <- make_mixed_MultiGraph()$g
    df <- extractFromTo(mg)$e1
    mgEdgeDataDefaults(mg, "e1", attr = "color") <- "unknown"
    mgEdgeData(mg, "e1", from = c("a"), to = c("b"), attr = "color") <- df
    current  <- mgEdgeData(mg, "e1", attr = "color")
     
    nms <- paste( c("b", "a", "a", "a", "b"), c("a","b", "c", "d", "d"),sep ="|")
    target <- structure( list( "unknown", df, "unknown", "unknown", "unknown"), names = nms)
    checkEquals(target, current)

    mgEdgeData(mg, "e1", to = "d",  attr = "color") <- matrix(1)
    current <- mgEdgeData(mg, "e1", attr = "color") 
    target <- structure( list("unknown", df, "unknown", matrix(1), matrix(1)), names = nms)
    checkEquals(target, current)
   
    mgEdgeData(mg, "e1", from = "b", attr = "color") <- matrix(0)
    current <- mgEdgeData(mg, "e1", attr = "color") 
    target <- structure( list(matrix(0), df, "unknown", matrix(1), matrix(0)), names = nms)
    checkEquals(target, current)
   
    checkException(mgEdgeData(mg, "e1", attr = "class"))
    checkException(mgEdgeData(mg, "e9", attr = "color"))
 
}

test_MultiGraph_edgeAttributes_undirected  <- function() {
    ## Check undirected edge
    mg <- make_mixed_MultiGraph()$g

    mgEdgeDataDefaults(mg, "e3", attr = "color") <- "unknown"
    mgEdgeData(mg, "e3", from = c("a"), to = c("b"), attr = "color") <- "red"
    current  <- mgEdgeData(mg, "e3", to = "a",  attr = "color")
    nms <- paste( c("b", "c", "x"), rep("a", 3), sep = "|")
    target <- structure(list("red", "unknown", "unknown"), names = nms)
    checkEquals(target, current)


    mgEdgeData(mg, "e3", to = "c",  attr = "color") <- "green"
    current <- mgEdgeData(mg, "e3", to = "c",  attr = "color") 
    nms <- paste(c("a", "x"), c("c", "c"), sep = "|")
    target <- structure( list("green", "green"), names = nms)
    checkEquals(target, current)
   
    mgEdgeData(mg, "e3", from = "b", attr = "color") <- c("pink")
    current <- mgEdgeData(mg, "e3", attr = "color") 
    nms <- paste(c("a", "a", "a", "c", "x", "b", "c", "x", "x", "y"),
            c("b", "c", "x", "x", "y", "a", "a", "a", "c", "x"), sep ="|")
    target <- structure( list("pink", "green", "unknown", "green", "unknown", 
                "pink", "green", "unknown", "green", "unknown"), names = nms)
    checkEquals(target, current)
}

test_MultiGraph_edgeAttributes_undirected_S4  <- function() {
    ## Check undirected edge
    mg <- make_mixed_MultiGraph()$g

    mgEdgeDataDefaults(mg, "e3", attr = "color") <- "unknown"
    mgEdgeData(mg, "e3", from = c("a"), to = c("b"), attr = "color") <- matrix(1)
    current  <- mgEdgeData(mg, "e3", to = "a",  attr = "color")
    nms <- paste( c("b", "c", "x"), rep("a", 3), sep = "|")
    target <- structure(list(matrix(1), "unknown", "unknown"), names = nms)
    checkEquals(target, current)

    df <- data.frame(1:2)
    mgEdgeData(mg, "e3", to = "c",  attr = "color") <- df
    current <- mgEdgeData(mg, "e3", to = "c",  attr = "color") 
    nms <- paste(c("a", "x"), c("c", "c"), sep = "|")
    target <- structure( list(df, df), names = nms)
    checkEquals(target, current)
   
    mgEdgeData(mg, "e3", from = "b", attr = "color") <- matrix(0)
    current <- mgEdgeData(mg, "e3", attr = "color") 
    nms <- paste(c("a", "a", "a", "c", "x", "b", "c", "x", "x", "y"),
            c("b", "c", "x", "x", "y", "a", "a", "a", "c", "x"), sep ="|")
    target <- structure( list(matrix(0), df, "unknown", df, "unknown", matrix(0),
                    df, "unknown", df, "unknown"), names = nms)
    checkEquals(target, current)
}


test_MultiGraph_edgeAttributes_empty <- function() {

    mg <- make_mixed_MultiGraph()$g
    checkException(mgEdgeData(mg, "e4", from ="a", attr = "color") <- "green")
    checkException(mgEdgeData(mg, "e4", attr = "color"))
}

test_MultiGraph_edgeAttributes_subGraph <- function() {

    mg <- make_mixed_MultiGraph()$g
    mgEdgeDataDefaults(mg, "e1", attr = "color") <- "violet"
    mgEdgeData(mg, "e1", attr = "color") <- "red"
    mgEdgeData(mg, "e1", from = c("a", "a"), to = c("b", "c"), attr = "color") <- c("yellow", "pink")
    
    current <- mgEdgeData(mg,"e1", attr ="color")
    nms <- paste(c("b", "a", "a", "a", "b"), c("a", "b", "c", "d", "d"), sep = "|")
    target <- structure(list("red", "yellow", "pink", "red", "red"), names = nms)
    checkEquals(target, current)
 
    sg <- subGraph(c("a","c", "d"), mg)
    current <- mgEdgeData(sg, "e1", attr = "color")
    nms <- paste(c("a", "a"), c("c", "d"), sep ="|")
    target <- structure( list("pink", "red"), names = nms)
    checkEquals(target, current)
}

test_MultiGraph_Intersection_Attributes <- function(use.factors=TRUE){

    setClass("myType", representation = representation(typ ="character")) 
    myType <- function(typ){ new("myType", typ = typ)}

    colorFun <- function(x,y) {
        if(x@typ =="low" || y@typ == "med")
            return("low")
        else
            return("high")
    }
    typeFun <- function(x,y) {
        if(is(x, "myType")  && is(y, "myType")){
            if(x@typ =="low" || y@typ == "med")
                return("low")
            else
                return("high")
        }
        else {return (NA)}

    }
    edgeFun <- structure(rep(list(structure( list(colorFun), names = "myType")),2),
                         names = c("e1", "e2"))

    ft1 <- data.frame(from=c("a", "b", "b", "f"),
                      to=c("b", "c", "d", "g"),
                      weight=c(1, 2, 3, 4),
                      stringsAsFactors = use.factors)

    ft2 <- data.frame(from =c("e", "f", "a", "b", "f"),
                      to =c( "a", "a", "b", "c", "c"),
                      weight =c(2, 3, 4, 5, 9),
                      stringsAsFactors= use.factors)
    esets <- list(e1=ft1, e2=ft2)
    g1 <- MultiGraph(esets)
    
    mgEdgeDataDefaults(g1, "e1", attr = "color") <- "violet"
    mgEdgeData(g1, "e1", from = c("b", "b", "f"), 
               to = c("c", "d", "g"),
               attr = "color") <- c ("red", "green", "green")
    
    mgEdgeDataDefaults(g1, "e2", attr = "color") <- "violet"
    mgEdgeData(g1, "e2", from = c("b", "f"), 
               to = c("c", "c"),
               attr = "color") <- c("red", "green")
    mgEdgeDataDefaults(g1, "e1", attr = "myType") <- "missing"
    mgEdgeData(g1, "e1", from = c("b", "f"), 
               to = c("c", "g"),
               attr = "myType") <- c(myType("low"), myType("high"))
    
    mgEdgeDataDefaults(g1, "e2", attr = "myType") <- "missing"
    mgEdgeData(g1, "e2", from = c("a", "f"), 
               to = c("b", "c"),
               attr = "myType") <-  c(myType("high"), myType("low"))

    ft1 <- data.frame(from=c("b", "f"),
                      to=c("c", "g"),
                      weight=c(2, 3),
                      stringsAsFactors = use.factors)

    ft2 <- data.frame(from =c("a", "f"),
                      to =c( "b", "c"),
                      weight =c(4, 6),
                      stringsAsFactors= use.factors)
    esets <- list(e1=ft1, e2=ft2)
    g2 <- MultiGraph(esets)
    mgEdgeDataDefaults(g2, "e1", attr = "color") <- "violet"
    mgEdgeData(g2, "e1", from = c("b", "f"), 
               to = c("c", "g"),
               attr = "color") <- c ("red", "green")
    mgEdgeDataDefaults(g2, "e2", attr = "color") <- "violet"
    mgEdgeData(g2, "e2", from = c("a", "f"), 
               to = c("b", "c"), 
               attr = "color") <- c("red", "green")
   
    mgEdgeDataDefaults(g2, "e1", attr="myType") <- "unknown"
    mgEdgeData(g2, "e1", from = c("b", "f"), 
               to = c("c", "g"),
               attr = "myType") <- c(myType("med"), myType("high"))
    mgEdgeDataDefaults(g2, "e2", attr="myType") <- "unknown"
    mgEdgeData(g2, "e2", from = c("a", "f"), 
               to = c("b", "c"),
               attr = "myType") <- c(myType("high"), myType("med"))
    nodeDataDefaults(g1, attr = "color") <- "violet"
    nodeDataDefaults(g2, attr = "color") <- "violet"
    nodeDataDefaults(g1, attr= "type") <- "unknown"
    nodeDataDefaults(g2, attr = "type") <- "unknown"
    nodeData(g1,n = c("a", "b", "c"), attr ="color") <- c("red", "green", "blue")
    nodeData(g1,n = c("b", "c"), attr ="type") <- c(myType("low"), myType("high"))
    nodeData(g2,n = c("a", "b", "c"), attr ="color") <- c("red", "green", "red")
    nodeData(g2,n = c("b", "c"), attr ="type") <- c(myType("med"), myType("low"))

    res <- graphIntersect(g1, g2, nodeFun = list(type =typeFun), edgeFun = edgeFun)

    current <- mgEdgeData(res,"e1", attr = "weight")
    target <- structure(list(2, as.numeric(NA)), names = paste(c("b", "f"), c("c", "g"), sep="|"))
    checkEquals(target, current)

    current <- mgEdgeData(res,"e1", attr = "color")
    target <- structure(list("red", "green"), names = paste(c("b", "f"), c("c", "g"), sep="|"))
    checkEquals(target, current)

    current <- mgEdgeData(res,"e2", attr = "weight")
    target <- structure(list(4, as.numeric(NA)), names = paste(c("a", "f"), c("b", "c"), sep="|"))
    checkEquals(target, current)

    current <- mgEdgeData(res,"e2", attr = "color")
    target <- structure(list(as.character(NA), "green"), 
                        names = paste(c("a", "f"), c("b", "c"), sep="|"))
    checkEquals(target, current)

    current <- mgEdgeData(res,"e1", attr = "myType")
    target <- structure(list("low", "high"), 
                        names = paste(c("b", "f"), c("c", "g"), sep="|"))
    checkEquals(target, current)

    current <- mgEdgeData(res,"e2", attr = "myType")
    target <- structure(list("high", "low"), 
                        names = paste(c("a", "f"), c("b", "c"), sep="|"))
    checkEquals(target, current)


    nodeColor <- nodeData(res, attr = "color")
    target <-  as.list(structure(c("red", "green", NA, "violet", "violet"), 
                                 names = c("a", "b", "c", "f", "g")))
    checkEquals(target, nodeColor)

    nodeType <- nodeData(res, attr = "type")
    target <-  as.list(structure(c("unknown", "low", "high", "unknown",
                                   "unknown"), 
                                 names = c("a", "b", "c", "f", "g")))
    checkEquals(target, nodeType)

}

test_MultiGraph_Union_Attributes <- function(use.factors=TRUE){

    setClass("myType", representation = representation(typ ="character")) 
    myType <- function(typ){ new("myType", typ = typ)}
    typeFun <- function(x,y) {
        if(is(x, "myType")  && is(y, "myType")){
            if(x@typ =="low" || y@typ == "med")
                return("low")
            else
                return("high")
        }
        else {return (NA)}

    }

    funList <- structure(rep(list(structure( list(typeFun), names = "myType")),2),
                         names = c("e1", "e2"))

    ft1 <- data.frame(from=c("a", "b", "b", "f"),
                      to=c("b", "c", "d", "g"),
                      weight=c(1, 2, 3, 4),
                      stringsAsFactors = use.factors)

    ft2 <- data.frame(from =c("e", "f", "a", "b", "f"),
                      to =c( "a", "a", "b", "c", "c"),
                      weight =c(2, 3, 4, 5, 9),
                      stringsAsFactors= use.factors)
    esets <- list(e1=ft1, e2=ft2)
    g1 <- MultiGraph(esets)
    
    mgEdgeDataDefaults(g1, "e1", attr = "color") <- "violet"
    mgEdgeData(g1, "e1", from = c("b", "b", "f"), to = c("c", "d", "g"),
               attr = "color") <- c ("red", "green", "green")
    
    mgEdgeDataDefaults(g1, "e2", attr = "color") <- "violet"
    mgEdgeData(g1, "e2", from = c("b", "f"), 
               to = c("c", "c"),
               attr = "color") <- c("red", "green")

    mgEdgeDataDefaults(g1, "e1", attr ="myType") <- "unknown"
    mgEdgeData(g1, "e1", from = c("b", "f"), 
               to = c("c", "g"),
               attr = "myType") <- c(myType("low"), myType("high"))
    
    mgEdgeDataDefaults(g1, "e2", attr ="myType") <- "unknown"
    mgEdgeData(g1, "e2", from = c("a", "f"), 
               to = c("b", "c"),
               attr = "myType") <-  c(myType("high"), myType("low"))


    ft1 <- data.frame(from=c("b", "f"),
                      to=c("c", "g"),
                      weight=c(2, 3),
                      stringsAsFactors = use.factors)

    ft2 <- data.frame(from =c("a", "f"),
                      to =c( "b", "c"),
                      weight =c(4, 6),
                      stringsAsFactors= use.factors)
    esets <- list(e1=ft1, e2=ft2)
    g2 <- MultiGraph(esets)
    mgEdgeDataDefaults(g2, "e1", attr = "color") <- "violet"
    mgEdgeData(g2, "e1", from = c("b", "f"), 
               to = c("c", "g"),
               attr = "color") <- c ("red", "green")
    
    mgEdgeDataDefaults(g2, "e2", attr = "color") <- "violet"
    mgEdgeData(g2, "e2", from = c("a", "f"), 
               to = c("b", "c"), 
               attr = "color") <- c("red", "green")
    
    mgEdgeDataDefaults(g2, "e1", attr = "myType") <- "unknown"
    mgEdgeData(g2, "e1", from = c("b", "f"), 
               to = c("c", "g"),
               attr = "myType") <- c(myType("med"), myType("high"))

    mgEdgeDataDefaults(g2, "e2", attr = "myType") <- "unknown"
    mgEdgeData(g2, "e2", from = c("a", "f"), 
               to = c("b", "c"),
               attr = "myType") <- c(myType("high"), myType("med"))

    res <- graphUnion(g1, g2, edgeFun = funList)

    current <- mgEdgeData(res,"e1", attr = "weight")
    target <- structure(list(1, 2, 3, as.numeric(NA)), 
                        names = paste(c("a", "b", "b", "f"), c( "b", "c", "d", "g"), sep="|"))
    checkEquals(target, current)

    current <- mgEdgeData(res,"e1", attr = "color")
    target <- structure(list("violet", "red", "green", "green"), 
                        names = paste(c("a", "b", "b", "f"), c("b", "c", "d", "g"), sep="|"))
    checkEquals(target, current)

    current <- mgEdgeData(res,"e2", attr = "weight")
    target <- structure(list(2, 3, 4, 5, as.numeric(NA)),
                        names = paste(c("e", "f", "a", "b", "f"), c("a", "a", "b", "c", "c"), sep="|"))
    checkEquals(target, current)

    current <- mgEdgeData(res,"e2", attr = "color")
    target <- structure(list( "violet", "violet", as.character(NA), "red", "green"), 
                        names = paste(c("e", "f", "a", "b", "f"), c("a", "a", "b", "c", "c"), sep="|"))
    checkEquals(target, current)

    current <- mgEdgeData(res,"e1", attr = "myType")
    target <- structure(list("unknown", "low", "unknown", "high"), 
                        names = paste(c("a", "b", "b", "f"), c("b", "c", "d", "g"), sep="|"))
    checkEquals(target, current)

    current <- mgEdgeData(res,"e2", attr = "myType")
    target <- structure(list("unknown", "unknown", "high", "unknown","low"), 
                        names = paste(c("e", "f", "a", "b", "f"), c("a", "a", "b", "c", "c"), sep="|"))
    checkEquals(target, current)
}


test_MultiGraph_nodeUnion_Attributes <- function(use.factors=TRUE){

    setClass("myType", representation = representation(typ ="character")) 
    myType <- function(typ){ new("myType", typ = typ)}
    testFun <- function(x,y) {
        if(is(x, "myType")  && is(y, "myType")){
            if(x@typ =="aa" || y@typ == "ac")
                return("ax")
            else 
                return("notype")
        }
        else
            return("ab")
    }
    funList <- structure(list(testFun), names ="gene")

    ft1 <- data.frame(from=c("a", "a", "a", "b", "b"),
                      to  =c("b", "c", "d", "a", "d"),
                      weight=c(1, 3.1, 5.4, 1, 2.2),
                      stringsAsFactors = use.factors)
    ft2 <- data.frame(from=c("a", "a"),
                      to=c("b", "c"),
                      weight=c(3.4, 2.6),
                      stringsAsFactors = use.factors)
    ft3 <- data.frame(from=c("a", "a"),
                      to  =c("d", "b"),
                      weight=c(1,2),
                      stringsAsFactors = use.factors)
    esets <- list(e1=ft1, e2=ft2, e3=ft3, e4=ft2[FALSE, ],
                  e5=ft3[FALSE, ])
    g1 <- MultiGraph(esets, directed = c(TRUE, FALSE, TRUE, TRUE, FALSE))
   
    nodeDataDefaults(g1, attr = "color") <- "violet" 
    nodeData(g1, n = c("a", "b", "c") , attr = "color") <- c("red", "green", "blue")
    nodeDataDefaults(g1, attr = "type") <- "unknown"
    nodeData(g1, n = c("a", "b"), attr = "type") <- c("low", "high")
    nodeDataDefaults(g1, attr = "kp") <- "kpunknown"
    nodeData(g1, n = c("a", "b"), attr = "kp") <- c("kplow", "kphigh")
    nodeDataDefaults(g1, attr = "gene") <- "XX"
    nodeData(g1, n = c("a", "b"), attr = "gene") <- c(myType("aa"), myType("bt"))

    ft1 <- data.frame(from=c("a", "a", "b"),
                      to=c("b", "x", "z"),
                      weight=c(6, 5, 2),
                      stringsAsFactors = use.factors)

    ft2 <- data.frame(from=c("a", "a", "a"),
                      to=c("a", "x", "y"),
                      weight=c(1, 2, 3),
                      stringsAsFactors = use.factors)

    esets <- list(e1=ft1, e2=ft2)
    g2 <- MultiGraph(esets, directed = c(TRUE, FALSE))
    nodeDataDefaults(g2, attr ="color") <- "violet"
    nodeData(g2, n = c("a", "b", "x", "y", "z") , attr = "color") <- c("red", "red", "green", "pink", "yellow")
    nodeDataDefaults(g2, attr ="type") <- "unknown"
    nodeData(g2, n = c("a", "b"), attr = "type") <- c("low", "high")
    nodeDataDefaults(g2, attr ="gene") <- "XX"
    nodeData(g2, n = c("a", "b"), attr = "gene") <- c(myType("at"), myType("kt"))

    res <- graphUnion(g1, g2, nodeFun = funList )

    current <- nodeData(res, attr = "color")
    target <- as.list( structure(c("red", NA, "blue", "violet", "green", "pink", "yellow"), 
                                 names = c("a", "b", "c", "d", "x", "y", "z")))
    checkEquals(target, current)

    current <- nodeData(res, attr = "type")
    target <- as.list( structure(c("low", "high", "unknown", "unknown",
                                   "unknown", "unknown", "unknown"), 
                                 names = c("a", "b", "c", "d", "x", "y", "z")))
    checkEquals(target, current)

    current <- nodeData(res, attr = "kp")
    cn <- "kpunknown"
    target <- as.list( structure(c("kplow", "kphigh", cn, cn, cn, cn, cn), 
                                 names = c("a", "b", "c", "d", "x", "y", "z")))
    checkEquals(target, current)
    current <- nodeData(res, n =c("a", "b"), attr ="gene")
    target <- as.list( structure(c("ax", "notype"), 
                                 names = c("a", "b")))
    checkEquals(target, current)

    current <- nodeData(res, n = c("c", "d", "x", "y", "z"), attr ="gene")
    target <- as.list( structure(c( "XX", "XX", "XX", "XX", "XX"), 
                                 names = c("c", "d", "x", "y", "z")))
    checkEquals(target, current)

}

