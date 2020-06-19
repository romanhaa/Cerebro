test_ES_subset <- function() {
    data(sample.ExpressionSet)
    eset <- sample.ExpressionSet
    map <- getAnnMap("SYMBOL", annotation(eset))

    gss <- getBroadSets(system.file("extdata", "Broad.xml",
                                    package="GSEABase"))
    es <- eset[gss[[2]],]

    ids <- mget(geneIds(gss[[2]]), revmap(map), ifnotfound=NA)
    ids <- ids[!is.na(ids)]

    checkEquals(sum(unique(unlist(ids)) %in% featureNames(eset)),
                unname(nrow(es)))
    m <- mapIdentifiers(gss[[2]], AnnotationIdentifier(annotation(es)))
    checkTrue(all(featureNames(es) %in% geneIds(m)))
    checkTrue(all.equal(es, sample.ExpressionSet[m,]))
}
