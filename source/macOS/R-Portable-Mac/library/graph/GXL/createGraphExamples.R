gxlFiles <- list.files(pattern="graphExample-[0-9]+\\.gxl\\.gz$")
graphExamples <- list()
i <- 1
for (gf in gxlFiles) {
    con <- gzfile(gf, open="rb")
    graphExamples[[i]] <- graph::fromGXL(con)
    i <- i + 1
    close(con)
}
    
