## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = interactive())

## ---- eval = FALSE-------------------------------------------------------
#  chooseCRANmirror()
#  install.packages("BiocManager")

## ---- eval = FALSE-------------------------------------------------------
#  BiocManager::install(c("GenomicRanges", "Organism.dplyr"))

## ---- eval = FALSE-------------------------------------------------------
#  BiocManager::install()

## ------------------------------------------------------------------------
#  BiocManager::version()

## ------------------------------------------------------------------------
#  BiocManager::valid()

## ------------------------------------------------------------------------
#  avail <- BiocManager::available()
#  length(avail)                               # all CRAN & Bioconductor packages
#  BiocManager::available("BSgenome.Hsapiens") # BSgenome.Hsapiens.* packages

## ---- eval = FALSE-------------------------------------------------------
#  BiocManager::install()

## ---- eval = FALSE-------------------------------------------------------
#  BiocManager::valid()

## ---- eval = FALSE-------------------------------------------------------
#  BiocManager::install(version="3.7")

## ---- eval = FALSE-------------------------------------------------------
#  .libPaths()

## ---- eval = FALSE-------------------------------------------------------
#  options(
#      repos = "file:///path/to/CRAN-mirror",
#      BioC_mirror = "file:///path/to/Bioc-mirror"
#  )

## ---- eval = FALSE-------------------------------------------------------
#  options(
#      BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = FALSE
#  )

## ---- eval = FALSE-------------------------------------------------------
#  install.package(c("BiocManager", "BiocVersion"))

## ---- eval = TRUE--------------------------------------------------------
sessionInfo()

