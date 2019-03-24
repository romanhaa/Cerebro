## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval = FALSE---------------------------------------------------------
#  chooseCRANmirror()
#  install.packages("BiocManager")

## ---- eval = FALSE---------------------------------------------------------
#  BiocManager::install(c("GenomicRanges", "Organism.dplyr"))

## ---- eval = FALSE---------------------------------------------------------
#  BiocManager::install()

## --------------------------------------------------------------------------
BiocManager::version()

## --------------------------------------------------------------------------
BiocManager::valid()

## --------------------------------------------------------------------------
avail <- BiocManager::available()
length(avail)
BiocManager::available("BSgenome.Hsapiens")

## ---- eval = FALSE---------------------------------------------------------
#  BiocManager::install(version="3.7")

## ---- eval = FALSE---------------------------------------------------------
#  .libPaths()

