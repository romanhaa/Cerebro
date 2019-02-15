## ----setup, include = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
# increase the screen width
options(width = 90)
# reduce the minimum number of characters for the tibble column titles
options(pillar.min_title_chars = 8)

## ----install-package, eval=FALSE--------------------------------------------------------
#  install.packages("msigdbr")

## ----load-package, message=FALSE--------------------------------------------------------
library(msigdbr)

## ----show-species-----------------------------------------------------------------------
msigdbr_show_species()

## ----get-human-all----------------------------------------------------------------------
m_df = msigdbr(species = "Homo sapiens")
head(m_df)

## ----get-mouse-h------------------------------------------------------------------------
m_df = msigdbr(species = "Mus musculus", category = "H")
head(m_df)

## ----get-mouse-c2-----------------------------------------------------------------------
m_df = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
head(m_df)

## ----get-mouse-h-filter-----------------------------------------------------------------
m_df = msigdbr(species = "Mus musculus") %>% dplyr::filter(gs_cat == "H")
head(m_df)

## ----cp-entrez, eval=FALSE--------------------------------------------------------------
#  m_t2g = m_df %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
#  enricher(gene = genes_entrez, TERM2GENE = m_t2g, ...)

## ----cp-symbols, eval=FALSE-------------------------------------------------------------
#  m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
#  enricher(gene = genes_symbols, TERM2GENE = m_t2g, ...)

## ----fgsea, eval=FALSE------------------------------------------------------------------
#  m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#  fgsea(pathways = m_list, ...)

