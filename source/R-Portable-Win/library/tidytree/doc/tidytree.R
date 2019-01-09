## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
		   message = FALSE)

## ----echo=FALSE, results="hide", message=FALSE---------------------------
library("ape")
library("tidytree")

## ------------------------------------------------------------------------
library(ape)
set.seed(2017)
tree <- rtree(4)
tree
x <- as_data_frame(tree)
x

## ------------------------------------------------------------------------
as.phylo(x)

## ------------------------------------------------------------------------
d <- tibble(label = paste0('t', 1:4),
            trait = rnorm(4))

y <- full_join(x, d, by = 'label')

## ------------------------------------------------------------------------
as.treedata(y)

## ------------------------------------------------------------------------
y %>% as.treedata %>% as_data_frame

## ------------------------------------------------------------------------
child(y, 5)
parent(y, 2)
offspring(y, 5)
ancestor(y, 2)
sibling(y, 2)
MRCA(y, 2, 3)

## ------------------------------------------------------------------------
nwk <- '(((((((A:4,B:4):6,C:5):8,D:6):3,E:21):10,((F:4,G:12):14,H:8):13):13,((I:5,J:2):30,(K:11,L:11):2):17):4,M:56);'
tree <- read.tree(text=nwk)

groupClade(as_data_frame(tree), c(17, 21))

## ------------------------------------------------------------------------
## the input nodes can be node ID or label
groupOTU(x, c('t1', 't4'), group_name = "fake_group")

## ------------------------------------------------------------------------
cls <- list(c1=c("A", "B", "C", "D", "E"),
            c2=c("F", "G", "H"),
            c3=c("L", "K", "I", "J"),
            c4="M")

as_data_frame(tree) %>% groupOTU(cls)

