## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
		   message = FALSE)

## ----echo=FALSE, results="hide", message=FALSE---------------------------
library('jsonlite')
library("treeio")

## ----comment=NA----------------------------------------------------------
nhxfile <- system.file("extdata/NHX", "phyldog.nhx", package="treeio")
nhx <- read.nhx(nhxfile)
# write.beast(nhx, file = "phyldog.tree")
write.beast(nhx)

## ----comment=NA----------------------------------------------------------
mlcfile <- system.file("extdata/PAML_Codeml", "mlc", package="treeio")
ml <- read.codeml_mlc(mlcfile)
# write.beast(ml, file = "codeml.tree")
write.beast(ml)

## ----comment=NA----------------------------------------------------------
phylo <- as.phylo(nhx)
## print the newick text
write.tree(phylo)

N <- Nnode2(phylo)
fake_data <- data_frame(node = 1:N, fake_trait = rnorm(N), another_trait = runif(N))
fake_tree <- treedata(phylo = phylo, data = fake_data)
write.beast(fake_tree)

## ------------------------------------------------------------------------
## combine tree object with data
tree_with_data <- full_join(nhx, fake_data, by = "node")
tree_with_data

## merge two tree object
tree2 <- merge_tree(nhx, fake_tree)
tree2

identical(tree_with_data, tree2)

## ----comment=NA----------------------------------------------------------
write.beast(tree2)

## ------------------------------------------------------------------------
outfile <- tempfile(fileext = ".tree")
write.beast(tree2, file = outfile)
read.beast(outfile)

## ----comment=NA----------------------------------------------------------
write.jtree(tree2)

## ----comment=NA----------------------------------------------------------
jtree_file <- tempfile(fileext = '.jtree')
write.jtree(tree2, file = jtree_file)
jsonlite::fromJSON(jtree_file)

## ------------------------------------------------------------------------
read.jtree(jtree_file)

