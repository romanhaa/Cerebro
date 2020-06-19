## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
		   message = FALSE)

## ----echo=FALSE, results="hide", message=FALSE---------------------------
library(tidyr)
library(dplyr)
library(tidytree)
library(ggplot2)

library("treeio")

CRANpkg <- function (pkg) {
    cran <- "https://CRAN.R-project.org/package"
    fmt <- "[%s](%s=%s)"
    sprintf(fmt, pkg, cran, pkg)
}

Biocpkg <- function (pkg) {
    sprintf("[%s](http://bioconductor.org/packages/%s)", pkg, pkg)
}

## ----treeio-function, echo=F, message=FALSE------------------------------
ff <- matrix(c(
  "read.beast"      , "parsing output of BEAST",
  "read.codeml"     , "parsing output of CodeML (rst and mlc files)",
  "read.codeml_mlc" , "parsing mlc file (output of CodeML)",
  "read.hyphy"      , "parsing output of HYPHY",
  "read.jplace"     , "parsing jplace file including output of EPA and pplacer",
  "read.mrbayes"    , "parsing output of MrBayes",
  "read.newick"     , "parsing newick string, with ability to parse node label as support values",
  "read.nhx"        , "parsing NHX file including output of PHYLDOG and RevBayes",
  "read.paml_rst"   , "parsing rst file (output of BaseML or CodeML)",
  "read.phylip"     , "parsing phylip file (phylip alignment + newick string)",
  "read.r8s"        , "parsing output of r8s",
  "read.raxml"      , "parsing output of RAxML"
  ), ncol=2, byrow=TRUE)
ff <- as.data.frame(ff)
colnames(ff) <- c("Parser function", "Description")
knitr::kable(ff, caption = "Parser functions defined in treeio", booktabs = T)

## ------------------------------------------------------------------------
file <- system.file("extdata/BEAST", "beast_mcc.tree", package="treeio")
beast <- read.beast(file)
beast

## ------------------------------------------------------------------------
get.fields(beast)

## ------------------------------------------------------------------------
file <- system.file("extdata/MrBayes", "Gq_nxs.tre", package="treeio")
read.mrbayes(file)

## ----fig.width=12, fig.height=10, warning=FALSE, fig.align="center"------
brstfile <- system.file("extdata/PAML_Baseml", "rst", package="treeio")
brst <- read.paml_rst(brstfile)
brst

## ------------------------------------------------------------------------
crstfile <- system.file("extdata/PAML_Codeml", "rst", package="treeio")
## type can be one of "Marginal" or "Joint"
crst <- read.paml_rst(crstfile, type = "Joint")
crst

## ------------------------------------------------------------------------
mlcfile <- system.file("extdata/PAML_Codeml", "mlc", package="treeio")
mlc <- read.codeml_mlc(mlcfile)
mlc

## ------------------------------------------------------------------------
## tree can be one of "rst" or "mlc" to specify
## using tree from which file as base tree in the object
ml <- read.codeml(crstfile, mlcfile, tree = "mlc")
ml

## ----warning=FALSE-------------------------------------------------------
ancseq <- system.file("extdata/HYPHY", "ancseq.nex", package="treeio")
read.hyphy.seq(ancseq)

## ----warning=FALSE-------------------------------------------------------
nwk <- system.file("extdata/HYPHY", "labelledtree.tree", package="treeio")
tipfas <- system.file("extdata", "pa.fas", package="treeio")
hy <- read.hyphy(nwk, ancseq, tipfas)
hy

## ----fig.width=4, fig.height=6, width=60, warning=FALSE, fig.align="center"----
r8s <- read.r8s(system.file("extdata/r8s", "H3_r8s_output.log", package="treeio"))
r8s

## ----fig.width=12, fig.height=10, width=60, warning=FALSE, fig.align="center"----
raxml_file <- system.file("extdata/RAxML", "RAxML_bipartitionsBranchLabels.H3", package="treeio")
raxml <- read.raxml(raxml_file)
raxml

## ------------------------------------------------------------------------
nhxfile <- system.file("extdata/NHX", "phyldog.nhx", package="treeio")
nhx <- read.nhx(nhxfile)
nhx

## ------------------------------------------------------------------------
phyfile <- system.file("extdata", "sample.phy", package="treeio")
phylip <- read.phylip(phyfile)
phylip

## ------------------------------------------------------------------------
jpf <- system.file("extdata/EPA.jplace",  package="treeio")
jp <- read.jplace(jpf)
print(jp)

## ------------------------------------------------------------------------
jtree_file <- tempfile(fileext = '.jtree')
write.jtree(beast, file = jtree_file)
read.jtree(file = jtree_file)

## ------------------------------------------------------------------------
x <- data_frame(label = as.phylo(beast)$tip.label, trait = rnorm(Ntip(beast)))
full_join(beast, x, by="label")

N <- Nnode2(beast)
y <- data_frame(node = 1:N, fake_trait = rnorm(N), another_trait = runif(N))
full_join(beast, y, by="node")

## ------------------------------------------------------------------------
beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
rst_file <- system.file("examples/rst", package="ggtree")
mlc_file <- system.file("examples/mlc", package="ggtree")
beast_tree <- read.beast(beast_file)
codeml_tree <- read.codeml(rst_file, mlc_file)

merged_tree <- merge_tree(beast_tree, codeml_tree)
merged_tree

## ----warning=FALSE, fig.width=9, fig.height=3----------------------------
library(tidytree)
library(ggplot2)

as_tibble(merged_tree) %>%
    dplyr::select(dN_vs_dS, dN, dS, rate) %>%
    subset(dN_vs_dS >=0 & dN_vs_dS <= 1.5) %>%
    tidyr::gather(type, value, dN_vs_dS:dS) %>%
    ggplot(aes(rate, value)) + geom_hex() +
    facet_wrap(~factor(type, levels = c('dN_vs_dS', 'dN', 'dS')),
               scale='free_y') +
    ylab(NULL)

## ------------------------------------------------------------------------
phylo <- as.phylo(beast_tree)
N <- Nnode2(phylo)
d <- data_frame(node = 1:N, fake_trait = rnorm(N), another_trait = runif(N))
fake_tree <- treedata(phylo = phylo, data = d)
triple_tree <- merge_tree(merged_tree, fake_tree)
triple_tree

## ------------------------------------------------------------------------
# or get.tree
as.phylo(beast_tree)

## ------------------------------------------------------------------------
get.fields(beast_tree)

## ------------------------------------------------------------------------
get.data(beast_tree)

## ------------------------------------------------------------------------
beast_tree[, c("node", "height")]
head(beast_tree[["height_median"]])

