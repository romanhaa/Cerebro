## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
		   message = FALSE)

## ----echo=FALSE, results="hide", message=FALSE---------------------------
library("ape")
library("grid")
library("ggplot2")
library("cowplot")
library("treeio")
library("ggtree")


CRANpkg <- function (pkg) {
    cran <- "https://CRAN.R-project.org/package"
    fmt <- "[%s](%s=%s)"
    sprintf(fmt, pkg, cran, pkg)
}

Biocpkg <- function (pkg) {
    sprintf("[%s](http://bioconductor.org/packages/%s)", pkg, pkg)
}

## ------------------------------------------------------------------------
library("treeio")
library("ggtree")

nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
ggtree(tree, color="firebrick", size=1, linetype="dotted")

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
ggtree(tree, ladderize=FALSE)

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
ggtree(tree, branch.length="none")

## ----eval=F--------------------------------------------------------------
#  library(ggtree)
#  set.seed(2017-02-16)
#  tr <- rtree(50)
#  ggtree(tr)
#  ggtree(tr, layout="slanted")
#  ggtree(tr, layout="circular")
#  ggtree(tr, layout="fan", open.angle=120)
#  ggtree(tr, layout="equal_angle")
#  ggtree(tr, layout="daylight")
#  ggtree(tr, branch.length='none')
#  ggtree(tr, branch.length='none', layout='circular')
#  ggtree(tr, layout="daylight", branch.length = 'none')

## ----echo=F,  fig.width=8, fig.height=8, message=FALSE-------------------
library(ggtree)
set.seed(2017-02-16)
tr <- rtree(50)
library(cowplot)
theme_layout <- theme(plot.title = element_text(hjust = 0.5))
plot_grid(
    ggtree(tr) + ggtitle("rectangular (phylogram)")+ theme_layout,
    ggtree(tr, layout="slanted") + ggtitle("slanted (phylogram)")+theme_layout,
    ggtree(tr, layout="circular") + ggtitle("circular (phylogram)")+theme_layout,
    ggtree(tr, layout="fan", open.angle=120) + ggtitle("fan (phylogram)")+theme_layout,
    ggtree(tr, layout="equal_angle")+ ggtitle("equal angle (unrooted)")+theme_layout,
    ggtree(tr, layout="daylight")+ ggtitle("daylight (unrooted)")+theme_layout,
    ggtree(tr, branch.length='none')+ ggtitle("rectangular (cladogram)")+theme_layout,
    ggtree(tr, branch.length='none', layout='circular')+ ggtitle("circular (cladogram)")+theme_layout,
    ggtree(tr, layout="daylight", branch.length = 'none')+ ggtitle("daylight (cladogram)")+theme_layout,
    ncol=3)

## ----eval=FALSE----------------------------------------------------------
#  ggtree(tr) + scale_x_reverse()
#  ggtree(tr) + coord_flip()
#  ggtree(tr) + scale_x_reverse() + coord_flip()
#  print(ggtree(tr), newpage=TRUE, vp = grid::viewport(angle=-30, width=.9, height=.9))
#  ggtree(tr, layout='slanted') + coord_flip()
#  ggtree(tr, layout='slanted', branch.length='none') +
#      coord_flip() + scale_y_reverse() +scale_x_reverse()
#  ggtree(tr, layout='circular') + xlim(-10, NA)
#  ggtree(tr) + scale_x_reverse() + coord_polar(theta='y')
#  ggtree(tr) + scale_x_reverse(limits=c(10, 0)) + coord_polar(theta='y')

## ----fig.keep='none', echo=FALSE, warning=FALSE--------------------------
tree_angle <- grid::grid.grabExpr(print(ggtree(tr), newpage=TRUE, vp = grid::viewport(angle=-30, width=.9, height=.9)))

## ----fig.width=8, fig.height = 8, echo=FALSE, warning=FALSE--------------
plot_grid(
    ggtree(tr) + scale_x_reverse(),
    ggtree(tr) + coord_flip(),
    ggtree(tr) + scale_x_reverse() + coord_flip(),
    tree_angle,
    ggtree(tr, layout='slanted') + coord_flip(),
    ggtree(tr, layout='slanted', branch.length='none') + coord_flip() + scale_y_reverse() +scale_x_reverse(),
    ggtree(tr, layout='circular') + xlim(-10, NA),
    ggtree(tr) + scale_x_reverse() + coord_polar(theta='y'),
    ggtree(tr) + scale_x_reverse(limits=c(15, 0)) + coord_polar(theta='y'),
    ncol = 3, labels = LETTERS[1:9])

## ----fig.width=8, fig.height=4, fig.align="center"-----------------------
tree2d <- read.beast(system.file("extdata", "twoD.tree", package="treeio"))
ggtree(tree2d, mrsd = "2014-05-01") + theme_tree2()

## ----fig.width=9, fig.height=4, fig.align="center"-----------------------
ggtree(tree2d, mrsd = "2014-05-01",
       yscale="NGS", yscale_mapping=c(N2=2, N3=3, N4=4, N5=5, N6=6, N7=7)) +
           theme_classic() + theme(axis.line.x=element_line(), axis.line.y=element_line()) +
               theme(panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3),
                     panel.grid.major.y=element_blank()) +
                         scale_y_continuous(labels=paste0("N", 2:7))

## ----fig.width=4, fig.height=4, fig.align="center"-----------------------
ggtree(tree) + geom_treescale()

## ----eval=F--------------------------------------------------------------
#  ggtree(tree) + geom_treescale(x=0, y=12, width=6, color='red')
#  ggtree(tree) + geom_treescale(fontsize=8, linesize=2, offset=-1)

## ----fig.width=8, fig.height=4, fig.align="center", echo=F---------------
plot_grid(
    ggtree(tree)+geom_treescale(x=0, y=12, width=6, color='red'),
    ggtree(tree)+geom_treescale(fontsize=8, linesize=2, offset=-1),
    ncol = 2, labels = LETTERS[1:2])

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
ggtree(tree) + theme_tree2()

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
ggtree(tree)+geom_point(aes(shape=isTip, color=isTip), size=3)

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
p <- ggtree(tree) + geom_nodepoint(color="#b5e521", alpha=1/4, size=10)
p + geom_tippoint(color="#FDAC4F", shape=8, size=3)

## ----fig.width=3, fig.height=3, warning=FALSE, fig.align="center"--------
p + geom_tiplab(size=3, color="purple")

## ----fig.width=6, fig.height=6, warning=FALSE, fig.align="center"--------
ggtree(tree, layout="circular") + geom_tiplab(aes(angle=angle), color='blue')

## ----fig.width=6, fig.height=6, warning=FALSE, fig.align="center"--------
ggtree(tree, layout="circular") + geom_tiplab2(color='blue')

## ----fig.width=4, fig.height=3, warning=FALSE, fig.align="center"--------
p + geom_tiplab(aes(x=branch), size=3, color="purple", vjust=-0.3)

## ----fig.width=3, fig.height=3, fig.align="center"-----------------------
p %<% rtree(50)

## ----eval=F--------------------------------------------------------------
#  ggtree(rtree(30), color="red") + theme_tree("steelblue")
#  ggtree(rtree(20), color="white") + theme_tree("black")

## ----fig.width=8, fig.height=3, fig.align="center", echo=F---------------
cowplot::plot_grid(
    ggtree(rtree(30), color="red") + theme_tree("steelblue"),
    ggtree(rtree(20), color="purple") + theme_tree("black"),
    ncol=2)

## ----fig.width=12, fig.height=4------------------------------------------
trees <- lapply(c(10, 20, 40), rtree)
class(trees) <- "multiPhylo"
ggtree(trees) + facet_wrap(~.id, scale="free") + geom_tiplab()

## ----fig.width=20, fig.height=20-----------------------------------------
btrees <- read.tree(system.file("extdata/RAxML", "RAxML_bootstrap.H3", package="treeio"))
ggtree(btrees) + facet_wrap(~.id, ncol=10)

## ------------------------------------------------------------------------
p <- ggtree(btrees, layout="rectangular",   color="lightblue", alpha=.3)

best_tree <- read.tree(system.file("extdata/RAxML", "RAxML_bipartitionsBranchLabels.H3", package="treeio"))
df <- fortify(best_tree, branch.length = 'none')
p+geom_tree(data=df, color='firebrick')

## ----fig.width=10, fig.height=5------------------------------------------
library("treeio")
beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)
beast_tree
p1 <- ggtree(beast_tree, mrsd='2013-01-01') + theme_tree2() +
    ggtitle("Divergence time")
p2 <- ggtree(beast_tree, branch.length = 'rate') + theme_tree2() +
    ggtitle("Substitution rate")

library(cowplot)
plot_grid(p1, p2, ncol=2)

## ----fig.width=10, fig.height=5------------------------------------------
mlcfile <- system.file("extdata/PAML_Codeml", "mlc", package="treeio")
mlc_tree <- read.codeml_mlc(mlcfile)
p1 <- ggtree(mlc_tree) + theme_tree2() +
    ggtitle("nucleotide substitutions per codon")
p2 <- ggtree(mlc_tree, branch.length='dN_vs_dS') + theme_tree2() +
    ggtitle("dN/dS tree")
plot_grid(p1, p2, ncol=2)

## ------------------------------------------------------------------------
beast_tree2 <- rescale_tree(beast_tree, branch.length = 'rate')
ggtree(beast_tree2) + theme_tree2()

## ----fig.width=9, fig.height=5, fig.align="center"-----------------------
library("ape")
data(chiroptera)
library("ggtree")
gzoom(chiroptera, grep("Plecotus", chiroptera$tip.label))

## ----fig.width=9, fig.height=5, message=FALSE, warning=FALSE-------------
groupInfo <- split(chiroptera$tip.label, gsub("_\\w+", "", chiroptera$tip.label))
chiroptera <- groupOTU(chiroptera, groupInfo)
p <- ggtree(chiroptera, aes(color=group)) + geom_tiplab() + xlim(NA, 23)
gzoom(p, grep("Plecotus", chiroptera$tip.label), xmax_adjust=2)

## ----fig.width=5, fig.height=5-------------------------------------------
ggtree(beast_tree, aes(color=rate)) +
    scale_color_continuous(low='darkgreen', high='red') +
    theme(legend.position="right")

