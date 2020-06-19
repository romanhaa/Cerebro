## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
		   message = FALSE)

## ----echo=FALSE, results="hide", message=FALSE---------------------------
library("ape")
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

inset <- ggtree::inset

## ------------------------------------------------------------------------
set.seed(2015-12-21)
tree <- rtree(30)
p <- ggtree(tree) + xlim(NA, 6)

p + geom_cladelabel(node=45, label="test label") +
    geom_cladelabel(node=34, label="another clade")

## ------------------------------------------------------------------------
p + geom_cladelabel(node=45, label="test label", align=TRUE, offset=.5) +
    geom_cladelabel(node=34, label="another clade", align=TRUE, offset=.5)

## ------------------------------------------------------------------------
p + geom_cladelabel(node=45, label="test label", align=T, color='red') +
    geom_cladelabel(node=34, label="another clade", align=T, color='blue')

## ------------------------------------------------------------------------
p + geom_cladelabel(node=45, label="test label", align=T, angle=270, hjust='center', offset.text=.5) +
    geom_cladelabel(node=34, label="another clade", align=T, angle=45)

## ------------------------------------------------------------------------
p + geom_cladelabel(node=45, label="test label", align=T, angle=270, hjust='center', offset.text=.5, barsize=1.5) +
    geom_cladelabel(node=34, label="another clade", align=T, angle=45, fontsize=8)

## ------------------------------------------------------------------------
p + geom_cladelabel(node=34, label="another clade", align=T, geom='label', fill='lightblue')

## ----fig.wdith=7, fig.height=7, fig.align='center', warning=FALSE, message=FALSE----
pg <- ggtree(tree, layout="daylight")
pg + geom_cladelabel2(node=45, label="test label", angle=10) +
    geom_cladelabel2(node=34, label="another clade", angle=305)

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)
ggtree(tree) + geom_tiplab() + 
  geom_strip(5, 7, barsize=2, color='red') + 
  geom_strip(6, 12, barsize=2, color='blue')

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
ggtree(tree) + geom_hilight(node=21, fill="steelblue", alpha=.6) +
    geom_hilight(node=17, fill="darkgreen", alpha=.6)

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
ggtree(tree, layout="circular") + geom_hilight(node=21, fill="steelblue", alpha=.6) +
    geom_hilight(node=23, fill="darkgreen", alpha=.6)

## ----fig.width=4, fig.height=5, fig.align='center', warning=FALSE--------
ggtree(tree) +
  geom_balance(node=16, fill='steelblue', color='white', alpha=0.6, extend=1) +
  geom_balance(node=19, fill='darkgreen', color='white', alpha=0.6, extend=1)

## ----fig.width=5, fig.height=5, fig.align='center', warning=FALSE, message=FALSE----
pg + geom_hilight_encircle(node=45) + geom_hilight_encircle(node=34, fill='darkgreen')

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
ggtree(tree) + geom_tiplab() + geom_taxalink('A', 'E') + 
  geom_taxalink('F', 'K', color='red', arrow=grid::arrow(length=grid::unit(0.02, "npc")))

## ----warning=FALSE, fig.width=5, fig.height=5, fig.align='center'--------
file <- system.file("extdata/BEAST", "beast_mcc.tree", package="treeio")
beast <- read.beast(file)
ggtree(beast, aes(color=rate))  +
    geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
    geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
    scale_color_continuous(low="darkgreen", high="red") +
    theme(legend.position=c(.1, .8))

## ------------------------------------------------------------------------
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)
p <- ggtree(tree)

dd <- data.frame(taxa = LETTERS[1:13],
                 place = c(rep("GZ", 5), rep("HK", 3), rep("CZ", 4), NA),
                 value = round(abs(rnorm(13, mean=70, sd=10)), digits=1))
## you don't need to order the data
## data was reshuffled just for demonstration
dd <- dd[sample(1:13, 13), ]
row.names(dd) <- NULL

## ----eval=FALSE----------------------------------------------------------
#  print(dd)

## ----echo=FALSE, results='asis'------------------------------------------
knitr::kable(dd)

## ----fig.width=6, fig.height=5, warning=FALSE, fig.align="center"--------
p <- p %<+% dd + geom_tiplab(aes(color=place)) +
       geom_tippoint(aes(size=value, shape=place, color=place), alpha=0.25)
p + theme(legend.position="right")

## ----fig.width=6, fig.height=5, warning=FALSE, fig.align="center"--------
p + geom_text(aes(color=place, label=place), hjust=1, vjust=-0.4, size=3) +
    geom_text(aes(color=place, label=value), hjust=1, vjust=1.4, size=3)

## ----fig.width=8, fig.height=6, fig.align="center", warning=FALSE, message=FALSE----
beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)

genotype_file <- system.file("examples/Genotype.txt", package="ggtree")
genotype <- read.table(genotype_file, sep="\t", stringsAsFactor=F)
colnames(genotype) <- sub("\\.$", "", colnames(genotype))
p <- ggtree(beast_tree, mrsd="2013-01-01") + geom_treescale(x=2008, y=1, offset=2)
p <- p + geom_tiplab(size=2)
gheatmap(p, genotype, offset=5, width=0.5, font.size=3, colnames_angle=-45, hjust=0) +
    scale_fill_manual(breaks=c("HuH3N2", "pdm", "trig"), values=c("steelblue", "firebrick", "darkgreen"))

## ----fig.width=8, fig.height=6, fig.align="center", warning=FALSE--------
p <- ggtree(beast_tree, mrsd="2013-01-01") + geom_tiplab(size=2, align=TRUE, linesize=.5) + theme_tree2()
pp <- (p + scale_y_continuous(expand=c(0, 0.3))) %>%
    gheatmap(genotype, offset=8, width=0.6, colnames=FALSE) %>%
        scale_x_ggtree()
pp + theme(legend.position="right")

## ----fig.width=8, fig.height=6, fig.align='center', warning=FALSE--------
fasta <- system.file("examples/FluA_H3_AA.fas", package="ggtree")
msaplot(ggtree(beast_tree), fasta)

## ----fig.width=7, fig.height=7, fig.align='center', warning=FALSE--------
msaplot(ggtree(beast_tree), fasta, window=c(150, 200)) + coord_polar(theta='y')

## ----warning=F, fig.width=10, fig.height=6-------------------------------
tr <- rtree(30)

d1 <- data.frame(id=tr$tip.label, val=rnorm(30, sd=3))
p <- ggtree(tr)

p2 <- facet_plot(p, panel="dot", data=d1, geom=geom_point, aes(x=val), color='firebrick')
d2 <- data.frame(id=tr$tip.label, value=abs(rnorm(30, mean=100, sd=50)))

facet_plot(p2, panel='bar', data=d2, geom=geom_segment, aes(x=0, xend=value, y=y, yend=y), size=3, color='steelblue') + theme_tree2()

