## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
		   message = FALSE)

## ----echo=FALSE, results="hide", message=FALSE---------------------------
library("ape")
library("ggplot2")
library("cowplot")
library("ggtree")
expand <- ggtree::expand
rotate <- ggtree::rotate
flip <- ggtree::flip
collapse <- ggtree:::collapse.ggtree

## ------------------------------------------------------------------------
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)
ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

## ------------------------------------------------------------------------
MRCA(tree, tip=c('A', 'E'))
MRCA(tree, tip=c('H', 'G'))

p <- ggtree(tree)
MRCA(p, tip=c('A', 'E'))

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
viewClade(p+geom_tiplab(), node=21)

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
tree <- groupClade(tree, .node=21)
ggtree(tree, aes(color=group, linetype=group))

## ----eval=FALSE----------------------------------------------------------
#  ggtree(read.tree(nwk)) %>% groupClade(.node=21) + aes(color=group, linetype=group)

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
tree <- groupClade(tree, .node=c(21, 17))
ggtree(tree, aes(color=group, linetype=group)) + geom_tiplab(aes(subset=(group==2)))

## ------------------------------------------------------------------------
tree <- groupOTU(tree, .node=c("D", "E", "F", "G"))

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
ggtree(tree, aes(color=group)) + geom_tiplab()

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
cls <- list(c1=c("A", "B", "C", "D", "E"),
            c2=c("F", "G", "H"),
            c3=c("L", "K", "I", "J"),
            c4="M")

tree <- groupOTU(tree, cls)
library("colorspace")
ggtree(tree, aes(color=group, linetype=group)) + geom_tiplab() +
     scale_color_manual(values=c("black", rainbow_hcl(4))) + theme(legend.position="right")

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
p <- ggtree(tree)
groupOTU(p, LETTERS[1:5]) + aes(color=group) + geom_tiplab() + scale_color_manual(values=c("black", "firebrick"))

## ----fig.width=6, fig.height=6-------------------------------------------
library("ape")
data(chiroptera)
groupInfo <- split(chiroptera$tip.label, gsub("_\\w+", "", chiroptera$tip.label))
chiroptera <- groupOTU(chiroptera, groupInfo)
ggtree(chiroptera, aes(color=group), layout='circular') + geom_tiplab(size=1, aes(angle=angle))

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
cp <- collapse(p, node=21)
cp + geom_point2(aes(subset=(node == 21)), size=5, shape=23, fill="steelblue")

## ----fig.width=5, fig.height=5, fig.align="center", warning=FALSE--------
cp %>% expand(node=21)

## ----fig.width=15, fig.height=3, warning=FALSE---------------------------
p1 <- ggtree(tree)
p2 <- collapse(p1, 21) + geom_point2(aes(subset=(node==21)), size=5, shape=23, fill="blue")
p3 <- collapse(p2, 17) + geom_point2(aes(subset=(node==17)), size=5, shape=23, fill="red")
p4 <- expand(p3, 17)
p5 <- expand(p4, 21)

library(cowplot)
plot_grid(p1, p2, p3, p4, p5, ncol=5)

## ----fig.width=12, fig.height=6, warning=F-------------------------------
plot_grid(ggtree(tree) + geom_hilight(21, "steelblue"),
          ggtree(tree) %>% scaleClade(21, scale=0.3) + geom_hilight(21, "steelblue"),
          ncol=2)

## ----fig.width=12, fig.height=6, warning=F-------------------------------
plot_grid(ggtree(tree) + geom_hilight(17, fill="steelblue") +
                 geom_hilight(21, fill="darkgreen"),
          ggtree(tree) %>% scaleClade(17, scale=2) %>% scaleClade(21, scale=0.3) +
                 geom_hilight(17, "steelblue") + geom_hilight(21, fill="darkgreen"),
          ncol=2)

## ----fig.width=8, fig.height=4-------------------------------------------
tree <- groupClade(tree, c(21, 17))
p <- ggtree(tree, aes(color=group)) + scale_color_manual(values=c("black", "firebrick", "steelblue"))
p2 <- rotate(p, 21) %>% rotate(17)
plot_grid(p, p2, ncol=2)

## ----eval=F--------------------------------------------------------------
#  set.seed(2016-05-29)
#  p <- ggtree(tree <- rtree(50)) + geom_tiplab()
#  for (n in reorder(tree, 'postorder')$edge[,1] %>% unique) {
#      p <- rotate(p, n)
#      print(p + geom_point2(aes(subset=(node == n)), color='red'))
#  }

## ----fig.width=8, fig.height=4-------------------------------------------
plot_grid(p, flip(p, 17, 21), ncol=2)

## ----eval=FALSE----------------------------------------------------------
#  set.seed(123)
#  tr <- rtree(50)
#  
#  p <- ggtree(tr, layout='circular') + geom_tiplab2()
#  
#  for (angle in seq(0, 270, 10)) {
#      print(open_tree(p, angle=angle) + ggtitle(paste("open angle:", angle)))
#  }

## ----eval=FALSE----------------------------------------------------------
#  for (angle in seq(0, 270, 10)) {
#      print(rotate_tree(p, angle) + ggtitle(paste("rotate angle:", angle)))
#  }

