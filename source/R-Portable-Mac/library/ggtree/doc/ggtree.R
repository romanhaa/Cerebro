## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
		   message = FALSE)

## ----echo=FALSE----------------------------------------------------------
CRANpkg <- function (pkg) {
    cran <- "https://CRAN.R-project.org/package"
    fmt <- "[%s](%s=%s)"
    sprintf(fmt, pkg, cran, pkg)
}

Biocpkg <- function (pkg) {
    sprintf("[%s](http://bioconductor.org/packages/%s)", pkg, pkg)
}


## ----treeman, echo=FALSE, out.extra='', message=FALSE--------------------
treeman <- matrix(c(
  "collapse", "collapse a selecting clade",
  "expand", "expand collapsed clade",
  "flip", "exchange position of 2 clades that share a parent node",
  "groupClade", "grouping clades",
  "groupOTU", "grouping OTUs by tracing back to most recent common ancestor",
  "identify", "interactive tree manipulation",
  "rotate", "rotating a selected clade by 180 degree",
  "rotate_tree", "rotating circular layout tree by specific angle",
  "scaleClade", "zoom in or zoom out selecting clade",
  "open_tree", "convert a tree to fan layout by specific open angle"
), ncol=2, byrow=TRUE)
treeman <- as.data.frame(treeman)
colnames(treeman) <- c("Function", "Descriptiotn")
knitr::kable(treeman, caption = "Tree manipulation functions.", booktabs = T)

## ----geoms, echo=FALSE, message=FALSE------------------------------------
geoms <- matrix(c(
  "geom_balance", "highlights the two direct descendant clades of an internal node",
  "geom_cladelabel", "annotate a clade with bar and text label",
  "geom_cladelabel2", "annotate a clade with bar and text label for unrooted layout",
  "geom_hilight", "highlight a clade with rectangle",
  "geom_hilight_encircle", "highlight a clade with xspline for unrooted layout",
  "geom_label2", "modified version of geom_label, with subsetting supported",
  "geom_nodelab", "layer for node labels, which can be text or image",
  "geom_nodepoint", "annotate internal nodes with symbolic points",
  "geom_point2", "modified version of geom_point, with subsetting supported",
  "geom_range", "bar layer to present uncertainty of evolutionary inference",
  "geom_rootpoint", "annotate root node with symbolic point",
  "geom_segment2", "modified version of geom_segment, with subsetting supported",
  "geom_strip", "annotate associated taxa with bar and (optional) text label",
  "geom_taxalink", "associate two related taxa by linking them with a curve",
  "geom_text2", "modified version of geom_text, with subsetting supported",
  "geom_tiplab", "layer of tip labels, which can be text or image",
  "geom_tiplab2", "layer of tip labels for circular layout",
  "geom_tippoint", "annotate external nodes with symbolic points",
  "geom_tree", "tree structure layer, with multiple layout supported",
  "geom_treescale", "tree branch scale legend"
), ncol=2, byrow=TRUE)
geoms <- as.data.frame(geoms)
colnames(geoms) <- c("Layer", "Description")
knitr::kable(geoms, caption = "Geom layers defined in ggtree.", booktabs = T)

## ----echo=FALSE----------------------------------------------------------
sessionInfo()

