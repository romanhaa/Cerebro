### R code from vignette source 'byChroms.Rnw'

###################################################
### code chunk number 1: loaddata
###################################################

 library("annotate")
 library("hu6800.db")
 lens <- unlist(eapply(hu6800CHR, length))

 table(lens)
 wh2 = mget(names(lens)[lens==2], env = hu6800CHR)

 wh2[1]


###################################################
### code chunk number 2: fixdata
###################################################
chrs2 <- unlist(eapply(hu6800CHR, function(x) x[1]))
chrs2 <- factor(chrs2)
length(chrs2)
 table(unlist(chrs2))


###################################################
### code chunk number 3: strandloc
###################################################

 strand <- as.list(hu6800CHRLOC)

 splits <- split(strand, chrs2)
 length(splits)
 names(splits)



###################################################
### code chunk number 4: chrloc
###################################################

 newChrClass <- buildChromLocation("hu6800")



###################################################
### code chunk number 5: cPlot
###################################################

  library(geneplotter)

  cPlot(newChrClass)



