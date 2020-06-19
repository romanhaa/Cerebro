### R code from vignette source 'annotate.Rnw'

###################################################
### code chunk number 1: loadlibs
###################################################
library("annotate")
library("hgu95av2.db")
ls("package:hgu95av2.db")


###################################################
### code chunk number 2: qc
###################################################
hgu95av2()


###################################################
### code chunk number 3: locusid
###################################################
hgu95av2ENTREZID


###################################################
### code chunk number 4: getting
###################################################

get("1000_at", env=hgu95av2ENTREZID)
hgu95av2ENTREZID[["1000_at"]]
hgu95av2ENTREZID$"1000_at"



###################################################
### code chunk number 5: tolist
###################################################
LLs = as.list(hgu95av2ENTREZID)
length(LLs)
names(LLs)[1:10]



###################################################
### code chunk number 6: annotate.Rnw:216-217
###################################################
sessionInfo()


