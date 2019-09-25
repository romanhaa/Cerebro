### R code from vignette source 'useProbeInfo.Rnw'

###################################################
### code chunk number 1: loadlibs
###################################################
library("annotate")
library("rae230a.db")
library("rae230aprobe")


###################################################
### code chunk number 2: selprobe
###################################################

ps = names(as.list(rae230aACCNUM))

myp = ps[1001]

myA = get(myp, rae230aACCNUM)

wp = rae230aprobe$Probe.Set.Name == myp
myPr = rae230aprobe[wp,]



###################################################
### code chunk number 3: getACC
###################################################

myseq = getSEQ(myA)
nchar(myseq)

library("Biostrings")
mybs = DNAString(myseq)

match1 = matchPattern(as.character(myPr[1,1]), mybs)
match1
as.matrix(ranges(match1))
myPr[1,5]


###################################################
### code chunk number 4: getRev
###################################################

myp = ps[100]

myA = get(myp, rae230aACCNUM)

wp = rae230aprobe$Probe.Set.Name == myp

myPr = rae230aprobe[wp,]

myseq = getSEQ(myA)

mybs = DNAString(myseq)

Prstr = as.character(myPr[1,1])

match2 = matchPattern(Prstr, mybs)

## expecting 0 (no match)
length(match2)

match2 = matchPattern(reverseComplement(DNAString(Prstr)), mybs)

nchar(match2)

nchar(myseq) - as.matrix(ranges(match2))
myPr[1,5]


###################################################
### code chunk number 5: useProbeInfo.Rnw:159-160
###################################################
sessionInfo()


