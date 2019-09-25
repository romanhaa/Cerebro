### R code from vignette source 'GOusage.Rnw'

###################################################
### code chunk number 1: Setup
###################################################
library("Biobase")
library("annotate")
library("xtable")
require("Rgraphviz", quietly=TRUE)
library("hgu95av2.db")
library("GO.db")


###################################################
### code chunk number 2: parentrel
###################################################
 GOTERM$"GO:0003700"

 GOMFPARENTS$"GO:0003700"
 GOMFCHILDREN$"GO:0003700"



###################################################
### code chunk number 3: locusid
###################################################

 ll1 = hgu95av2GO[["39613_at"]]
 length(ll1)
 sapply(ll1, function(x) x$Ontology)



###################################################
### code chunk number 4: getmappings
###################################################

getOntology(ll1, "BP")
getEvidence(ll1)
zz = dropECode(ll1)
getEvidence(zz)



###################################################
### code chunk number 5: sizeofonts
###################################################

 zz = Ontology(GOTERM)
 table(unlist(zz))



###################################################
### code chunk number 6: isa-partof
###################################################

 BPisa = eapply(GOBPPARENTS, function(x) names(x))
 table(unlist(BPisa))

 MFisa = eapply(GOMFPARENTS, function(x) names(x))
 table(unlist(MFisa))

 CCisa = eapply(GOCCPARENTS, function(x) names(x))
 table(unlist(CCisa))



###################################################
### code chunk number 7: finding these
###################################################
 goterms = unlist(Term(GOTERM))
 whmf = grep("fertilization", goterms)


###################################################
### code chunk number 8: subsetGT
###################################################
 goterms[whmf]



###################################################
### code chunk number 9: getMF
###################################################
affyGO = eapply(hgu95av2GO, getOntology)
table(sapply(affyGO, length))



###################################################
### code chunk number 10: getEvidence
###################################################
affyEv = eapply(hgu95av2GO, getEvidence)

table(unlist(affyEv, use.names=FALSE))



###################################################
### code chunk number 11: dropOneEvidence
###################################################
test1 = eapply(hgu95av2GO, dropECode, c("IEA", "NR"))

table(unlist(sapply(test1, getEvidence),
             use.names=FALSE))


