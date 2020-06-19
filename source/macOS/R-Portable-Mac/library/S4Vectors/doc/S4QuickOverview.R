### R code from vignette source 'S4QuickOverview.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: setup
###################################################
options(width=60)
library(Matrix)
library(IRanges)
library(ShortRead)
library(graph)


###################################################
### code chunk number 2: S4_object_in_dataset
###################################################
library(graph)
data(apopGraph)
apopGraph


###################################################
### code chunk number 3: S4_object_from_constructor
###################################################
library(IRanges)
IRanges(start=c(101, 25), end=c(110, 80))


###################################################
### code chunk number 4: S4_object_from_ceorcion
###################################################
library(Matrix)
m <- matrix(3:-4, nrow=2)
as(m, "Matrix")


###################################################
### code chunk number 5: S4_object_from_high_level_IO_function
###################################################
library(ShortRead)
path_to_my_data <- system.file(
    package="ShortRead",
    "extdata", "Data", "C1-36Firecrest", "Bustard", "GERALD")
lane1 <- readFastq(path_to_my_data, pattern="s_1_sequence.txt")
lane1


###################################################
### code chunk number 6: S4_object_inside_another_object
###################################################
sread(lane1)


###################################################
### code chunk number 7: getters_and_setters
###################################################
ir <- IRanges(start=c(101, 25), end=c(110, 80))
width(ir)
width(ir) <- width(ir) - 5
ir


###################################################
### code chunk number 8: specialized_methods
###################################################
qa1 <- qa(lane1, lane="lane1")
class(qa1)


###################################################
### code chunk number 9: showMethods
###################################################
showMethods("qa")


###################################################
### code chunk number 10: showClass
###################################################
class(lane1)
showClass("ShortReadQ")


###################################################
### code chunk number 11: setClass
###################################################
setClass("SNPLocations",
    slots=c(
      genome="character",  # a single string
      snpid="character",   # a character vector of length N
      chrom="character",   # a character vector of length N
      pos="integer"        # an integer vector of length N
    )
)


###################################################
### code chunk number 12: SNPLocations
###################################################
SNPLocations <- function(genome, snpid, chrom, pos)
    new("SNPLocations", genome=genome, snpid=snpid, chrom=chrom, pos=pos)


###################################################
### code chunk number 13: test_SNPLocations
###################################################
snplocs <- SNPLocations("hg19",
             c("rs0001", "rs0002"),
             c("chr1", "chrX"),
             c(224033L, 1266886L))


###################################################
### code chunk number 14: length
###################################################
setMethod("length", "SNPLocations", function(x) length(x@snpid))


###################################################
### code chunk number 15: test_length
###################################################
length(snplocs)  # just testing


###################################################
### code chunk number 16: genome
###################################################
setGeneric("genome", function(x) standardGeneric("genome"))
setMethod("genome", "SNPLocations", function(x) x@genome)


###################################################
### code chunk number 17: snpid
###################################################
setGeneric("snpid", function(x) standardGeneric("snpid"))
setMethod("snpid", "SNPLocations", function(x) x@snpid)


###################################################
### code chunk number 18: chrom
###################################################
setGeneric("chrom", function(x) standardGeneric("chrom"))
setMethod("chrom", "SNPLocations", function(x) x@chrom)


###################################################
### code chunk number 19: pos
###################################################
setGeneric("pos", function(x) standardGeneric("pos"))
setMethod("pos", "SNPLocations", function(x) x@pos)


###################################################
### code chunk number 20: test_slot_getters
###################################################
genome(snplocs)  # just testing
snpid(snplocs)   # just testing


###################################################
### code chunk number 21: show
###################################################
setMethod("show", "SNPLocations",
    function(object)
        cat(class(object), "instance with", length(object),
            "SNPs on genome", genome(object), "\n")
)


###################################################
### code chunk number 22: S4QuickOverview.Rnw:374-375
###################################################
snplocs  # just testing


###################################################
### code chunk number 23: validity
###################################################
setValidity("SNPLocations",
    function(object) {
        if (!is.character(genome(object)) ||
            length(genome(object)) != 1 || is.na(genome(object)))
            return("'genome' slot must be a single string")
        slot_lengths <- c(length(snpid(object)),
                          length(chrom(object)),
                          length(pos(object)))
        if (length(unique(slot_lengths)) != 1)
            return("lengths of slots 'snpid', 'chrom' and 'pos' differ")
        TRUE
    }
)


###################################################
### code chunk number 24: set_chrom
###################################################
setGeneric("chrom<-", function(x, value) standardGeneric("chrom<-"))
setReplaceMethod("chrom", "SNPLocations",
    function(x, value) {x@chrom <- value; validObject(x); x})


###################################################
### code chunk number 25: test_slot_setters
###################################################
chrom(snplocs) <- LETTERS[1:2]  # repair currently broken object


###################################################
### code chunk number 26: setAs
###################################################
setAs("SNPLocations", "data.frame",
    function(from)
        data.frame(snpid=snpid(from), chrom=chrom(from), pos=pos(from))
)


###################################################
### code chunk number 27: test_coercion
###################################################
as(snplocs, "data.frame")  # testing


###################################################
### code chunk number 28: AnnotatedSNPs
###################################################
setClass("AnnotatedSNPs",
    contains="SNPLocations",
    slots=c(
        geneid="character"  # a character vector of length N
    )
)


###################################################
### code chunk number 29: slot_inheritance
###################################################
showClass("AnnotatedSNPs")


###################################################
### code chunk number 30: AnnotatedSNPs
###################################################
AnnotatedSNPs <- function(genome, snpid, chrom, pos, geneid)
{
    new("AnnotatedSNPs",
        SNPLocations(genome, snpid, chrom, pos),
        geneid=geneid)
}


###################################################
### code chunk number 31: method_inheritance
###################################################
snps <- AnnotatedSNPs("hg19",
             c("rs0001", "rs0002"),
             c("chr1", "chrX"),
             c(224033L, 1266886L),
             c("AAU1", "SXW-23"))


###################################################
### code chunk number 32: method_inheritance
###################################################
snps


###################################################
### code chunk number 33: as_data_frame_is_not_right
###################################################
as(snps, "data.frame")  # the 'geneid' slot is ignored


###################################################
### code chunk number 34: S4QuickOverview.Rnw:527-530
###################################################
is(snps, "AnnotatedSNPs")     # 'snps' is an AnnotatedSNPs object
is(snps, "SNPLocations")      # and is also a SNPLocations object
class(snps)                   # but is *not* a SNPLocations *instance*


###################################################
### code chunk number 35: automatic_coercion_method
###################################################
as(snps, "SNPLocations")


###################################################
### code chunk number 36: incremental_validity_method
###################################################
setValidity("AnnotatedSNPs",
    function(object) {
        if (length(object@geneid) != length(object))
            return("'geneid' slot must have the length of the object")
        TRUE
    }
)


