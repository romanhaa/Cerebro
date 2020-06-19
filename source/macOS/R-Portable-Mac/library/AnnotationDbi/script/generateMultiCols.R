## This script is used to pre-compute information about which cols()
## values have a many to one relationship with the _id cols in each
## table, and how bad that relationship is.  The values are scanned
## across all the org packages and then collected into a vector so
## that they can be used by AnntoationDbi to make sensible
## recommendations about what cols() a user should ask for at once, or
## if multiple cols should be requested together.

## The final output is a named integer vector.  The names of this
## vector are the values for cols().  And the values of this vector is
## a number that represents the max value ever seen for any given
## field, and the names are the names used by cols().  So this number
## is meant to be an indicator for how bad something can be when
## cols() is used.

## If cols() indicates potential danger, then the plan is to message()
## the user (or use warning("",immediately=TRUE)) and tell them that
## what they are doing is going to maybe take a long time.  Herve has
## suggested that I could just use count() at the front of the same
## query and just pre-call that to get a message that tells them how
## bad it is (so that they can escape if they wish).  Paul asked me to
## add a parameter (in that case) so that users who know what they are
## doing can switch that safety off.  Martin mentioned that being able
## to hit Ctrl-C is probably sufficient.  I agree with Martin in this
## case.




##############################################################################
##  ACK Combinatorics!
##############################################################################
## As for the actual problem.  I am just going to put together a
## hand-made blacklist of cols that I know have many to one
## relationships based on this script here.  This script will just
## learn which fields do that (not try to rate them).


## Then if the user uses more than a few of them, I will issue a
## warning as described above...  The test can trap this warning maybe
## to make sure it goes out?



## So lets just do it for HUMAN
## library(org.Hs.eg.db)
## x <- org.Hs.eg.db


getManyToOneStatus <- function(x){
    cols <- cols(x)
    testManyToOne <- function(c, x){
        k <- keys(x,"ENTREZID")
        res <- select(x, cols=c, keys=k, keytype="ENTREZID")
        if(length(unique(res[["ENTREZID"]])) <
           length(res[["ENTREZID"]])){
            return(TRUE)
        }else{
            return(FALSE)
        }
    }
    res <- unlist(lapply(cols, testManyToOne, x))
    names(res) <- cols
    res
}

## Tests
## This should be TRUE
## testManyToOne("PFAM", x)
## This should be FALSE:
## testManyToOne("SYMBOL", x)

## res <- getManyToOneStatus(x)


## so then do the above for all the org packages.
require("org.Hs.eg.db")
require("org.Mm.eg.db")
require("org.At.tair.db")
require("org.Bt.eg.db")
require("org.Cf.eg.db")
require("org.Gg.eg.db")
require("org.Dm.eg.db")
require("org.Rn.eg.db")
require("org.Ce.eg.db")
require("org.Xl.eg.db")
require("org.Sc.sgd.db")
require("org.Ss.eg.db")
require("org.Dr.eg.db")
require("org.EcK12.eg.db")
require("org.EcSakai.eg.db")


pkgs <- c(org.Hs.eg.db,
          org.Mm.eg.db,
          ##           org.At.tair.db, ## exclude b/c tair is (atypical)
          org.Bt.eg.db,
          org.Cf.eg.db,
          org.Gg.eg.db,
          org.Dm.eg.db,
          org.Rn.eg.db,
          org.Ce.eg.db,
          org.Xl.eg.db,
          ##        org.Sc.sgd.db,   ## There is a problem with this one.
          org.Ss.eg.db,
          org.Dr.eg.db,
          org.EcK12.eg.db,
          org.EcSakai.eg.db)

res <- lapply(pkgs, getManyToOneStatus)
many2Ones = res
save(many2Ones, file="many2Ones.Rda")


## then combine all these vectors into one vector.  We want to have
## the vectors added together and then call unique, but we ALSO want
## to give preference to things being TRUE.  So if you are true for
## one, you are always true from then on...

## So we just do this to sort so that the trues are in front (and will
## get grabbed by the final filtering)
blackList <- sort(unlist(res), decreasing=TRUE)
## And then we call unique (which grabs the ones it sees 1st)
blackList <- blackList[unique(names(blackList))]
## Then we just keep the ones that are marked as TRUE.
blackList <- names(blackList[blackList])
## and save it...
save(blackList, file="manyToOneBlackList.Rda")






## TROUBLESHOOTING THE strange yeast issue:
### problem with yeast is just this one:
## x <- org.Sc.sgd.db
## k <- keys(x,"ENTREZID")

## debug(AnnotationDbi:::.select)
## debug(AnnotationDbi:::.extractData)
## res <- select(x, cols="ORF", keys=k, keytype="ENTREZID")

## gives us this:
## SELECT  genes.gene_id,gene2systematic.systematic_name,sgd.sgd_id  FROM genes LEFT JOIN  gene2systematic USING ( systematic_name ) LEFT JOIN  sgd USING ( systematic_name );
## plus a huge where clause like  WHERE  genes.gene_id in ( '9164990','9

## and this here 
## res <- select(x, cols="COMMON", keys=k, keytype="ENTREZID")
## give us this:
## SELECT  genes.gene_id,gene2systematic.gene_name,sgd.sgd_id  FROM genes LEFT JOIN  gene2systematic USING ( systematic_name ) LEFT JOIN  sgd USING ( systematic_name );

## which has the same problem (a bad join)


## OK, so I have a couple options:
## 1) add _id to gene2systematic and then simplify code

## BUT problem with adding _id is just that some of the rows will not
## have one?  So: does that matter?  would those rows EVER be joined
## in any meaningful way?  It turns out that only in the context of
## keys() is that information worth anything.  IOW you can't link it
## to anything, because the sgd table doesn't have rows where there
## isn't an _id.  Is that true?  - NO! Look at this (for example):

## select * from gene2systematic LEFT JOIN sgd USING (systematic_name) WHERE systematic_name IN ('AIP5','AGS1');

## output shows that AIP5 has a _id and an sgd_id, so these exist for
## many of the systematic_name values even if they are not represented
## by a gene_name...

## AND actually: ALL of the systematic_name vals are in both sgd and
## gene2systematic.  What are different are the gene_name vals.  There
## are 10960 distinct gene names in gene2systematic but only 5872 of
## those are associated in sgd with an _id or a sgd id etc.

## So really it's just the gene_name field that can be a problem.  But
## that data is only useful in the sense of "keys", since it can only
## connect if it is associated with one of the 8699 systematic_name
## fields (or one of the _id values)



## lets test this:
## add the _id col:  It should be an integer and it should be NULL by default
## ALTER TABLE gene2systematic ADD COLUMN _id INTEGER NULL;
## then add the index on this col (it can't be a primary key for this table).

## CREATE INDEX g2s_id ON gene2systematic (_id);
## then an insert:

## And actually I am going to really want an insert that looks more like this:
## INSERT INTO gene2systematic SELECT sgd._id, g2s.gene_name, g2s.systematic_name  FROM gene2systematic as g2s LEFT JOIN sgd USING (systematic_name) ;






## OR (ruled this out already)

## and for the sake of completeness, lets work out #2 1st...
## 2) figure out why sgd has fewer gene_name values than
## gene2systematic (and possibly fix that / adjust code to just use
## sgd table instead of gene2systematic)
## It's because table sgd has a not null constraint on sgd.  So rows
## that have a gene_name and no sgd_id are excluded.  So I can't do
## this approach (or shouldn't)


