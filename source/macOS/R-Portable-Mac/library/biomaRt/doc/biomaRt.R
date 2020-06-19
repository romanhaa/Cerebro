## ----setup, cache = F, echo = FALSE----------------------------------------
knitr::opts_chunk$set(error = TRUE)

## ----annotate,echo=FALSE----------------------------------------------------------------------------------------------
## library("annotate")
options(width=120)

## ----biomaRt----------------------------------------------------------------------------------------------------------
library("biomaRt")
listMarts()

## ----putenv, eval = FALSE---------------------------------------------------------------------------------------------
#  Sys.setenv("http_proxy" = "http://my.proxy.org:9999")

## ----rCurlOptions, eval = FALSE---------------------------------------------------------------------------------------
#  options(RCurlOptions = list(proxy="uscache.kcc.com:80",proxyuserpwd="------:-------"))

## ----ensembl1---------------------------------------------------------------------------------------------------------
ensembl=useMart("ensembl")

## ----listDatasets-----------------------------------------------------------------------------------------------------
listDatasets(ensembl)

## ----ensembl2, eval=TRUE----------------------------------------------------------------------------------------------
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

## ----ensembl3---------------------------------------------------------------------------------------------------------
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

## ----filters----------------------------------------------------------------------------------------------------------
filters = listFilters(ensembl)
filters[1:5,]

## ----attributes-------------------------------------------------------------------------------------------------------
attributes = listAttributes(ensembl)
attributes[1:5,]

## ----getBM1, echo=TRUE,eval=TRUE--------------------------------------------------------------------------------------
affyids=c("202763_at","209310_s_at","207500_at")
getBM(attributes=c('affy_hg_u133_plus_2', 'entrezgene'), 
      filters = 'affy_hg_u133_plus_2', 
      values = affyids, 
      mart = ensembl)

## ----searchDatasets, echo = TRUE, eval = TRUE-------------------------------------------------------------------------
searchDatasets(mart = ensembl, pattern = "hsapiens")

## ----searchAttributes, echo = TRUE, eval = TRUE-----------------------------------------------------------------------
searchAttributes(mart = ensembl, pattern = "hgnc")

## ----searchFilters, echo = TRUE, eval = TRUE--------------------------------------------------------------------------
searchFilters(mart = ensembl, pattern = "ensembl.*id")

## ----task1, echo=TRUE,eval=TRUE---------------------------------------------------------------------------------------
affyids=c("202763_at","209310_s_at","207500_at")
getBM(attributes = c('affy_hg_u133_plus_2', 'hgnc_symbol', 'chromosome_name',
                   'start_position', 'end_position', 'band'),
      filters = 'affy_hg_u133_plus_2', 
      values = affyids, 
      mart = ensembl)

## ----task2, echo=TRUE,eval=TRUE---------------------------------------------------------------------------------------
entrez=c("673","837")
goids = getBM(attributes = c('entrezgene', 'go_id'), 
              filters = 'entrezgene', 
              values = entrez, 
              mart = ensembl)
head(goids)

## ----task3, echo=TRUE,eval=TRUE---------------------------------------------------------------------------------------
 go=c("GO:0051330","GO:0000080","GO:0000114","GO:0000082")
 chrom=c(17,20,"Y")
 getBM(attributes= "hgnc_symbol",
        filters=c("go","chromosome_name"),
        values=list(go, chrom), mart=ensembl)

## ----task4, echo=TRUE,eval=TRUE---------------------------------------------------------------------------------------
refseqids = c("NM_005359","NM_000546")
ipro = getBM(attributes=c("refseq_mrna","interpro","interpro_description"), 
             filters="refseq_mrna",
             values=refseqids, 
             mart=ensembl)
ipro

## ----task5, eval = TRUE-----------------------------------------------------------------------------------------------
getBM(attributes = c('affy_hg_u133_plus_2','ensembl_gene_id'), 
      filters = c('chromosome_name','start','end'),
      values = list(16,1100000,1250000), 
      mart = ensembl)

## ----task6, echo=TRUE, eval = TRUE------------------------------------------------------------------------------------
getBM(attributes = c('entrezgene','hgnc_symbol'), 
      filters = 'go', 
      values = 'GO:0004707', 
      mart = ensembl)

## ----task7, eval=TRUE-------------------------------------------------------------------------------------------------
entrez=c("673","7157","837")
getSequence(id = entrez, 
            type="entrezgene",
            seqType="coding_gene_flank",
            upstream=100, 
            mart=ensembl) 

## ----task8, echo=TRUE,eval=TRUE---------------------------------------------------------------------------------------
utr5 = getSequence(chromosome=3, start=185514033, end=185535839,
                   type="entrezgene",
                   seqType="5utr", 
                   mart=ensembl)
utr5

## ----task9, echo=TRUE, eval=TRUE--------------------------------------------------------------------------------------
protein = getSequence(id=c(100, 5728),
                      type="entrezgene",
                      seqType="peptide", 
                      mart=ensembl)
protein

## ----task10, echo=TRUE, eval=TRUE-------------------------------------------------------------------------------------
snpmart = useMart(biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp")

## ----task10b----------------------------------------------------------------------------------------------------------
getBM(attributes = c('refsnp_id','allele','chrom_start','chrom_strand'), 
      filters = c('chr_name','start','end'), 
      values = list(8,148350,148612), 
      mart = snpmart)

## ----getLDS-----------------------------------------------------------------------------------------------------------
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
getLDS(attributes = c("hgnc_symbol","chromosome_name", "start_position"),
       filters = "hgnc_symbol", values = "TP53",mart = human,
      attributesL = c("refseq_mrna","chromosome_name","start_position"), martL = mouse)

## ----archiveMarts1----------------------------------------------------------------------------------------------------
listMarts(archive = TRUE)

## ----archiveMarts2, echo = TRUE, eval = TRUE--------------------------------------------------------------------------
ensembl = useMart("ensembl_mart_46", dataset="hsapiens_gene_ensembl", archive = TRUE)

## ----archiveMarts, echo = TRUE, eval = TRUE---------------------------------------------------------------------------
listEnsemblArchives()

## ----archiveMarts3, echo = TRUE, eval = TRUE--------------------------------------------------------------------------
listMarts(host = 'may2009.archive.ensembl.org')
ensembl54 <- useMart(host='may2009.archive.ensembl.org', 
                     biomart='ENSEMBL_MART_ENSEMBL', 
                     dataset='hsapiens_gene_ensembl')

## ----wormbase, echo=TRUE, eval=TRUE-----------------------------------------------------------------------------------
listMarts(host = "parasite.wormbase.org")
wormbase = useMart(biomart = "parasite_mart", 
                   host = "https://parasite.wormbase.org", 
                   port = 443)
listDatasets(wormbase)
wormbase <- useDataset(mart = wormbase, dataset = "wbps_gene")
head(listFilters(wormbase))
head(listAttributes(wormbase))
getBM(attributes = c("external_gene_id", "wbps_transcript_id", "transcript_biotype"), 
      filters="gene_name", 
      values=c("unc-26","his-33"), 
      mart=wormbase)
     

## ----filterType-------------------------------------------------------------------------------------------------------
filterType("with_affy_hg_u133_plus_2",ensembl)

## ----filterOptions----------------------------------------------------------------------------------------------------
filterOptions("biotype",ensembl)

## ----attributePages---------------------------------------------------------------------------------------------------
pages = attributePages(ensembl)
pages

## ----listAttributes---------------------------------------------------------------------------------------------------
head(listAttributes(ensembl, page="feature_page"))

## ----localCopy, eval = FALSE------------------------------------------------------------------------------------------
#  listMarts(host="www.myLocalHost.org", path="/myPathToWebservice/martservice")
#  mart=useMart("nameOfMyMart",dataset="nameOfMyDataset",host="www.myLocalHost.org", path="/myPathToWebservice/martservice")

## ----columnsAndKeyTypes-----------------------------------------------------------------------------------------------
mart <- useMart(dataset="hsapiens_gene_ensembl",biomart='ensembl')
head(keytypes(mart), n=3)
head(columns(mart), n=3)

## ----keys1------------------------------------------------------------------------------------------------------------
k = keys(mart, keytype="chromosome_name")
head(k, n=3)

## ----keys2------------------------------------------------------------------------------------------------------------
k = keys(mart, keytype="chromosome_name", pattern="LRG")
head(k, n=3)

## ----select-----------------------------------------------------------------------------------------------------------
affy=c("202763_at","209310_s_at","207500_at")
select(mart, keys=affy, columns=c('affy_hg_u133_plus_2','entrezgene'),
  keytype='affy_hg_u133_plus_2')

## ----sessionInfo------------------------------------------------------------------------------------------------------
sessionInfo()
warnings()

