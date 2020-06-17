# How to generate gene ID/name conversion tables

What is needed: Gene annotation as `GTF` file for human and mouse.
We use the one from GENCODE: <https://www.gencodegenes.org/>

## Human (hg38)

```r
library('tidyverse')

gtf <- rtracklayer::import('gencode.v27.annotation.gtf')

t <- gtf %>%
  as_tibble() %>%
  select(gene_id, havana_gene, gene_name, gene_type) %>%
  rename(
    gencode = gene_id,
    havana = havana_gene,
    symbol = gene_name,
    type = gene_type
  ) %>%
  mutate(ensembl = gsub(gencode, pattern = '\\.[0-9]{1,3}', replacement = '')) %>%
  distinct() %>%
  select(gencode, ensembl, havana, symbol, type)

write_tsv(t, 'hg38_gene_ID_name.tsv.gz')
```

## Mouse (mm10)

```r
library('tidyverse')

gtf <- rtracklayer::import('gencode.vM16.comprehensive.CHR.gtf')

t <- gtf %>%
  as_tibble() %>%
  select(gene_id, havana_gene, gene_name, gene_type) %>%
  rename(
    gencode = gene_id,
    havana = havana_gene,
    symbol = gene_name,
    type = gene_type
  ) %>%
  mutate(ensembl = gsub(gencode, pattern = '\\.[0-9]{1,3}', replacement = '')) %>%
  distinct() %>%
  select(gencode, ensembl, havana, symbol, type)

write_tsv(t, 'mm10_gene_ID_name.tsv.gz')
```
