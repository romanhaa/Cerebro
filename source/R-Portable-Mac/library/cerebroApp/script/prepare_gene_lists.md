# How to generate lists of ribosomal and mitochondrial genes

What is needed: Gene annotation as `GTF` file for human and mouse.
We use the one from GENCODE: <https://www.gencodegenes.org/>

```r
library('tidyverse')
```

## Human (hg38)

```r
gtf <- rtracklayer::import('gencode.v27.annotation.gtf')
```

### Mitochondrial genes

```r
t <- gtf %>%
  as_tibble() %>%
  filter(
    type == 'gene',
    grepl(gene_name, pattern = '^MT-')
  )

t %>%
select(gene_name) %>%
distinct() %>%
write_tsv('genes_mt_hg_name.tsv.gz', col_names = FALSE)

t %>%
select(gene_id) %>%
distinct() %>%
write_tsv('genes_mt_hg_gencode_v27.tsv.gz', col_names = FALSE)

t %>%
select(gene_id) %>%
mutate(gene_id = gsub(gene_id, pattern = '\\.[0-9]{1,3}', replacement = '')) %>%
distinct() %>%
write_tsv('genes_mt_hg_ensembl.tsv.gz', col_names = FALSE)
```

### Ribosomal genes

```r
t <- gtf %>%
  as_tibble() %>%
  filter(
    type == 'gene',
    grepl(gene_name, pattern = '^RP[SL]')
  )

t %>%
select(gene_name) %>%
distinct() %>%
write_tsv('genes_ribo_hg_name.tsv.gz', col_names = FALSE)

t %>%
select(gene_id) %>%
distinct() %>%
write_tsv('genes_ribo_hg_gencode_v27.tsv.gz', col_names = FALSE)

t %>%
select(gene_id) %>%
mutate(gene_id = gsub(gene_id, pattern = '\\.[0-9]{1,3}', replacement = '')) %>%
distinct() %>%
write_tsv('genes_ribo_hg_ensembl.tsv.gz', col_names = FALSE)
```

## Mouse (mm10)

```r
gtf <- rtracklayer::import('gencode.vM16.comprehensive.CHR.gtf')
```

### Mitochondrial genes

```r
t <- gtf %>%
  as_tibble() %>%
  filter(
    type == 'gene',
    grepl(gene_name, pattern = '^mt-')
  )

t %>%
select(gene_name) %>%
distinct() %>%
write_tsv('genes_mt_mm_name.tsv.gz', col_names = FALSE)

t %>%
select(gene_id) %>%
distinct() %>%
write_tsv('genes_mt_mm_gencode_vM16.tsv.gz', col_names = FALSE)

t %>%
select(gene_id) %>%
mutate(gene_id = gsub(gene_id, pattern = '\\.[0-9]{1,3}', replacement = '')) %>%
distinct() %>%
write_tsv('genes_mt_mm_ensembl.tsv.gz', col_names = FALSE)
```

### Ribosomal genes

```r
t <- gtf %>%
  as_tibble() %>%
  filter(
    type == 'gene',
    grepl(gene_name, pattern = '^Rp[sl]')
  )

t %>%
select(gene_name) %>%
distinct() %>%
write_tsv('genes_ribo_mm_name.tsv.gz', col_names = FALSE)

t %>%
select(gene_id) %>%
distinct() %>%
write_tsv('genes_ribo_mm_gencode_vM16.tsv.gz', col_names = FALSE)

t %>%
select(gene_id) %>%
mutate(gene_id = gsub(gene_id, pattern = '\\.[0-9]{1,3}', replacement = '')) %>%
distinct() %>%
write_tsv('genes_ribo_mm_ensembl.tsv.gz', col_names = FALSE)
```
