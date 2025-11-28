# Effect of high fat diet on the liver transcriptome - RNA-sequencing

Time-series data of the mouse liver transcriptome measured under normal
chow (NC) and high fat diet (HFD). The RNA-seq data aligned using STAR
to mm9 genome and quantified using featureCounts (see GSE108688 for
details) and annotated by Ensembl Gene ID. The sample specification are
present in the column names in a `<group>`*ZT`<time>`*`<replicate>`
format.

## Usage

``` r
high_fat_diet_rnaseq
```

## Format

An object of class `matrix` (inherits from `array`) with 37310 rows and
36 columns.

## Source

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108688>
