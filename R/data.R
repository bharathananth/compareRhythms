#' Effect of high fat diet on the liver transcriptome - Microarray
#'
#' Time-series data of the mouse liver transcriptome measured under normal chow
#' read using custom CDF file from Brainarray version 24.0.0 annotated by
#' Ensembl Gene ID followed by RMA normalization. The expression data are in the
#' log2 scale. The sample specification are present in the column names in a
#' `<group>`_ZT`<time>` format.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52333}
"high_fat_diet_ma"


#' Effect of high fat diet on the liver transcriptome - RNA-sequencing
#'
#' Time-series data of the mouse liver transcriptome measured under normal chow
#' (NC) and high fat diet (HFD). The RNA-seq data aligned using STAR to mm9
#' genome and quantified using featureCounts (see GSE108688 for details) and
#' annotated by Ensembl Gene ID. The sample specification are present in the
#' column names in a `<group>`_ZT`<time>`_`<replicate>` format.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108688}
"high_fat_diet_rnaseq"
