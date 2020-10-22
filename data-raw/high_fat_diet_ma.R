## code to prepare `high_fat_diet_ma` dataset goes here
suppressPackageStartupMessages(library(GEOquery))
library(magrittr)
library(oligo)
library(pd.mogene.1.0.st.v1)
library(mogene10sttranscriptcluster.db)

GEO <- "GSE52333"
gse <- getGEO(GEO)[[1]]

base_dir <- paste(base::getwd(), "data-raw", sep="/")

if (!base::dir.exists(base_dir)) base::dir.create(base_dir)

files <- GEOquery::getGEOSuppFiles(GEO,
                                   baseDir = base_dir,
                                   filter_regex = "RAW",
                                   fetch_files = FALSE)
raw_data_file <- paste(base_dir, GEO, files$fname, sep = "/")

if (!all(file.exists(raw_data_file))) {
  GEOquery::getGEOSuppFiles(GEO,
                            baseDir = base_dir,
                            filter_regex = "RAW",
                            fetch_files = TRUE)
}

cel_files <- as.character(gse[["supplementary_file"]])

err_file_name <- cel_files[grepl("GSM1263252", cel_files)]
cel_files[grepl("GSM1263252", cel_files)] <- sub("pli[\\w\\-\\.]+", "CEL.gz", sub("bA", "bB", err_file_name), perl=TRUE)

cel_files <- sapply(cel_files, basename)

cel_files_full <- paste(dirname(raw_data_file), cel_files, sep="/")

if (!all(file.exists(cel_files_full))) {

  utils::untar(raw_data_file,
               files = cel_files,
               exdir = dirname(raw_data_file))
}


raw_data <- oligo::read.celfiles(cel_files_full)

eset <- oligo::rma(raw_data)

expr <- Biobase::exprs(eset)

exp_design <- Biobase::pData(gse)

colnames(expr) <- exp_design %$%
                      paste(`diet:ch1`, `harvest timepoint:ch1`, sep="_") %>%
                      gsub("normal chow", "NC", .) %>% gsub("high fat diet", "HFD", .)

collapser <- function(x){
  x %>% unique %>% sort %>% paste(collapse = "|")
}

annots <- AnnotationDbi::select(
  x       = mogene10sttranscriptcluster.db,
  keys    = rownames(expr),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
) %>%
  group_by(PROBEID) %>%
  summarize(across(.fns = collapser)) %>%
  ungroup



usethis::use_data(high_fat_diet_ma, overwrite = TRUE)
