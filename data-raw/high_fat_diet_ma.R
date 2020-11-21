suppressPackageStartupMessages(library(GEOquery))
library(magrittr)
library(oligo)
library(pd.mogene.1.0.st.v1)
library(mogene10sttranscriptcluster.db)
library(dplyr)

GEO <- "GSE52333"
chip <- "mogene10st"
org <- "mm"
annot <- "ensg"

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

platform_design <- paste0("pd.", chip, ".", org, ".", annot)

if (!require(platform_design, character.only = TRUE, quietly = TRUE,
             attach.required = FALSE, warn.conflicts = FALSE)) {
  link_to_cdf <- paste0("http://mbni.org/customcdf/",
                        version, "/", annot, ".download/",
                        platform_design, "_", version, ".tar.gz")

  utils::install.packages(link_to_cdf,
                          repos = NULL,
                          type = "source",
                          quiet = TRUE,
                          warn.conflicts = FALSE,
                          verbose = FALSE)
}


if (require(platform_design, character.only = TRUE, quietly = TRUE,
            attach.required = FALSE, warn.conflicts = FALSE)) {
  return(platform_design)
} else {
  stop("Installation failed.")
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


raw_data <- oligo::read.celfiles(cel_files_full,
                                 pkgname = platform_design)

eset <- oligo::rma(raw_data)

expr <- Biobase::exprs(eset)

exp_design <- Biobase::pData(gse) %>%
              group_by(`diet:ch1`, `harvest timepoint:ch1`) %>%
              mutate(rep=row_number())

colnames(expr) <- exp_design %$%
                      paste(`diet:ch1`, `harvest timepoint:ch1`, rep, sep="_") %>%
                      gsub("normal chow", "NC", .) %>% gsub("high fat diet", "HFD", .)

rownames(expr) <- gsub("_at", "", rownames(expr))

high_fat_diet_ma <- expr

usethis::use_data(high_fat_diet_ma, overwrite = TRUE)
