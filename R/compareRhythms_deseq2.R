#' Run differential rhythmicity analysis for RNA-seq data using DESeq2
#' @param counts A count matrix with gene in the rows and samples in columns
#' @param exp_design A data.frame of the experimental design with at least
#'   columns sample name, time point and group
#' @param lengths Average transcript length for the gene in each sample.
#'   Obtained from tximport (default = NULL)
#' @param period The period of rhythm being tested (default = 24)
#' @param rhythm_fdr The false discovery cutoff for finding rhythmic time series
#'   (default = 0.05)
#' @param compare_fdr The false discovery cutoff for the comparison of rhythms
#'   (default = 0.05)
#' @param amp_cutoff The minimum peak-to-trough amp in log2 scale considered
#'   biologically relevant
#' @param just_classify Logical to select whether p-values, amplitudes and
#'   phases must be supressed in the results

compareRhythms_deseq2 <- function(counts, exp_design, lengths, period,
                                  rhythm_fdr, compare_fdr, amp_cutoff,
                                  just_classify) {

  exp_design_aug <- base::cbind(exp_design,
                                inphase = cos(2 * pi * exp_design$time / period),
                                outphase = sin(2 * pi * exp_design$time / period))

  if ("batch" %in% colnames(exp_design)) {

    design <- stats::model.matrix(~group + group:inphase + group:outphase + batch,
                                  data = exp_design_aug)

  } else {

    design <- stats::model.matrix(~group + group:inphase + group:outphase,
                                  data = exp_design_aug)
  }

  colnames(design) <- gsub("group", "", colnames(design))
  colnames(design) <- gsub(":", "_", colnames(design))

  mode(counts) <- "integer"

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                        colData = exp_design_aug, design = design)

  if (!is.null(lengths)) {
    SummarizedExperiment::assays(dds)[["avgTxLength"]] <- lengths
  }

  dds <- DESeq2::DESeq(dds, test = "LRT",
                       reduced = design[, !grepl("phase", colnames(design))],
                       quiet = TRUE)

  group_id <- base::levels(exp_design$group)

  rhythmic_in_either <- DESeq2::results(dds, independentFiltering = FALSE)

  results <- compute_model_params(coef(dds), group_id, type = "coef")

  results <- data.frame(results)

  results$symbol <- rownames(results)

  rownames(results) <- NULL

  results$max_amp <- pmax(results[, paste0(group_id[1], "_amp")],
                          results[, paste0(group_id[2], "_amp")])

  results$adj_p_val_A_or_B <- rhythmic_in_either$padj

  results <- results[(results$adj_p_val_A_or_B < rhythm_fdr) &
                       (results$max_amp > amp_cutoff), ]

  assertthat::assert_that(assertthat::not_empty(results),
                          msg = "Sorry no rhythmic genes in either dataset for the thresholds provided.")

  results$max_amp <- NULL

  if ("batch" %in% colnames(exp_design)) {

    design <- stats::model.matrix(~group + group*(inphase + outphase) + batch,
                                  data = exp_design_aug, contrasts=list(group="contr.treatment"))

  } else {

    design <- stats::model.matrix(~group + group*(inphase + outphase),
                                  data = exp_design_aug, contrasts=list(group="contr.treatment"))
  }


  diff_rhy_fit <- DESeq2::DESeq(dds, test = "LRT",
                                reduced = design[, !grepl(":\\w+phase", colnames(design))],
                                quiet = TRUE)


  diff_rhy_results <- DESeq2::results(diff_rhy_fit, independentFiltering = FALSE)

  diff_rhy_results <- diff_rhy_results[rownames(diff_rhy_results) %in% results$symbol, ]

  results$adj_p_val_DR <- stats::p.adjust(diff_rhy_results$pvalue,
                                          method = "BH")
  results$diff_rhythmic <- results$adj_p_val_DR < compare_fdr

  results$rhythmic_in_A <- results[, paste0(group_id[1], "_amp")] > amp_cutoff

  results$rhythmic_in_B <- results[, paste0(group_id[2], "_amp")] > amp_cutoff

  results$category <- base::mapply(categorize,
                                   results$rhythmic_in_A,
                                   results$rhythmic_in_B,
                                   results$diff_rhythmic)

  main_cols <- c("symbol", "rhythmic_in_A", "rhythmic_in_B",
                 "diff_rhythmic", "category")

  results <- results[, c(main_cols,
                         base::setdiff(colnames(results), main_cols))]

  if (just_classify) {
    results <- results[, main_cols]
  }

  rownames(results) <- NULL
  colnames(results) <- gsub("A", group_id[1], colnames(results))
  colnames(results) <- gsub("B", group_id[2], colnames(results))

  return(results)
}
