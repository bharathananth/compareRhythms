#' Run differential rhythmicity analysis for microarray using limma
#' @param eset A matrix with gene in the rows and samples in columns
#' @param exp_design A data.frame of the experimental design with at least
#'   columns sample name, time point and group
#' @param period The period of rhythm being tested (default = 24)
#' @param rhythm_fdr The false discovery cutoff for finding rhythmic time series
#'   (default = 0.05)
#' @param compare_fdr The false discovery cutoff for the comparison of rhythms
#'   (default = 0.05)
#' @param amp_cutoff The minimum peak-to-trough amp in log2 scale considered
#'   biologically relevant
#' @param just_classify Logical to select whether p-values, amplitudes and
#'   phases must be supressed in the results
#' @param rna_seq Indicates whether the data source is RNA-seq or microarray
#'    (default = "False")

compareRhythms_limma <- function(eset, exp_design, period, rhythm_fdr,
                                 compare_fdr, amp_cutoff, just_classify,
                                 rna_seq, robust) {

  group_id <- base::levels(exp_design$group)

  exp_design <- base::cbind(exp_design,
                            inphase = cos(2 * pi * exp_design$time / period),
                            outphase = sin(2 * pi * exp_design$time / period))

  if ("batch" %in% colnames(exp_design)) {

    design <- stats::model.matrix(~0 + group + group:inphase + group:outphase + batch,
                                  data = exp_design)
  } else {

    design <- stats::model.matrix(~0 + group + group:inphase + group:outphase,
                                  data = exp_design)
  }

  colnames(design) <- gsub("group", "", colnames(design))
  colnames(design) <- gsub(":", "_", colnames(design))

  fit <- limma::lmFit(eset, design)
  fit <- limma::eBayes(fit, robust = robust, trend = !rna_seq)

  rhythmic_in_either <- limma::topTable(fit,
                                        coef = grep("phase", colnames(design)),
                                        number = Inf,
                                        sort.by = "none")

  results <- compute_model_params(eset, group_id, design)

  results <- data.frame(results)

  results$symbol <- rownames(results)

  rownames(results) <- NULL

  results$max_amp <- pmax(results[, paste0(group_id[1], "_amp")],
                          results[, paste0(group_id[2], "_amp")])

  results$adj_p_val_A_or_B <- rhythmic_in_either$adj.P.Val

  results <- results[(results$adj_p_val_A_or_B < rhythm_fdr) &
                      (results$max_amp > amp_cutoff), ]

  assertthat::assert_that(assertthat::not_empty(results),
                          msg = "Sorry no rhythmic genes in either dataset for the thresholds provided.")

  results$max_amp <- NULL

  contrasts <- c(paste0(group_id, "_inphase", collapse = "-"),
                 paste0(group_id, "_outphase", collapse = "-"))

  diff_rhy_contrast <- limma::makeContrasts(contrasts = contrasts,
                                            levels = design)

  diff_rhy_fit <- limma::contrasts.fit(fit, diff_rhy_contrast)

  diff_rhy_fit <- limma::eBayes(diff_rhy_fit, robust = robust, trend = !rna_seq)

  diff_rhy_results <- limma::topTable(diff_rhy_fit, number = Inf,
                                      sort.by = "none")

  diff_rhy_results <- diff_rhy_results[rownames(diff_rhy_results) %in% results$symbol, ]

  results$adj_p_val_DR <- stats::p.adjust(diff_rhy_results$P.Value,
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
