#' Run differential rhythmicity analysis for RNA-seq data using edgeR
#'
#' @inheritParams compareRhythms
#' @keywords internal

compareRhythms_edgeR <- function(counts, exp_design, lengths, period,
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

  if (is.null(lengths)) {
    lengths <- matrix(1, NROW(counts), NCOL(counts))
  }

  normMat <- lengths
  normMat <- normMat/exp(base::rowMeans(log(normMat)))
  normCts <- counts/normMat

  eff.lib <- edgeR::calcNormFactors(normCts) * base::colSums(normCts)

  normMat <- sweep(normMat, 2, eff.lib, "*")
  normMat <- log(normMat)

  y <- edgeR::DGEList(counts)
  y <- edgeR::scaleOffset(y, normMat)

  y <- edgeR::estimateDisp(y, design, robust = TRUE)

  fit <- edgeR::glmQLFit(y, design, robust = TRUE)

  group_id <- base::levels(exp_design$group)

  qlf <- edgeR::glmQLFTest(fit, coef = grep("phase", colnames(design)))

  rhythmic_in_either <- edgeR::topTags(qlf, n = Inf, sort.by = "none")

  rhythmic_in_either <- data.frame(rhythmic_in_either)

  colnames(rhythmic_in_either) <- base::gsub("logFC.", "", colnames(rhythmic_in_either))

  results <- compute_model_params(rhythmic_in_either, group_id, type = "coef")

  results <- data.frame(results)

  results$id <- rownames(results)

  rownames(results) <- NULL

  results$max_amp <- pmax(results[, paste0(group_id[1], "_amp")],
                          results[, paste0(group_id[2], "_amp")])

  results$adj_p_val_A_or_B <- rhythmic_in_either$FDR

  results <- results[(results$adj_p_val_A_or_B < rhythm_fdr) &
                       (results$max_amp > amp_cutoff), ]

  assertthat::assert_that(assertthat::not_empty(results),
                          msg = "Sorry no rhythmic genes in either dataset for the thresholds provided.")

  results$max_amp <- NULL

  contrasts <- c(paste0(group_id, "_inphase", collapse = "-"),
                 paste0(group_id, "_outphase", collapse = "-"))

  diff_rhy_contrast <- limma::makeContrasts(contrasts = contrasts,
                                            levels = design)

  diff_rhy_fit <- edgeR::glmQLFTest(fit, contrast = diff_rhy_contrast)

  diff_rhy_results <- edgeR::topTags(diff_rhy_fit, n = Inf, sort.by = "none")

  diff_rhy_results <- data.frame(diff_rhy_results)

  diff_rhy_results <- diff_rhy_results[results$id, ]

  results$adj_p_val_DR <- stats::p.adjust(diff_rhy_results$PValue,
                                          method = "BH")
  results$diff_rhythmic <- results$adj_p_val_DR < compare_fdr

  results$rhythmic_in_A <- results[, paste0(group_id[1], "_amp")] > amp_cutoff

  results$rhythmic_in_B <- results[, paste0(group_id[2], "_amp")] > amp_cutoff

  results$category <- base::mapply(categorize,
                                   results$rhythmic_in_A,
                                   results$rhythmic_in_B,
                                   results$diff_rhythmic)

  main_cols <- c("id", "category", "rhythmic_in_A", "rhythmic_in_B",
                 "diff_rhythmic")

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
