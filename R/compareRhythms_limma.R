
#' @export
compareRhythms_limma <- function(y, exp_design, period=24, rhythm_fdr = 0.05,
                                 compare_fdr = 0.05, amp_cutoff = 0.5,
                                 include_pvals = FALSE) {

  group_id <- base::unique(exp_design$group)

  exp_design <- base::cbind(exp_design,
                            inphase = cos(2 * pi * exp_design$time / period),
                            outphase = sin(2 * pi * exp_design$time / period))

  design <- stats::model.matrix(~0 + group + group:inphase + group:outphase,
                                data = exp_design)
  colnames(design) <- gsub("group", "", colnames(design))
  colnames(design) <- gsub(":", "_", colnames(design))

  fit <- limma::lmFit(y, design)
  fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)

  rhythms_A_or_B <- limma::topTable(fit,
                                    coef = grep("phase", colnames(design)),
                                    number = Inf,
                                    sort.by = "none")

  rhythms_A_or_B <- rhythms_A_or_B[rhythms_A_or_B$adj.P.Val < rhythm_fdr, ]

  contrasts <- c(paste0(group_id, "_inphase", collapse = "-"),
                 paste0(group_id, "_outphase", collapse = "-"))

  diff_rhy_contrast <- limma::makeContrasts(contrasts = contrasts,
                                            levels = design)

  diff_rhy_fit <- limma::contrasts.fit(fit, diff_rhy_contrast)

  diff_rhy_fit <- limma::eBayes(diff_rhy_fit, robust = TRUE, trend = TRUE)

  diff_rhy_results <- limma::decideTests(diff_rhy_fit[rownames(rhythms_A_or_B), ],
                                         method = "global",
                                         p.value = compare_fdr,
                                         adjust.method = "BH")

  results <- rhythms_A_or_B
  results$diff_rhythmic <- base::rowSums(abs(diff_rhy_results)) > 0
  return(results)
}
