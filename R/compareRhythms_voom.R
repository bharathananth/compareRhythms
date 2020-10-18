#' Run differential rhythmicity analysis for RNA-seq data using limma-voom
#'
#' @param counts A matrix of  expression counts
#' @inheritParams compareRhythms
#' @keywords internal

compareRhythms_voom <- function(counts, exp_design, period, rhythm_fdr,
                                compare_fdr, amp_cutoff, just_classify,
                                robust, outliers) {

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

  y <- edgeR::DGEList(counts)

  y <- edgeR::calcNormFactors(y)

  if (outliers) {
    v <- limma::voomWithQualityWeights(y, design)
  } else {
    v <- limma::voom(y, design)
  }

  results <- compareRhythms_limma(v, exp_design, period = period, rhythm_fdr = rhythm_fdr,
                       compare_fdr = compare_fdr, amp_cutoff = amp_cutoff,
                       just_classify = just_classify, robust = robust, rna_seq = TRUE)

  return(results)
}
