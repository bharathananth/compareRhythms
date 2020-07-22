#' Run differential rhythmicity analysis for RNA-seq data using limma-voom
#' @param counts A count matrix with gene in the rows and samples in columns
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
#' @param outliers Boolean specifying if weights must be computed for each
#'   sample to account for outliers.

compareRhythms_voom <- function(counts, exp_design, period=24, rhythm_fdr = 0.05,
                                 compare_fdr = 0.05, amp_cutoff = 0.5,
                                 just_classify = TRUE, outliers = FALSE) {

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
                       just_classify = just_classify, rna_seq = TRUE)

  return(results)
}
