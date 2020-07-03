#' Run differential rhythmicity analysis on microarray or RNA-seq data
#' @param object A matrix of log2 expression values (if microarray), expression
#'   counts (RNA-seq) or a \linkS4class{SummarizedExperiment} object.
#' @param exp_design A data.frame of the experimental design with at least two
#'   columns: "time" and "group". Ignored if object is a
#'   \code{SummarizedExperiment} object.
#' @param method The method of analysis. It should be one of "mod_sel": model
#'   selection, "dodr": based on , limma", "voom".
#' @param period The period of rhythm being tested (default = 24)
#' @param rhythm_fdr The false discovery cutoff for finding rhythmic time series
#'   (default = 0.05)
#' @param compare_fdr The false discovery cutoff for the comparison of rhythms
#'   (default = 0.05)
#' @param amp_cutoff The minimum peak-to-trough amp in log2 scale considered
#'   biologically relevant (default = 0.5)
#' @param criterion The criterion used for model selection. These can be "aic"
#'   or "bic" (default = "bic"). This is only used for method = "mod_sel".
#' @param schwartz_wt_cutoff The conditional probability that the best model is
#'   the true model. Genes with a smaller conditional probability smaller than
#'   this cutoff are deemed unclassifiable. This is only used for method =
#'   "mod_sel". (default = 0.4)
#' @param just_classify Boolean specifying whether genes must only be classified
#'   (TRUE) or if the amplitude and phases of fits should also be returned
#'   (FALSE)
#' @param outliers Boolean specifying if weights must be computed for each
#'   sample to account for outliers.
#'
#' @export
compareRhythms <- function(object, exp_design=NULL, method = "mod_sel",
                           period=24, rhythm_fdr = 0.05,
                           compare_fdr = 0.05, amp_cutoff = 0.5,
                           criterion = "bic", schwartz_wt_cutoff = 0.6,
                           just_classify = TRUE, outliers = FALSE) {

  if (is.null(exp_design)) {
    assertthat::assert_that(
      attr(object, "class") %in% c("SummarizedExperiment", "RangedSummarizedExperiment"),
      msg = "Object does not belong to SummarizedExperiment.")
    assertthat::assert_that(
      assertthat::not_empty(object),
      assertthat::has_name(SummarizedExperiment::colData(object), c("time", "group")),
      assertthat::noNA(SummarizedExperiment::assay(object, 1))
    )
    assertthat::assert_that(length(unique(SummarizedExperiment::colData(object)$group)) == 2,
                       msg = "Data does not have exactly two groups")
    expr <- SummarizedExperiment::assay(object, 1)
    exp_design <- SummarizedExperiment::colData(object)
  } else {
    assertthat::assert_that(
      is.matrix(object),
      assertthat::not_empty(object),
      is.data.frame(exp_design),
      assertthat::not_empty(exp_design),
      assertthat::are_equal(ncol(object), nrow(exp_design)),
      assertthat::has_name(exp_design, c("time", "group")),
      assertthat::noNA(object),
      length(unique(exp_design$group)) == 2
    )
    expr <- object
  }

  switch (method,
    mod_sel = compareRhythms_model_select(expr = expr,
                                            exp_design = exp_design,
                                            period = period,
                                            amp_cutoff = amp_cutoff,
                                            criterion = criterion,
                                            schwartz_wt_cutoff = schwartz_wt_cutoff,
                                            just_classify = just_classify),
    dodr = compareRhythms_rain(expr = expr,
                               exp_design = exp_design, period = period,
                               rhythm_fdr = rhythm_fdr,
                               compare_fdr = compare_fdr,
                               amp_cutoff = amp_cutoff,
                               just_classify = just_classify),
    limma = compareRhythms_limma(eset = expr, exp_design = exp_design,
                                 period = period, rhythm_fdr = rhythm_fdr,
                                 amp_cutoff = amp_cutoff,
                                 compare_fdr = compare_fdr,
                                 just_classify = just_classify),
    voom = compareRhythms_voom(counts = expr, exp_design = exp_design,
                                period = period, rhythm_fdr = rhythm_fdr,
                                amp_cutoff = amp_cutoff,
                                just_classify = just_classify, outliers = outliers)
  )
}
