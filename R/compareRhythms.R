#' Run differential rhythmicity analysis
#'
#' The differential rhythmicity analysis is run with a call to this function. To
#' execute this function, the three necessary ingredients are the timeseries
#' data, the experimental design and parameters to choose and tune the method.
#'
#' @param data A matrix of log2 expression values (if microarray), expression
#'   counts (RNA-seq) or normalized data (see Details).
#' @param exp_design A data.frame of the experimental design with at least two
#'   columns: "time" and "group".
#' @param lengths A data.frame of average transcript lengths. Only used with
#'   methods "deseq" and "edgeR".
#' @param method The method of analysis. It should be one of "mod_sel" for model
#'   selection, "dodr" for analysis using [DODR::dodr], "limma" for
#'   linear-modeling approach based on \pkg{limma}, "voom" for linear-modeling
#'   approach for RNA-Seq using [limma::voom], "deseq2" for RNA-seq analysis
#'   using \pkg{DESeq2}, "edger" for RNA-seq analysis using \pkg{edgeR}, and
#'   "cosinor" for simple cosinor-based analysis both for independent samples
#'   and repeated samples.
#' @param period The period of rhythm being tested (default = 24)
#' @param rhythm_fdr The false discovery cutoff for finding rhythmic time series
#'   (default = 0.05)
#' @param compare_fdr The false discovery cutoff for the comparison of rhythms
#'   (default = 0.05)
#' @param amp_cutoff The minimum peak-to-trough amplitude in log2 scale
#'   considered biologically relevant (default = 0.5)
#' @param criterion The criterion used for model selection. These can be "aic"
#'   or "bic" (default = "bic"). Only used for method = "mod_sel".
#' @param schwarz_wt_cutoff The conditional probability that the best model is
#'   the true model. Genes with a conditional probability smaller than this
#'   cutoff are deemed unclassifiable. This is only used for method = "mod_sel".
#'   (default = 0.4)
#' @param just_classify Boolean specifying whether genes must only be classified
#'   (TRUE) or if the amplitude and phases of fits should also be returned
#'   (FALSE)
#' @param robust Boolean to turn on robust computation of statistics in
#'   different methods (default = TRUE).
#' @param outliers Boolean specifying if weights must be computed for each
#'   sample to account for outliers. Only used by method = "voom".
#' @param longitudinal Boolean specifying if repeated samples from one
#'   experimental unit. Only used by method = "cosinor".
#' @param just_rhythms Boolean specifying whether only rhythm analysis (True,
#'   default) or differential expression/magnitude analysis should also be
#'   performed.
#'
#' @return A *data.frame* with the names of the differentially rhythmic
#'   features, the category it is classified under and optionally the rhythm
#'   parameters of the features in each group. The differential rhythmicity
#'   categories are **gain** of, **loss** of, **change** of, or **same** rhythms
#'   (with respect to the reference/control group).
#'
#' @export
compareRhythms <- function(data, exp_design, lengths=NULL,
                           method = "mod_sel", period=24, rhythm_fdr = 0.05,
                           compare_fdr = 0.05, amp_cutoff = 0.5,
                           criterion = "bic", schwarz_wt_cutoff = 0.6,
                           just_classify = TRUE, robust = TRUE, outliers = FALSE,
                           longitudinal = FALSE, just_rhythms = TRUE
                           ) {

  assertthat::assert_that(
    is.matrix(data),
    assertthat::not_empty(data),
    !is.null(rownames(data)),
    is.data.frame(exp_design),
    assertthat::not_empty(exp_design),
    assertthat::are_equal(ncol(data), nrow(exp_design)),
    assertthat::has_name(exp_design, c("time", "group")),
    method %in% c("mod_sel", "limma", "dodr", "voom", "deseq2", "edger","cosinor"),
    is.factor(exp_design$group),
    is.numeric(exp_design$time),
    length(levels(exp_design$group)) == 2,
    assertthat::is.number(period),
    period > 0,
    assertthat::is.number(rhythm_fdr),
    rhythm_fdr <= 1.0 & rhythm_fdr >= 0,
    assertthat::is.number(compare_fdr),
    compare_fdr <= 1.0 & compare_fdr >= 0,
    assertthat::is.number(amp_cutoff),
    amp_cutoff >= 0,
    assertthat::is.string(criterion),
    criterion %in% c("bic", "aic"),
    assertthat::is.number(schwarz_wt_cutoff),
    schwarz_wt_cutoff >=0 & schwarz_wt_cutoff <=1.0,
    assertthat::is.flag(just_classify),
    assertthat::is.flag(robust),
    assertthat::is.flag(outliers),
    assertthat::is.flag(longitudinal),
    assertthat::is.flag(just_rhythms)
  )
  if (method != "cosinor"){
    assertthat::noNA(data)
  }
  if (method %in% c("deseq", "edger") && !is.null(lengths)) {
    assertthat::assert_that(all(lengths>0), msg = "All transcript lengths are not positive")
  }

  if (method == "cosinor" && longitudinal) {
    assertthat::assert_that(assertthat::has_name(exp_design, "ID"))
    if (assertthat::has_name(exp_design, "ID")) {
      assertthat::assert_that(is.factor(exp_design$ID))
    }
  }

  if (assertthat::has_name(exp_design, "batch")) {
    assertthat::assert_that(is.factor(exp_design$batch))
  }


  switch (method,
          mod_sel = compareRhythms_model_select(data = data,
                                                exp_design = exp_design,
                                                period = period,
                                                amp_cutoff = amp_cutoff,
                                                criterion = criterion,
                                                schwarz_wt_cutoff = schwarz_wt_cutoff,
                                                just_classify = just_classify),
          dodr = compareRhythms_dodr(expr = data,
                                     exp_design = exp_design,
                                     period = period,
                                     rhythm_fdr = rhythm_fdr,
                                     compare_fdr = compare_fdr,
                                     amp_cutoff = amp_cutoff,
                                     just_classify = just_classify),
          limma = compareRhythms_limma(eset = data,
                                       exp_design = exp_design,
                                       period = period,
                                       rhythm_fdr = rhythm_fdr,
                                       amp_cutoff = amp_cutoff,
                                       compare_fdr = compare_fdr,
                                       just_classify = just_classify,
                                       rna_seq = FALSE,
                                       robust = robust,
                                       just_rhythms=just_rhythms),
          voom = compareRhythms_voom(counts = data,
                                     exp_design = exp_design,
                                     period = period,
                                     rhythm_fdr = rhythm_fdr,
                                     amp_cutoff = amp_cutoff,
                                     compare_fdr = compare_fdr,
                                     just_classify = just_classify,
                                     robust = robust,
                                     outliers = outliers,
                                     just_rhythms = just_rhythms),
          deseq2 = compareRhythms_deseq2(counts = data,
                                         exp_design = exp_design,
                                         lengths = lengths,
                                         period = period,
                                         rhythm_fdr = rhythm_fdr,
                                         amp_cutoff = amp_cutoff,
                                         compare_fdr = compare_fdr,
                                         just_classify = just_classify,
                                         just_rhythms = just_rhythms),
          edger = compareRhythms_edgeR(counts = data,
                                       exp_design = exp_design,
                                       lengths = lengths,
                                       period = period,
                                       rhythm_fdr = rhythm_fdr,
                                       amp_cutoff = amp_cutoff,
                                       compare_fdr = compare_fdr,
                                       just_classify = just_classify,
                                       just_rhythms = just_rhythms),
          cosinor = compareRhythms_cosinor(data = data,
                                           exp_design = exp_design,
                                           period = period,
                                           rhythm_fdr = rhythm_fdr,
                                           compare_fdr = compare_fdr,
                                           amp_cutoff = amp_cutoff,
                                           just_classify = just_classify,
                                           longitudinal = longitudinal)
  )
}
