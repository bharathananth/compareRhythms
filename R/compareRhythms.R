#' Run differential rhythmicity analysis
#' @param object A matrix of log2 expression values (if microarray), expression
#'   counts (RNA-seq) or a \linkS4class{SummarizedExperiment} object.
#' @param exp_design A data.frame of the experimental design with at least two
#'   columns: "time" and "group". Ignored if object is a
#'   \code{SummarizedExperiment} object.
#' @param lengths A data.frame of average transcript lengths. Only used with
#'   methods "deseq" and "edgeR".
#' @param method The method of analysis. It should be one of "mod_sel" for model
#'   selection, "dodr" for analysis using \code{\link[DODR]{dodr}}, "limma" for
#'   linear-modeling approach based on \pkg{limma}, "voom" for linear-modeling
#'   approah for RNA-Seq using \code{\link[limma]{voom}}, "deseq" for RNA-seq
#'   analysis using DESeq2, and "edgeR" for RNA-seq analysis using edgeR.
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
#'   different method (default = TRUE).
#' @param outliers Boolean specifying if weights must be computed for each
#'   sample to account for outliers. Only used by method = "voom".
#'
#' @export
compareRhythms <- function(object, exp_design=NULL, lengths=NULL,
                           method = "mod_sel", period=24, rhythm_fdr = 0.05,
                           compare_fdr = 0.05, amp_cutoff = 0.5,
                           criterion = "bic", schwarz_wt_cutoff = 0.6,
                           just_classify = TRUE, robust = TRUE, outliers = FALSE,
                           rna_seq = FALSE) {

  if (is.null(exp_design)) {
    assertthat::assert_that(
      attr(object, "class") %in% c("SummarizedExperiment", "RangedSummarizedExperiment"),
      msg = "Object does not belong to SummarizedExperiment.")
    assertthat::assert_that(
      assertthat::not_empty(object),
      assertthat::has_name(SummarizedExperiment::colData(object), c("time", "group")),
      assertthat::noNA(SummarizedExperiment::assay(object, 1)),
      is.factor(object$group),
      is.numeric(object$time),
      length(levels(object$group)) == 2
    )
    assertthat::assert_that(length(unique(SummarizedExperiment::colData(object)$group)) == 2,
                            msg = "Data does not have exactly two groups")
    expr <- SummarizedExperiment::assay(object, 1)
    exp_design <- SummarizedExperiment::colData(object)
    if ("lengths" %in% SummarizedExperiment::assayNames(object)) {
      lengths <- SummarizedExperiment::assay(object, "lengths")
    }
  } else {
    assertthat::assert_that(
      is.matrix(object),
      assertthat::not_empty(object),
      is.data.frame(exp_design),
      assertthat::not_empty(exp_design),
      assertthat::are_equal(ncol(object), nrow(exp_design)),
      assertthat::has_name(exp_design, c("time", "group")),
      assertthat::noNA(object),
      is.factor(exp_design$group),
      is.numeric(exp_design$time),
      length(levels(exp_design$group)) == 2
    )
    if (method %in% c("deseq", "edger") && !is.null(lengths)) {
      assertthat::assert_that(all(lengths>0), msg = "All transcript lengths are not positive")
    }
    expr <- object
  }

  assertthat::assert_that(
    assertthat::is.number(period),
    assertthat::is.number(rhythm_fdr),
    rhythm_fdr <= 1.0 & rhythm_fdr > 0,
    assertthat::is.number(compare_fdr),
    compare_fdr <= 1.0 & compare_fdr > 0,
    assertthat::is.number(amp_cutoff),
    amp_cutoff >= 0,
    assertthat::is.string(criterion),
    criterion %in% c("bic", "aic"),
    assertthat::is.number(schwarz_wt_cutoff),
    schwarz_wt_cutoff >=0 & schwarz_wt_cutoff <=1.0,
    assertthat::is.flag(just_classify),
    assertthat::is.flag(robust),
    assertthat::is.flag(outliers)
  )

  switch (method,
          mod_sel = compareRhythms_model_select(expr = expr,
                                                exp_design = exp_design,
                                                period = period,
                                                amp_cutoff = amp_cutoff,
                                                criterion = criterion,
                                                schwarz_wt_cutoff = schwarz_wt_cutoff,
                                                just_classify = just_classify),
          dodr = compareRhythms_dodr(expr = expr,
                                     exp_design = exp_design,
                                     period = period,
                                     rhythm_fdr = rhythm_fdr,
                                     compare_fdr = compare_fdr,
                                     amp_cutoff = amp_cutoff,
                                     just_classify = just_classify),
          limma = compareRhythms_limma(eset = expr,
                                       exp_design = exp_design,
                                       period = period,
                                       rhythm_fdr = rhythm_fdr,
                                       amp_cutoff = amp_cutoff,
                                       compare_fdr = compare_fdr,
                                       just_classify = just_classify,
                                       rna_seq = rna_seq,
                                       robust = robust),
          voom = compareRhythms_voom(counts = expr,
                                     exp_design = exp_design,
                                     period = period,
                                     rhythm_fdr = rhythm_fdr,
                                     amp_cutoff = amp_cutoff,
                                     compare_fdr = compare_fdr,
                                     just_classify = just_classify,
                                     robust = robust,
                                     outliers = outliers),
          deseq = compareRhythms_deseq2(counts = expr, exp_design = exp_design,
                                        lengths = lengths, period = period, rhythm_fdr = rhythm_fdr,
                                        amp_cutoff = amp_cutoff, compare_fdr = compare_fdr,
                                        just_classify = just_classify),
          edger = compareRhythms_edger(counts = expr, exp_design = exp_design,
                                        lengths = lengths, period = period, rhythm_fdr = rhythm_fdr,
                                        amp_cutoff = amp_cutoff, compare_fdr = compare_fdr,
                                        just_classify = just_classify)
  )
}
