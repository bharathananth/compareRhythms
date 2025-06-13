#' Run differential rhythmicity analysis for RNA-seq data using DESeq2
#'
#' @inheritParams compareRhythms
#' @keywords internal

compareRhythms_deseq2 <- function(counts, exp_design, lengths, period,
                                  rhythm_fdr, compare_fdr, amp_cutoff,
                                  just_classify, just_rhythms) {

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
                                        colData = exp_design_aug,
                                        design = design)

  if (!is.null(lengths)) {
    SummarizedExperiment::assays(dds)[["avgTxLength"]] <- lengths
  }

  dds <- DESeq2::DESeq(dds, test = "LRT",
                       reduced = design[, !grepl("phase", colnames(design))],
                       quiet = TRUE)

  group_id <- base::levels(exp_design$group)

  rhythmic_in_either <- DESeq2::results(dds, independentFiltering = FALSE,
                                        cooksCutoff = FALSE)

  results <- compute_model_params(stats::coef(dds), group_id, type = "coef")

  results <- data.frame(results)

  results$id <- rownames(results)

  rownames(results) <- NULL

  results$max_amp <- pmax(results[, paste0(group_id[1], "_amp")],
                          results[, paste0(group_id[2], "_amp")])

  results$adj_p_val_A_or_B <- rhythmic_in_either$padj

  results <- results[(results$adj_p_val_A_or_B < rhythm_fdr) &
                       (results$max_amp >= amp_cutoff), ]

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


  diff_rhy_results <- DESeq2::results(diff_rhy_fit, independentFiltering = FALSE,
                                      cooksCutoff = FALSE)

  diff_rhy_results <- diff_rhy_results[results$id, ]

  results$adj_p_val_DR <- stats::p.adjust(diff_rhy_results$pvalue,
                                          method = "BH")
  results$diff_rhythmic <- results$adj_p_val_DR < compare_fdr

  results$rhythmic_in_A <- results[, paste0(group_id[1], "_amp")] > amp_cutoff

  results$rhythmic_in_B <- results[, paste0(group_id[2], "_amp")] > amp_cutoff

  assertthat::assert_that(
    assertthat::noNA(diff_rhy_results$pvalue),
    assertthat::noNA(diff_rhy_results$padj),
    assertthat::noNA(rhythmic_in_either$padj),
    msg = "DESeq2 is producing NA in p-values. Check if your count matrix has
    too many zero counts or too many outliers."
  )

  results$category <- base::mapply(categorize,
                                   results$rhythmic_in_A,
                                   results$rhythmic_in_B,
                                   results$diff_rhythmic)

  main_cols <- c("id", "category", "rhythmic_in_A", "rhythmic_in_B",
                 "diff_rhythmic")

  if (!just_rhythms) {
    group_id <- base::levels(exp_design$group)
    diff_exp_results <- DESeq2::results(dds,
                                        name = group_id[2],
                                        independentFiltering = FALSE,
                                        cooksCutoff = FALSE,
                                        alpha = rhythm_fdr)

    padj <- NULL
    diff_exp_results <- subset(diff_exp_results,
                               padj <= rhythm_fdr)

    diff_exp_results <- diff_exp_results[, c("log2FoldChange", "padj")]
    colnames(diff_exp_results) <- c("logFC_DE", "adj_p_val_DE")

    diff_exp_results$category_DE <- ifelse(diff_exp_results$logFC_DE>=0,
                                           "up-reg","down-reg")
    diff_exp_results$id <- rownames(diff_exp_results)

    rownames(diff_exp_results) <- NULL
    results <- merge(results, diff_exp_results, by="id", all=TRUE)

    main_cols <- c(main_cols, "category_DE")

    print(diff_exp_results)
  }

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
