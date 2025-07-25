context("compareRhythms_deseq2")

load("test_data_rnaseq.rda")

exp_design_batch <- cbind(exp_design, batch, stringsAsFactors=TRUE)

test_that("DESeq2 analysis works for default params", {
  results <- compareRhythms(countsFromAbundance, exp_design, method = "deseq2")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic"))
  results <- compareRhythms(counts, exp_design, lengths = lengths, method = "deseq2")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic"))

  results <- compareRhythms(counts, exp_design_batch, lengths = lengths, method = "deseq2")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic"))

  results <- compareRhythms(countsFromAbundance, exp_design_batch, method = "deseq2")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic"))
})

test_that("DESeq2 analysis works for different input params", {
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design, period = 12, method = "deseq2"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design, outliers = TRUE, method = "deseq2"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design, rhythm_fdr = 0.01, method = "deseq2"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design, compare_fdr = 0.01, method = "deseq2"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design, rhythm_fdr = 0.1, compare_fdr = 0.01, method = "deseq2"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design, amp_cutoff = 0, method = "deseq2"), "data.frame")
  expect_error(compareRhythms(countsFromAbundance, exp_design, rhythm_fdr = 0, method = "deseq2"))
  results <- compareRhythms(countsFromAbundance, exp_design, just_classify = FALSE, method = "deseq2")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic", "P66KO_amp",
                 "P66KO_phase", "WT_amp", "WT_phase", "adj_p_val_P66KO_or_WT",
                 "adj_p_val_DR"))
  results <- compareRhythms(countsFromAbundance, exp_design_batch, method = "deseq2", just_rhythms=FALSE)
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic", "category_DE"))
  results <- compareRhythms(countsFromAbundance, exp_design_batch, method = "deseq2", just_classify=FALSE, just_rhythms=FALSE)
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic", "category_DE", "P66KO_amp",
                 "P66KO_phase", "WT_amp", "WT_phase", "adj_p_val_P66KO_or_WT",
                 "adj_p_val_DR", "logFC_DE", "adj_p_val_DE"))
})
