context("compareRhythms_edgeR")

load("test_data_rnaseq.rda")

exp_design_batch <- cbind(exp_design, batch, stringsAsFactors=TRUE)

test_that("edger analysis works for default params", {
  results <- compareRhythms(countsFromAbundance, exp_design, method = "edger")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic"))
  results <- compareRhythms(counts, exp_design, lengths = lengths, method = "edger")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic"))
  results <- compareRhythms(countsFromAbundance, exp_design_batch, method = "edger")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic"))
  results <- compareRhythms(counts, exp_design_batch, lengths = lengths, method = "edger")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic"))
})

test_that("edger analysis works for different input params", {
  expect_s3_class(compareRhythms(counts, exp_design, lengths = lengths, period = 12, method = "edger"), "data.frame")
  expect_s3_class(compareRhythms(counts, exp_design, lengths = lengths, outliers = TRUE, method = "edger"), "data.frame")
  expect_s3_class(compareRhythms(counts, exp_design, lengths = lengths, rhythm_fdr = 0.01, method = "edger"), "data.frame")
  expect_s3_class(compareRhythms(counts, exp_design, lengths = lengths, compare_fdr = 0.01, method = "edger"), "data.frame")
  expect_s3_class(compareRhythms(counts, exp_design, lengths = lengths, rhythm_fdr = 0.1, compare_fdr = 0.01, method = "edger"), "data.frame")
  expect_s3_class(compareRhythms(counts, exp_design, lengths = lengths, amp_cutoff = 0, method = "edger"), "data.frame")
  expect_error(compareRhythms(counts, exp_design, lengths = lengths, rhythm_fdr = 0, method = "edger"))
  results <- compareRhythms(counts, exp_design, lengths = lengths, just_classify = FALSE, method = "edger")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic", "P66KO_amp",
                 "P66KO_phase", "WT_amp", "WT_phase", "adj_p_val_P66KO_or_WT",
                 "adj_p_val_DR"))
})

test_that("edger analysis works for different input params with batch", {
  expect_s3_class(compareRhythms(counts, exp_design_batch, lengths = lengths, period = 12, method = "edger"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design_batch, period = 12, method = "edger"), "data.frame")
  expect_s3_class(compareRhythms(counts, exp_design_batch, lengths = lengths, rhythm_fdr = 0.01, method = "edger"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design_batch, rhythm_fdr = 0.01, method = "edger"), "data.frame")
  expect_s3_class(compareRhythms(counts, exp_design_batch, lengths = lengths, rhythm_fdr = 0.01, method = "edger"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design_batch, rhythm_fdr = 0.01, method = "edger"), "data.frame")
  expect_s3_class(compareRhythms(counts, exp_design_batch, lengths = lengths, rhythm_fdr = 0.1, compare_fdr = 0.01, method = "edger"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design_batch, rhythm_fdr = 0.1, compare_fdr = 0.01, method = "edger"), "data.frame")
  expect_s3_class(compareRhythms(counts, exp_design_batch, lengths = lengths, amp_cutoff = 0, method = "edger"), "data.frame")
  results <- compareRhythms(counts, exp_design_batch, lengths = lengths, just_classify = FALSE, method = "edger")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic", "P66KO_amp",
                 "P66KO_phase", "WT_amp", "WT_phase", "adj_p_val_P66KO_or_WT",
                 "adj_p_val_DR"))
  results <- compareRhythms(counts, exp_design_batch, lengths = lengths, just_classify = FALSE, method = "edger")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic",
                 "P66KO_amp", "P66KO_phase", "WT_amp", "WT_phase", "adj_p_val_P66KO_or_WT",
                 "adj_p_val_DR"))
  results <- compareRhythms(countsFromAbundance, exp_design_batch, just_classify = FALSE, method = "edger")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic",
                 "P66KO_amp", "P66KO_phase", "WT_amp", "WT_phase", "adj_p_val_P66KO_or_WT",
                 "adj_p_val_DR"))
})
