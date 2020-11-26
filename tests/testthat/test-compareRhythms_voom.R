
context("compareRhythms_voom")

load("test_data_rnaseq.rda")

exp_design_batch <- cbind(exp_design, batch, stringsAsFactors=TRUE)

test_that("limma-voom analysis works for default params", {
  results <- compareRhythms(countsFromAbundance, exp_design, method = "voom")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic"))
  results <- compareRhythms(countsFromAbundance, exp_design_batch, method = "voom")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic"))
})

test_that("limma-voom analysis works for different input params", {
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design, period = 12, method = "voom"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design, outliers = TRUE, method = "voom"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design, rhythm_fdr = 0.01, method = "voom"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design, compare_fdr = 0.01, method = "voom"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design, rhythm_fdr = 0.1, compare_fdr = 0.01, method = "voom"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design, amp_cutoff = 0, method = "voom"), "data.frame")
  expect_error(compareRhythms(countsFromAbundance, exp_design, rhythm_fdr = 0, method = "voom"))
  results <- compareRhythms(countsFromAbundance, exp_design, just_classify = FALSE, method = "voom")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic",
                 "P66KO_amp", "P66KO_phase", "WT_amp", "WT_phase", "adj_p_val_P66KO_or_WT",
                 "adj_p_val_DR"))
})

test_that("limma-voom analysis works for different input params with batch", {
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design_batch, period = 12, method = "voom"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design_batch, outliers = TRUE, method = "voom"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design_batch, rhythm_fdr = 0.01, method = "voom"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design_batch, compare_fdr = 0.01, method = "voom"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design_batch, rhythm_fdr = 0.1, compare_fdr = 0.01, method = "voom"), "data.frame")
  expect_s3_class(compareRhythms(countsFromAbundance, exp_design_batch, amp_cutoff = 0, method = "voom"), "data.frame")
  expect_error(compareRhythms(countsFromAbundance, exp_design_batch, rhythm_fdr = 0, method = "voom"))
  results <- compareRhythms(countsFromAbundance, exp_design_batch, just_classify = FALSE, method = "voom")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_P66KO", "rhythmic_in_WT", "diff_rhythmic", "P66KO_amp",
                 "P66KO_phase", "WT_amp", "WT_phase", "adj_p_val_P66KO_or_WT",
                 "adj_p_val_DR"))
})
