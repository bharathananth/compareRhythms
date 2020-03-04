
context("compareRhythms_rain")

y <- readRDS("y.RDS")
exp_design <- readRDS("exp_design.RDS")

test_that("DODR analysis works for default params", {
  results <- compareRhythms_rain(y, exp_design)
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("symbol", "rhythmic_in_A", "rhythmic_in_B", "diff_rhythmic", "class"))
})

test_that("DODR analysis works for different input params", {
  expect_error(compareRhythms_rain(y, exp_design, period=12))
  expect_s3_class(compareRhythms_rain(y, exp_design, rhythm_fdr = 0.01), "data.frame")
  expect_s3_class(compareRhythms_rain(y, exp_design, compare_fdr = 0.01), "data.frame")
  expect_s3_class(compareRhythms_rain(y, exp_design, rhythm_fdr = 0.1, compare_fdr = 0.01), "data.frame")
  expect_s3_class(compareRhythms_rain(y, exp_design, amp_cutoff = 0), "data.frame")
  expect_error(compareRhythms_rain(y, exp_design, rhythm_fdr = 0))
  results <- compareRhythms_rain(y, exp_design, just_classify = FALSE)
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("symbol", "rhythmic_in_A", "rhythmic_in_B", "diff_rhythmic", "class", "CC_amp",
                 "CC_phase", "KD_amp", "KD_phase", "adj_p_val_A_or_B",
                 "adj_p_val_DR"))
})
