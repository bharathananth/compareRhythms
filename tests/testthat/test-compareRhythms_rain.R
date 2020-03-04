
context("compareRhythms_rain")

y <- readRDS("y.RDS")
exp_design <- readRDS("exp_design.RDS")

test_that("DODR analysis works for default params", {
  results <- compareRhythms_rain(y, exp_design)
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("symbol", "rhythmic_in_CC", "rhythmic_in_KD", "diff_rhythmic", "class"))
})

test_that("DODR analysis works for different input params", {
  expect_error(compareRhythms_rain(y, exp_design, period = 12))
  expect_s3_class(compareRhythms_rain(y, exp_design, rhythm_fdr = 0.01), "data.frame")
  expect_s3_class(compareRhythms_rain(y, exp_design, compare_fdr = 0.01), "data.frame")
  expect_s3_class(compareRhythms_rain(y, exp_design, rhythm_fdr = 0.1, compare_fdr = 0.01), "data.frame")
  expect_s3_class(compareRhythms_rain(y, exp_design, amp_cutoff = 0), "data.frame")
  expect_error(compareRhythms_rain(y, exp_design, rhythm_fdr = 0))
  results <- compareRhythms_rain(y, exp_design, just_classify = FALSE)
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("symbol", "rhythmic_in_CC", "rhythmic_in_KD", "diff_rhythmic", "class", "CC_amp",
                 "CC_phase", "KD_amp", "KD_phase", "adj_p_val_CC", "adj_p_val_KD", "adj_p_val_dodr"))
})

test_that("RAIN functionality is okay", {
  expect_error(compareRhythms_rain(y[, 1:15], exp_design = exp_design[1:15, ]),
               regexp = "Not enough")
})

test_that("DODR analysis runs for arrhythmic dataset", {
  dim_y <- dim(y)
  y_null <- matrix(0, nrow = dim_y[1], ncol = dim_y[2])
  colnames(y_null) <- colnames(y)
  rownames(y_null) <- rownames(y)
  expect_error(compareRhythms_limma(y_null, exp_design))
})
