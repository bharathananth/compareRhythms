
context("compareRhythms_limma")

load("test_data_ma.rda")

test_that("limma analysis works for default params", {
  results <- compareRhythms(expr, exp_design, method = "limma")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("symbol", "rhythmic_in_CC", "rhythmic_in_KD", "diff_rhythmic", "category"))
})

test_that("limma analysis works for different input params", {
  expect_error(compareRhythms(expr, exp_design, period = 12, method = "limma"))
  expect_s3_class(compareRhythms(expr, exp_design, rhythm_fdr = 0.01, method = "limma"), "data.frame")
  expect_s3_class(compareRhythms(expr, exp_design, compare_fdr = 0.01, method = "limma"), "data.frame")
  expect_s3_class(compareRhythms(expr, exp_design, rhythm_fdr = 0.1, compare_fdr = 0.01, method = "limma"), "data.frame")
  expect_s3_class(compareRhythms(expr, exp_design, amp_cutoff = 0, method = "limma"), "data.frame")
  expect_error(compareRhythms(expr, exp_design, rhythm_fdr = 0, method = "limma"))
  results <- compareRhythms(expr, exp_design, just_classify = FALSE, method = "limma")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("symbol", "rhythmic_in_CC", "rhythmic_in_KD", "diff_rhythmic", "category", "CC_amp",
                 "CC_phase", "KD_amp", "KD_phase", "adj_p_val_CC_or_KD",
                 "adj_p_val_DR"))
})

test_that("limma analysis runs for arrhythmic dataset", {
  dim_y <- dim(expr)
  y_null <- matrix(0, nrow = dim_y[1], ncol = dim_y[2])
  colnames(y_null) <- colnames(expr)
  rownames(y_null) <- rownames(expr)
  expect_error(compareRhythms(y_null, exp_design, method = "limma"))
})
