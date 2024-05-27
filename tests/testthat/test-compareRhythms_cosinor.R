context("compareRhythms_cosinor")

load("test_data_ma.rda")

test_that("cosinor analysis works for default params and independent sampling", {
  results <- compareRhythms(expr, exp_design, method = "cosinor")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_CC", "rhythmic_in_KD", "diff_rhythmic"))
})

test_that("cosinor analysis works for different input params", {
  expect_error(compareRhythms(expr, exp_design, period = 12, method = "cosinor"))
  expect_s3_class(compareRhythms(expr, exp_design, rhythm_fdr = 0.01, method = "cosinor"), "data.frame")
  expect_s3_class(compareRhythms(expr, exp_design, compare_fdr = 0.01, method = "cosinor"), "data.frame")
  expect_s3_class(compareRhythms(expr, exp_design, rhythm_fdr = 0.1, compare_fdr = 0.01, method = "cosinor"), "data.frame")
  expect_s3_class(compareRhythms(expr, exp_design, amp_cutoff = 0, method = "cosinor"), "data.frame")
  expect_error(compareRhythms(expr, exp_design, rhythm_fdr = 0, method = "cosinor"))
  results <- compareRhythms(expr, exp_design, just_classify = FALSE, method = "cosinor")
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_CC", "rhythmic_in_KD", "diff_rhythmic", "CC_amp",
                 "CC_phase", "KD_amp", "KD_phase", "adj_p_val_CC_or_KD", "adj_p_val_DR"))
})

test_that("cosinor analysis runs for arrhythmic dataset", {
  dim_y <- dim(expr)
  y_null <- matrix(0, nrow = dim_y[1], ncol = dim_y[2])
  colnames(y_null) <- colnames(expr)
  rownames(y_null) <- rownames(expr)
  expect_error(compareRhythms(y_null, exp_design, method = "cosinor"))
})

exp_design <- dplyr::mutate(dplyr::group_by(exp_design, group, time), ID = dplyr::row_number())
exp_design$ID <- factor(exp_design$ID)

test_that("longitudinal cosinor analysis works for default params and independent sampling", {
  results <- compareRhythms(expr, exp_design, method = "cosinor", longitudinal = TRUE)
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_CC", "rhythmic_in_KD", "diff_rhythmic"))
})

test_that("cosinor analysis works for different input params", {
  expect_error(compareRhythms(expr, exp_design, period = 12, method = "cosinor", longitudinal = TRUE))
  expect_s3_class(compareRhythms(expr, exp_design, rhythm_fdr = 0.01, method = "cosinor", longitudinal = TRUE), "data.frame")
  expect_s3_class(compareRhythms(expr, exp_design, compare_fdr = 0.01, method = "cosinor", longitudinal = TRUE), "data.frame")
  expect_s3_class(compareRhythms(expr, exp_design, rhythm_fdr = 0.1, compare_fdr = 0.01, method = "cosinor", longitudinal = TRUE), "data.frame")
  expect_s3_class(compareRhythms(expr, exp_design, amp_cutoff = 0, method = "cosinor", longitudinal = TRUE), "data.frame")
  expect_error(compareRhythms(expr, exp_design, rhythm_fdr = 0, method = "cosinor", longitudinal = TRUE))
  results <- compareRhythms(expr, exp_design, just_classify = FALSE, method = "cosinor", longitudinal = TRUE)
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("id", "category", "rhythmic_in_CC", "rhythmic_in_KD", "diff_rhythmic", "CC_amp",
                 "CC_phase", "KD_amp", "KD_phase", "adj_p_val_CC_or_KD", "adj_p_val_DR"))
})

test_that("cosinor analysis runs for arrhythmic dataset", {
  dim_y <- dim(expr)
  y_null <- matrix(0, nrow = dim_y[1], ncol = dim_y[2])
  colnames(y_null) <- colnames(expr)
  rownames(y_null) <- rownames(expr)
  expect_error(compareRhythms(y_null, exp_design, method = "cosinor", longitudinal = TRUE))
})


