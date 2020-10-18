context("compareRhythms_model_select")

load("test_data_ma.rda")

test_that("model selection works for default params", {
  results <- compareRhythms(expr, exp_design, method = "mod_sel")
  expect_s3_class(results, "data.frame")
  expect_named(results, c("symbol", "category"))
})

test_that("model selection works for different input params", {
  expect_s3_class(compareRhythms(expr, exp_design, period = 12, method = "mod_sel"), "data.frame")
  expect_s3_class(compareRhythms(expr, exp_design, amp_cutoff = 0, method = "mod_sel"), "data.frame")
  expect_s3_class(compareRhythms(expr, exp_design, criterion = "aic", method = "mod_sel"), "data.frame")
  results <- compareRhythms(expr, exp_design, just_classify = FALSE, method = "mod_sel")
  expect_s3_class(results, "data.frame")
  expect_named(results, c("symbol", "category", "CC_amp", "CC_phase", "KD_amp", "KD_phase", "weights"))
})

test_that("model selection runs for arrhythmic dataset", {
  dim_y <- dim(expr)
  y_null <- matrix(0, nrow = dim_y[1], ncol = dim_y[2])
  colnames(y_null) <- colnames(expr)
  rownames(y_null) <- rownames(expr)
  expect_error(compareRhythms(y_null, exp_design, method = "mod_sel"))
})
