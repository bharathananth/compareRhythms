context("compareRhythms_model_select")

y <- readRDS("y.RDS")
exp_design <- readRDS("exp_design.RDS")

test_that("model selection works for default params", {
  results <- compareRhythms_model_select(y, exp_design)
  expect_s3_class(results, "data.frame")
  expect_named(results,
               c("symbol", "class", "CC_amp", "CC_phase", "KD_amp", "KD_phase"))
})

test_that("model selection works for different input params", {
  expect_s3_class(compareRhythms_model_select(y, exp_design, period = 12), "data.frame")
  expect_s3_class(compareRhythms_model_select(y, exp_design, amp_cutoff = 0), "data.frame")
  expect_s3_class(compareRhythms_model_select(y, exp_design, criterion = "aic"), "data.frame")
})

test_that("model selection runs for arrhythmic dataset", {
  dim_y <- dim(y)
  y_null <- matrix(0, nrow = dim_y[1], ncol = dim_y[2])
  colnames(y_null) <- colnames(y)
  rownames(y_null) <- rownames(y)
  expect_error(compareRhythms_model_select(y_null, exp_design))
})
