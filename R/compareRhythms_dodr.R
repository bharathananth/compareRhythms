#' Run differential rhythmicity analysis defined in Thaben & Westermark
#'
#' @param expr A matrix of expression values with gene in the rows and samples in columns
#' @inheritParams compareRhythms
#' @keywords internal

compareRhythms_dodr <- function(expr, exp_design, period=24, rhythm_fdr = 0.05,
                                compare_fdr = 0.05, amp_cutoff = 0.5,
                                just_classify = TRUE) {

  exp_design <- base::cbind(exp_design,
                            col_number = base::seq(base::nrow(exp_design)))

  group_id <- base::levels(exp_design$group)
  assertthat::are_equal(length(group_id), 2)

  exp_design_A <- exp_design[exp_design$group == group_id[1], ]

  exp_design_A <- exp_design_A[base::order(exp_design_A$time), ]

  exp_design_B <- exp_design[exp_design$group == group_id[2], ]

  exp_design_B <- exp_design_B[base::order(exp_design_B$time), ]

  deltat_A <- min(diff(base::unique(exp_design_A$time)))

  time_A <- base::seq(min(exp_design_A$time),
                      max(exp_design_A$time),
                      by = deltat_A)

  unique_times_A <- base::table(exp_design_A$time)

  measure_sequence_A <- base::vapply(time_A,
                                     function(t) {
                                       ifelse(any(names(unique_times_A) == t),
                                              unique_times_A[names(unique_times_A) == t],
                                              0)
                                     },
                                     integer(1))


  deltat_B <- min(diff(base::unique(exp_design_B$time)))

  time_B <- base::seq(min(exp_design_B$time),
                      max(exp_design_B$time),
                      by = deltat_B)

  unique_times_B <- base::table(exp_design_B$time)

  measure_sequence_B <- base::vapply(time_B,
                                     function(t) {
                                       ifelse(any(names(unique_times_B) == t),
                                              unique_times_B[names(unique_times_B) == t],
                                              0)
                                     },
                                     integer(1))


  expr_A <- expr[, exp_design_A$col_number]
  expr_B <- expr[, exp_design_B$col_number]

  assertthat::assert_that((sum(measure_sequence_A) >= 12) && (sum(measure_sequence_B) >= 12),
                          msg = "Not enough samples to run RAIN rhythmicity analysis confidently.")

  rain_A <- rain::rain(t(expr_A), deltat_A, period,
                       measure.sequence = measure_sequence_A)

  rain_B <- rain::rain(t(expr_B), deltat_B, period,
                       measure.sequence = measure_sequence_B)

  rain_results <- data.frame(stats::p.adjust(rain_A[, "pVal"], method = "BH"),
                             stats::p.adjust(rain_B[, "pVal"], method = "BH"))
  colnames(rain_results) <- c("adj_p_val_A", "adj_p_val_B")


  circ_params_A <- compute_circ_params(expr_A, exp_design_A$time, period = period)

  circ_params_B <- compute_circ_params(expr_B, exp_design_B$time, period = period)

  rhythmic_in_A <- (rain_results$adj_p_val_A < rhythm_fdr) &
    (circ_params_A[, "amps"] > amp_cutoff)

  rhythmic_in_B <- (rain_results$adj_p_val_B < rhythm_fdr) &
    (circ_params_B[, "amps"] > amp_cutoff)

  rhythmic_in_either <- rhythmic_in_A | rhythmic_in_B

  assertthat::assert_that(sum(rhythmic_in_either) > 0,
                          msg = "Sorry no rhythmic genes in either dataset for the thresholds provided.")

  dodr_results <- DODR::robustDODR(t(expr_A[rhythmic_in_either, ]),
                                   t(expr_B[rhythmic_in_either, ]),
                                   times1 = exp_design_A$time,
                                   times2 = exp_design_B$time,
                                   norm = TRUE,
                                   period = period)
  dodr_results$adj_p_val <- stats::p.adjust(dodr_results$p.value, method = "BH")

  results <- data.frame(symbol = rownames(expr_A)[rhythmic_in_either],
                        rhythmic_in_A = rhythmic_in_A[rhythmic_in_either],
                        rhythmic_in_B = rhythmic_in_B[rhythmic_in_either],
                        diff_rhythmic = dodr_results$adj_p_val < compare_fdr,
                        stringsAsFactors = FALSE)
  rownames(results) <- NULL

  results$category <- base::mapply(categorize,
                                results$rhythmic_in_A,
                                results$rhythmic_in_B,
                                results$diff_rhythmic)

  if (!just_classify) {
    expand_results <- data.frame(
      A_amp = circ_params_A[rhythmic_in_either, "amps"],
      A_phase = circ_params_A[rhythmic_in_either, "phases"],
      B_amp = circ_params_B[rhythmic_in_either, "amps"],
      B_phase = circ_params_B[rhythmic_in_either, "phases"],
      adj_p_val_A = rain_results$adj_p_val_A[rhythmic_in_either],
      adj_p_val_B = rain_results$adj_p_val_B[rhythmic_in_either],
      adj_p_val_dodr = dodr_results$adj_p_val
    )
    results <- base::cbind(results, expand_results)
  }

  colnames(results) <- gsub("A", group_id[1], colnames(results))
  colnames(results) <- gsub("B", group_id[2], colnames(results))

  return(results)
}
