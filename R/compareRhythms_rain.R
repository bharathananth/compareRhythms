#' Run differential rhythmicity analysis defined in Thaben & Westermark
#'
#' @param y A matrix with gene in the rows and samples in columns
#' @param exp_design A data.frame of the experimental design with at least
#'   columns sample name, time point and group
#' @param period The period of rhythm being tested (default = 24)
#' @param rhythm_fdr The false discovery cutoff for finding rhythmic time series
#'   (default = 0.05)
#' @param compare_fdr The false discovery cutoff for the comparison of rhythms
#'   (default = 0.05)
#' @param amp_cutoff The minimum peak-to-trough amp in log2 scale considered
#'   biologically relevant
#' @param include_pvals A flag to select whether p-values, amplitudes and phases
#'   must be included with the results
#' @return A data.frame with the symbol, boolean results of the rhythmicity
#'   tests and (optionally) the p-values and circadian parameters.
#' @export

compareRhythms_rain <- function(y, exp_design, period=24, rhythm_fdr = 0.05,
                                compare_fdr = 0.05, amp_cutoff = 0.5,
                                include_pvals = FALSE) {

  base::stopifnot(
    ncol(y) == nrow(exp_design),
    any(base::colnames(exp_design) == "time"),
    any(base::colnames(exp_design) == "group"),
    all(base::Negate(base::is.na)(y))
  )

  exp_design <- base::cbind(exp_design,
                            col_number = base::seq(base::nrow(exp_design)))

  group_id <- base::unique(exp_design$group)

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


  y_A <- y[, exp_design_A$col_number]
  y_B <- y[, exp_design_B$col_number]

  rain_A <- rain::rain(t(y_A), deltat_A, period,
                       measure.sequence = measure_sequence_A)

  rain_B <- rain::rain(t(y_B), deltat_B, period,
                       measure.sequence = measure_sequence_B)

  rain_results <- data.frame(stats::p.adjust(rain_A[, "pVal"], method = "BH"),
                             stats::p.adjust(rain_B[, "pVal"], method = "BH"))
  colnames(rain_results) <- c("p_val_adjust_A", "p_val_adjust_B")


  circ_params_A <- compute_circ_params(y_A, exp_design_A$time, period = period)

  circ_params_B <- compute_circ_params(y_B, exp_design_B$time, period = period)

  rhythmic_in_A <- (rain_results$p_val_adjust_A < rhythm_fdr) &
    (circ_params_A[, "amps"] > amp_cutoff)

  rhythmic_in_B <- (rain_results$p_val_adjust_B < rhythm_fdr) &
    (circ_params_B[, "amps"] > amp_cutoff)

  rhythmic_in_either <- rhythmic_in_A | rhythmic_in_B

  if (sum(rhythmic_in_either) == 0) {
    stop("Sorry no rhythmic genes in either data set for the thresholds provided.")
  }

  dodr_results <- DODR::robustDODR(t(y_A[rhythmic_in_either, ]),
                                   t(y_B[rhythmic_in_either, ]),
                                   times1 = exp_design_A$time,
                                   times2 = exp_design_B$time,
                                   norm = TRUE,
                                   period = period)
  dodr_results$p_val_adjust = stats::p.adjust(dodr_results$p.value, method = "BH")

  results <- data.frame(symbol = rownames(y_A)[rhythmic_in_either],
                        rhythmic_in_A = rhythmic_in_A[rhythmic_in_either],
                        rhythmic_in_B = rhythmic_in_B[rhythmic_in_either],
                        diff_rhythmic = dodr_results$p_val_adjust < compare_fdr)
  rownames(results) <- NULL

  if (include_pvals) {
    results$p_val_adjust_A = p_val_adjust_A[rhythmic_in_either]
    results$p_val_adjust_B = p_val_adjust_B[rhythmic_in_either]
    results$p_val_adjust_dodr = dodr_results$p_value_adjust
    results$amp_A = circ_params_A[rhythmic_in_either, "amps"]
    results$amp_B = circ_params_B[rhythmic_in_either, "amps"]
    results$phase_A = circ_params_A[rhythmic_in_either, "phases"]
    results$phase_B = circ_params_B[rhythmic_in_either, "phases"]

  }

  return(results)
}
