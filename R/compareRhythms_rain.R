#' Run differential rhythmicity analysis defined in Thaben & Westermark
#'
#' @param y A matrix with gene in the rows and samples in columns
#' @param exp_design A data.frame of the experimental design with at least
#'   columns sample name, time point and group
#' @param period The period of rhythm being tested (default = 24)
#' @param rhythm_fdr The false discovery cutoff for finding rhythmic time series
#'   (default = 0.05)
#' @param amp_cutoff The minimum peak-to-trough amp in log2 scale
#'   considered biologically relevant
#' @param include_pvals
#' @export

compareRhythms_rain <- function(y, exp_design, period=24, rhythm_fdr = 0.05,
                                compare_fdr = 0.05, amp_cutoff = 0.5,
                                include_pvals = FALSE) {

  exp_design <- base::cbind(exp_design,
                            col_number = base::seq(base::nrow(exp_design)))
  group_ID <- base::unique(exp_design$group)

  exp_design_A <- base::subset(exp_design, group == group_ID[1])

  exp_design_A <- exp_design_A[base::order(exp_design_A$time), ]

  exp_design_B <- base::subset(exp_design, group == group_ID[2])

  exp_design_B <- exp_design_B[base::order(exp_design_B$time), ]

  deltat_A <- exp_design_A %$% base::unique(time) %>%
                  diff %>%
                  min


  time_A <- base::seq(min(exp_design_A$time),
                      max(exp_design_A$time),
                      by = deltat_A)

  measure_sequence_A <- base::table(exp_design_A$time) %>% {
                        base::sapply(time_A,
                                     function(t) ifelse(any(names(.) == t),
                                             .[names(.) == t],
                                             0))
                          }

  deltat_B <- exp_design_B %$% base::unique(time) %>%
              diff %>%
              min


  time_B <- base::seq(min(exp_design_B$time),
                      max(exp_design_B$time),
                      by = deltat_B)

  measure_sequence_B <- base::table(exp_design_B$time) %>% {
                        base::sapply(time_B,
                                     function(t) ifelse(any(names(.) == t),
                                                   .[names(.) == t], 0))
                        }


  y_A <- y[, exp_design_A$col_number]
  y_B <- y[, exp_design_B$col_number]

  rain_A <- rain::rain(t(y_A), deltat_A, period,
                       measure.sequence = measure_sequence_A)

  rain_B <- rain::rain(t(y_B), deltat_B, period,
                       measure.sequence = measure_sequence_B)

  rain_results <- data.frame(p_val_adjust_A = p.adjust(rain_A[, "pVal"],
                                                       method = "BH"),
                             p_val_adjust_B = p.adjust(rain_B[, "pVal"],
                                                       method = "BH"))


  circ_params_A <- compute_circ_params(y_A, time_A, period = period)

  circ_params_B <- compute_circ_params(y_B, time_B, period = period)

  rhythmic_in_A <- (rain_results$p_val_adjust_A < rhythm_fdr) &
                    (circ_params_A[, "amps"] > amp_cutoff)

  rhythmic_in_B <- (rain_results$p_val_adjust_B < rhythm_fdr) &
                    (circ_params_B[, "amps"] > amp_cutoff)

  rhythmic_in_either <- rhythmic_in_A | rhythmic_in_B

  if (sum(rhythmic_in_either)==0) {
    stop("Sorry no rhythmic genes in either data set for the thresholds provided.")
  }

  dodr_results <- DODR::dodr(t(y_A[rhythmic_in_either, ]),
                             t(y_B[rhythmic_in_either, ]),
                             times1 = time_A,
                             times2 = time_B,
                             norm = FALSE,
                             period = period) %>%
                  transform(p_value_adjust = p.adjust(p.value, method="BH"))

  results <- data.frame(symbol = rownames(y_A)[rhythmic_in_either],
                        rhythmic_in_A = rhythmic_in_A[rhythmic_in_either],
                        rhythmic_in_B = rhythmic_in_B[rhythmic_in_either],
                        diff_rhythmic = dodr_results$p_value_adjust < compare_fdr)

  if (include_pvals) {
    results <- transform(results,
                         p_val_adjust_A = p_val_adjust_A[rhythmic_in_either],
                         p_val_adjust_B = p_val_adjust_B[rhythmic_in_either],
                         p_val_adjust_dodr = dodr_results$p_value_adjust,
                         amp_A = circ_params_A[rhythmic_in_either, "amps"],
                         amp_B = circ_params_B[rhythmic_in_either, "amps"],
                         phase_A = circ_params_A[rhythmic_in_either, "phases"],
                         phase_B = circ_params_B[rhythmic_in_either, "phases"])

  }

 return(results)
}
