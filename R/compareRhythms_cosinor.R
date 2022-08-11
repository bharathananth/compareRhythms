#' Run differential rhythmicity analysis for microarray using limma
#'
#' @param eset A matrix of expression values with gene in the rows and samples in columns
#' @inheritParams compareRhythms
#' @keywords internal
#' @export

compareRhythms_cosinor <- function(data, exp_design, period, rhythm_fdr,
                                 compare_fdr, amp_cutoff, just_classify, longitudinal) {

  group_id <- base::levels(exp_design$group)

  exp_design <- base::cbind(exp_design,
                            inphase = cos(2 * pi * exp_design$time / period),
                            outphase = sin(2 * pi * exp_design$time / period))

  lmer_control <- lme4::lmerControl(check.conv.singular = lme4::.makeCC(action = "ignore", tol = formals(lme4::isSingular)$tol))

  if ("batch" %in% colnames(exp_design)) {

    if (longitudinal) {
      fit <- lapply(1:nrow(data),
                    function(i) list(lme4::lmer(data[i,]~(1|ID) + group + group:inphase + group:outphase + batch, data = exp_design, REML=FALSE, control = lmer_control),
                                     lme4::lmer(data[i,]~(1|ID) + group + inphase + outphase + batch, data = exp_design, REML=FALSE, control = lmer_control),
                                     lme4::lmer(data[i,]~(1|ID) + group + batch, data = exp_design, REML=FALSE, control = lmer_control)))
    } else {
      fit <- lapply(1:nrow(data),
                    function(i) list(lm(data[i,]~0 + group + group:inphase + group:outphase + batch, data = exp_design),
                                     lm(data[i,]~0 + group + inphase + outphase + batch, data = exp_design),
                                     lm(data[i,]~0 + group + batch, data = exp_design)))
    }

    fit_coeffs <- vapply(fit, function(f){
      coefficients <- if (longitudinal) lme4::fixef(f[[1]]) else coef(f[[1]])
      names(coefficients) <- gsub("group", "", names(coefficients))
      names(coefficients) <- gsub(":", "_", names(coefficients))
      return(coefficients)
    }, FUN.VALUE = double(7L))
  } else {

   if (longitudinal) {
     fit <- lapply(1:nrow(data),
                   function(i) list(lme4::lmer(data[i,]~(1|ID) + group + group:inphase + group:outphase, data = exp_design, REML=FALSE, control = lmer_control, na.action = na.omit),
                                    lme4::lmer(data[i,]~(1|ID) + group + inphase + outphase, data = exp_design, REML=FALSE, control = lmer_control, na.action = na.omit),
                                    lme4::lmer(data[i,]~(1|ID) + group, data = exp_design, REML=FALSE, control = lmer_control, na.action = na.omit)))
   } else {
     fit <- lapply(1:nrow(data),
                   function(i) list(lm(data[i,]~0 + group + group:inphase + group:outphase, data = exp_design, na.action = na.omit),
                                    lm(data[i,]~0 + group + inphase + outphase, data = exp_design, na.action = na.omit),
                                    lm(data[i,]~0 + group, data = exp_design, na.action = na.omit)))
   }

    fit_coeffs <- vapply(fit, function(f){
      coefficients <- if (longitudinal) lme4::fixef(f[[1]]) else coef(f[[1]])
      names(coefficients) <- gsub("group", "", names(coefficients))
      names(coefficients) <- gsub(":", "_", names(coefficients))
      return(coefficients)
    }, FUN.VALUE = double(6L))

  }

  fit_coeffs <- t(fit_coeffs)
  rownames(fit_coeffs) <- rownames(data)

  rhythmic_in_either <- vapply(fit, function(f) {
                                  d <- anova(f[[3]], f[[1]], test=ifelse(longitudinal,"LRT", "F"))
                                  ifelse(longitudinal, d$`Pr(>Chisq)`[2], d$`Pr(>F)`[2])
                               }, FUN.VALUE = double(1L))

  names(rhythmic_in_either) <- rownames(data)

  adj_pval <- p.adjust(rhythmic_in_either, method = "fdr")

  results <- compute_model_params(fit_coeffs, group_id, type = "coef")

  results <- data.frame(results)

  results$id <- rownames(results)

  rownames(results) <- NULL

  results$max_amp <- pmax(results[, paste0(group_id[1], "_amp")],
                          results[, paste0(group_id[2], "_amp")])

  results$adj_p_val_A_or_B <- adj_pval

  results <- results[(results$adj_p_val_A_or_B < rhythm_fdr) &
                      (results$max_amp >= amp_cutoff), ]

  assertthat::assert_that(assertthat::not_empty(results),
                          msg = "Sorry no rhythmic genes in either dataset for the thresholds provided.")

  results$max_amp <- NULL

  diff_rhy_results <- vapply(fit,
                             function(f) {
                               d <- anova(f[[2]], f[[1]], test=ifelse(longitudinal,"LRT", "F"))
                               ifelse(longitudinal, d$`Pr(>Chisq)`[2], d$`Pr(>F)`[2])
                             }, FUN.VALUE = double(1L))
  names(diff_rhy_results) <- rownames(data)

  diff_rhy_results <- diff_rhy_results[results$id]

  results$adj_p_val_DR <- stats::p.adjust(diff_rhy_results,
                                          method = "BH")
  results$diff_rhythmic <- results$adj_p_val_DR < compare_fdr

  results$rhythmic_in_A <- results[, paste0(group_id[1], "_amp")] > amp_cutoff

  results$rhythmic_in_B <- results[, paste0(group_id[2], "_amp")] > amp_cutoff

  results$category <- base::mapply(categorize,
                                results$rhythmic_in_A,
                                results$rhythmic_in_B,
                                results$diff_rhythmic)


  main_cols <- c("id", "category", "rhythmic_in_A", "rhythmic_in_B",
                 "diff_rhythmic")

  results <- results[, c(main_cols,
                         base::setdiff(colnames(results), main_cols))]

  if (just_classify) {
    results <- results[, main_cols]
  }

  rownames(results) <- NULL
  colnames(results) <- gsub("A", group_id[1], colnames(results))
  colnames(results) <- gsub("B", group_id[2], colnames(results))

  return(results)
}
