#' Run differential rhythmicity analysis using linear model selection
#'
#' This function runs the information criteria-based model selection proposed by
#' Atger et al. (2015) for identifying timeseries with different rhythms in the
#' two datasets.
#'
#' @param expr A matrix with gene in the rows and samples in columns containing
#'   data from both datasets
#' @param exp_design A data.frame of the experimental design with at least
#'   columns sample name, time point and group
#' @param period The period of rhythm being tested (default = 24)
#' @param amp_cutoff The minimum peak-to-trough amp in log2 scale considered
#'   biologically relevant (default = 0.5)
#' @param criterion The criterion used for model selection. These can be "aic"
#'   or "bic" (default = "bic")
#' @param schwartz_wt_cutoff The conditional probability that the best models is
#'   the true model. Genes with a smaller conditional probability smaller than
#'   this cutoff are deemed unclassifiable (default = 0.4).
#' @param just_classify Logical to select whether p-values, amplitudes and
#'   phases must be supressed in the results
#' @return A data.frame with symbol, best linear model (termed category) and
#'   estimates of the amplitudes and phases for the two datasets

compareRhythms_model_select <- function(expr, exp_design, period = 24,
                                         amp_cutoff = 0.5,
                                         criterion = "bic",
                                         schwartz_wt_cutoff = 0.4,
                                         just_classify = TRUE) {

  group_id <- base::unique(exp_design$group)

  exp_design <- base::cbind(exp_design,
                            inphase = cos(2 * pi * exp_design$time / period),
                            outphase = sin(2 * pi * exp_design$time / period))

  design <- stats::model.matrix(~0 + group + group:inphase + group:outphase,
                                data = exp_design)
  colnames(design) <- gsub("group", "", colnames(design))
  colnames(design) <- gsub(":", "_", colnames(design))

  design_DR <- design

  design_AR <- design[, !grepl(paste0(group_id[2], "_"),
                               colnames(design))]

  design_BR <- design[, !grepl(paste0(group_id[1], "_"),
                               colnames(design))]

  design_ABR <- stats::model.matrix(~0 + group + inphase + outphase,
                                    data = exp_design)

  colnames(design_ABR) <- base::gsub("group", "", colnames(design_ABR))

  design_noR <- stats::model.matrix(~0 + group, data = exp_design)

  colnames(design_noR) <- gsub("group", "", colnames(design_noR))

  design_list <- list(noR = design_noR,
                      AR = design_AR,
                      BR = design_BR,
                      ABR = design_ABR,
                      DR = design_DR)

  model_selection <- limma::selectModel(expr, design_list, criterion = criterion)

  model_post_prob <- base::apply(model_selection$IC, 1,
                                    function(x) base::max(exp(-0.5*x)/sum(exp(-0.5*x))))

  model_assignment <- model_selection$pref[model_post_prob >= schwartz_wt_cutoff]

  assertthat::assert_that(assertthat::not_empty(model_assignment),
                          msg = "Sorry no rhythmic genes in either dataset for the thresholds provided.")

  model_circ_params <- base::lapply(design_list[-1],
                              compute_model_params, expr, group_id)

  model_circ_params[["noR"]] <- matrix(0, nrow = nrow(expr), ncol = 4,
                                       dimnames = dimnames(model_circ_params[["ABR"]]))

  circ_params <- base::vapply(model_assignment,
                              function(m) {
                                model_circ_params[[base::as.character(m)]][names(m), ]
                              },
                              FUN.VALUE = double(4))

  results <- data.frame(t(circ_params))
  results <- base::cbind(symbol = names(model_assignment),
                         results,
                         category = unname(model_assignment))

  results$max_amp <- pmax(results[, paste0(group_id[1], "_amp")],
                          results[, paste0(group_id[2], "_amp")])

  for (i in seq(nrow(results))) {
    if (results$max_amp[i] < amp_cutoff) {
      results[i, "category"] <- "noR"
      results[i, paste0(group_id[1], "_amp")] <- 0
      results[i, paste0(group_id[1], "_phase")] <- 0
      results[i, paste0(group_id[2], "_amp")] <- 0
      results[i, paste0(group_id[2], "_phase")] <- 0
    } else if (results[i, "category"] == "DR") {
      if (results[i, paste0(group_id[2], "_amp")] < amp_cutoff) {
        results[i, "category"] == "AR"
        results[i, paste0(group_id[2], "_amp")] <- 0
        results[i, paste0(group_id[2], "_phase")] <- 0
      }

      if (results[i, paste0(group_id[1], "_amp")] < amp_cutoff) {
        results[i, "category"] == "BR"
        results[i, paste0(group_id[1], "_amp")] <- 0
        results[i, paste0(group_id[1], "_phase")] <- 0
      }
    }
  }

  results$max_amp <- NULL
  rownames(results) <- NULL
  colnames(results) <- gsub("A", group_id[1], colnames(results))
  colnames(results) <- gsub("B", group_id[2], colnames(results))

  main_cols <- c("symbol", "category")
  if (just_classify) {
    results <- results[, main_cols]
  } else {
    results <- results[, c(main_cols,
                           base::setdiff(colnames(results), main_cols))]
  }

  return(results)
}
