#' Run differential rhythmicity analysis using linear model comparison
#'
#' This function runs the information criteria-based model selection proposed by
#' Atger et al. (2015) for identifying timeseries with different rhythms in the
#' two datasets.
#'
#' @param y A matrix with gene in the rows and samples in columns containing
#'   data from both datasets
#' @param exp_design A data.frame of the experimental design with at least
#'   columns sample name, time point and group
#' @param period The period of rhythm being tested (default = 24)
#' @param amp_cutoff The minimum peak-to-trough amp in log2 scale considered
#'   biologically relevant (default = 0.5)
#' @param criterion The criterion used for model selection. These can be "aic"
#'   or "bic" (default = "bic")
#' @return A data.frame with symbol, best linear model (termed class) and
#'   estimates of the amplitudes and phases for the two datasets
#'
#' @export

compareRhythms_model_compare <- function(y, exp_design, period = 24,
                                         amp_cutoff = 0.5,
                                         criterion = "bic") {
  input_check(y, exp_design)

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

  model_selection <- limma::selectModel(y, design_list, criterion = criterion)

  model_assignment <- model_selection$pref[model_selection$pref != "noR"]

  model_circ_params <- base::lapply(design_list[-1],
                              compute_model_params, y, group_id)

  circ_params <- base::vapply(model_assignment,
                              function(m) {
                                model_circ_params[[base::as.character(m)]][names(m), ]
                              },
                              FUN.VALUE = double(4))

  results <- data.frame(t(circ_params))
  results <- base::cbind(symbol = names(model_assignment),
                         results,
                         class = unname(model_assignment))

  results$max_amp <- pmax(results[, paste0(group_id[1], "_amp")],
                          results[, paste0(group_id[2], "_amp")])

  results <- results[results$max_amp > amp_cutoff, ]
  results$max_amp <- NULL

  for (i in seq(nrow(results))) {
    if (results[i, "class"] == "DR") {
      if (results[i, paste0(group_id[2], "_amp")] < amp_cutoff) {
        results[i, "class"] == "AR"
        results[i, paste0(group_id[2], "_amp")] <- 0
        results[i, paste0(group_id[2], "_phase")] <- 0
      }

      if (results[i, paste0(group_id[1], "_amp")] < amp_cutoff) {
        results[i, "class"] == "BR"
        results[i, paste0(group_id[1], "_amp")] <- 0
        results[i, paste0(group_id[1], "_phase")] <- 0
      }
    }
  }

  rownames(results) <- NULL

  return(results)
}
