#' Run differential rhythmicity analysis using linear model comparison
#'
#' This function runs the information criteria-based model selection proposed by
#' Atger et al. (2015) for identifying timeseries with different rhythms in the two
#' datasets.
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
#' @return A data.frame with symbol, best linear model and estimates of the
#'   amplitudes and phases for the two datasets
#'
#' @export

compareRhythms_model_compare <- function(y, exp_design, period = 24,
                                         amp_cutoff = 0.5,
                                         criterion = "bic") {
  base::stopifnot(
    ncol(y) == nrow(exp_design),
    any(base::colnames(exp_design) == "time"),
    any(base::colnames(exp_design) == "group"),
    all(base::Negate(base::is.na)(y))
  )

  group_id <- base::unique(exp_design$group)

  exp_design <- exp_design %>%
    base::transform(inphase = cos(2 * pi * time / period),
                    outphase = sin(2 * pi * time / period))

  design <- stats::model.matrix(~0 + group + group:inphase + group:outphase,
                                data = exp_design) %>%
    magrittr::set_colnames(gsub("group", "", colnames(.))) %>%
    magrittr::set_colnames(gsub(":", "_", colnames(.)))

  colnames(design) <- gsub("group", "", colnames(design)) %>% {
    gsub(":", "_", .)
  }

  design_DR <- design

  design_AR <- design[, !grepl(paste0(group_id[2], "_"),
                               colnames(design))]

  design_BR <- design[, !grepl(paste0(group_id[1], "_"),
                               colnames(design))]

  design_ABR <- stats::model.matrix(~0 + group + inphase + outphase,
                                    data = exp_design) %>%
    magrittr::set_colnames(gsub("group", "", colnames(.)))

  design_noR <- stats::model.matrix(~0 + group, data = exp_design) %>%
    magrittr::set_colnames(gsub("group", "", colnames(.)))

  design_list <- list(noR = design_noR,
                      AR = design_AR,
                      BR = design_BR,
                      ABR = design_ABR,
                      DR = design_DR)

  model_selection <- limma::selectModel(y, design_list, criterion = criterion)

  model_assignment <- with(model_selection, pref[pref != "noR"])

  model_circ_params <- lapply(design_list[-1],
                              compute_model_params, y, group_id)

  circ_params <- base::vapply(model_assignment,
                              function(m) {
                                model_circ_params[[base::as.character(m)]][names(m), ]
                              },
                              FUN.VALUE = double(4)) %>% t()

  results <- data.frame(circ_params) %>%
             {base::cbind(symbol = names(model_assignment),
                         best_model = unname(model_assignment), .,
                         max_amp = pmax(.[, paste0(group_id[1], "_amp")],
                                        .[, paste0(group_id[2], "_amp")]))} %>%
             base::subset(max_amp > amp_cutoff) %>%
             base::subset(select = -max_amp)

  for (i in seq(nrow(results))) {
    if (results[i, "best_model"] == "DR") {
      if (results[i, paste0(group_id[2], "_amp")] < amp_cutoff) {
        results[i, "best_model"] == "AR"
        results[i, paste0(group_id[2], "_amp")] <- 0
        results[i, paste0(group_id[2], "_phase")] <- 0
      }

      if (results[i, paste0(group_id[1], "_amp")] < amp_cutoff) {
        results[i, "best_model"] == "BR"
        results[i, paste0(group_id[1], "_amp")] <- 0
        results[i, paste0(group_id[1], "_phase")] <- 0
      }
    }
  }

  return(results)
}

compute_model_params <- function(d, y, group_id) {
  fit <- limma::lmFit(y, d)
  coeffs <- fit$coefficients
  if (any(base::grepl(paste0(group_id[1], "_"),
                      colnames(coeffs)))) {
    rhy_params <- coeffs[, base::paste(group_id[1],
                                       c("inphase", "outphase"),
                                       sep = "_")]
    amps_A <- 2 * sqrt(base::rowSums(rhy_params^2))
    phases_A <- base::atan2(rhy_params[, 2], rhy_params[, 1])
  } else {
    amps_A <- 0
    phases_A <- 0
  }

  if (any(base::grepl(paste0(group_id[2], "_"),
                      colnames(coeffs)))) {
    rhy_params <- coeffs[, base::paste(group_id[2],
                                       c("inphase", "outphase"),
                                       sep = "_")]
    amps_B <- 2 * sqrt(base::rowSums(rhy_params^2))
    phases_B <- base::atan2(rhy_params[, 2], rhy_params[, 1])
  } else {
    amps_B <- 0
    phases_B <- 0
  }

  if (all(base::is.element(c("inphase", "outphase"),
                           colnames(coeffs)))) {
    rhy_params <- coeffs[, c("inphase", "outphase")]
    amps <- 2 * sqrt(base::rowSums(rhy_params^2))
    phases <- base::atan2(rhy_params[, 2], rhy_params[, 1])
    model_params <- base::cbind(amps, phases, amps, phases)
  } else {
    model_params <- base::cbind(amps_A, phases_A, amps_B, phases_B)
  }

  colnames(model_params) <- base::paste(rep(group_id, each = 2),
                                        c("amp", "phase"), sep = "_")

  return(model_params)
}
