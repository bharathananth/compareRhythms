#' Run differential rhythmicity analysis using linear model comparison defined
#' in Atger et al.
#'
#' @export

compareRhythms_model_compare <- function(y, exp_design, period=24,
                                         amp_cutoff = 0.5,
                                         criterion = "bic") {
  base::stopifnot(
    ncol(y) == nrow(exp_design),
    any(base::colnames(exp_design) == "time"),
    any(base::colnames(exp_design) == "group"),
    all(base::Negate(base::is.na)(y))
  )

  group_ID <- base::unique(exp_design$group)

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

  design_AR <- design[, !grepl(paste0(group_ID[2], "_"),
                               colnames(design))]

  design_BR <- design[, !grepl(paste0(group_ID[1], "_"),
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

  model_circ_params <- lapply(design_list[-1],
                              compute_model_params, y, group_ID)

  circ_params <- base::vapply(model_selection$pref[model_selection$pref != "noR"],
                              function(i) model_circ_params[[i]][names(i),],
                              double(4))

  return(circ_params)

}

compute_model_params <- function(d, y, group_ID) {
  fit <- limma::lmFit(y, d)
  coeffs <- fit$coefficients
  if (any(base::grepl(paste0(group_ID[1], "_"),
                  colnames(coeffs)))) {
    rhy_params <- coeffs[, base::paste(group_ID[1],
                                       c("inphase", "outphase"),
                                       sep="_")]
    amps_A <- 2 * sqrt(base::rowSums(rhy_params^2))
    phases_A <- base::atan2(rhy_params[2], rhy_params[1])
  } else {
    amps_A <- 0
    phases_A <- 0
  }

  if (any(base::grepl(paste0(group_ID[2], "_"),
                      colnames(coeffs)))) {
    rhy_params <- coeffs[, base::paste(group_ID[2],
                                       c("inphase", "outphase"),
                                       sep="_")]
    amps_B <- 2 * sqrt(base::rowSums(rhy_params^2))
    phases_B <- base::atan2(rhy_params[2], rhy_params[1])
  } else {
    amps_B <- 0
    phases_B <- 0
  }

  if (all(base::is.element(c("inphase", "outphase"),
                           colnames(coeffs)))) {
    rhy_params <- coeffs[, c("inphase", "outphase")]
    amps <- 2 * sqrt(base::rowSums(rhy_params^2))
    phases <- base::atan2(rhy_params[2], rhy_params[1])
    model_params <- base::data.frame(amps, phases, amps, phases)
  } else {
    model_params <- base::data.frame(amps_A, phases_A, amps_B, phases_B)
  }

  colnames(model_params) <- base::paste(rep(group_ID, each = 2),
                                        c("amp","phase"), sep = "_")

  return(model_params)
}
