compute_circ_params <- function(y, t, period) {

  inphase <- cos(2 * pi * t / period)
  outphase <- sin(2 * pi * t / period)

  X <- stats::model.matrix(~inphase + outphase)

  fit <- stats::lm.fit(X, t(y))

  amps <- 2 * sqrt(base::colSums(fit$coefficients[-1, ]^2))

  phases <- (atan2(fit$coefficients[3, ], fit$coefficients[2, ]) %% (2*pi))

  return(base::cbind(amps = amps, phases = phases))

}

compute_model_params <- function(y, group_id, d=NULL, type="fit") {

  if (type == "fit") {
    fit <- limma::lmFit(y, d)
    coeffs <- fit$coefficients
  }
  else if (type == "coef") {
    coeffs <- y
  }

  if (any(base::grepl(paste0(group_id[1], "_"),
                      colnames(coeffs)))) {
    rhy_params <- coeffs[, base::paste(group_id[1],
                                       c("inphase", "outphase"),
                                       sep = "_")]
    amps_A <- 2 * sqrt(base::rowSums(rhy_params^2))
    phases_A <- base::atan2(rhy_params[, 2], rhy_params[, 1]) %% (2*pi)
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
    phases_B <- base::atan2(rhy_params[, 2], rhy_params[, 1])  %% (2*pi)
  } else {
    amps_B <- 0
    phases_B <- 0
  }

  if (all(base::is.element(c("inphase", "outphase"),
                           colnames(coeffs)))) {
    rhy_params <- coeffs[, c("inphase", "outphase")]
    amps <- 2 * sqrt(base::rowSums(rhy_params^2))
    phases <- base::atan2(rhy_params[, 2], rhy_params[, 1])  %% (2*pi)
    model_params <- base::cbind(amps, phases, amps, phases)
  } else {
    model_params <- base::cbind(amps_A, phases_A, amps_B, phases_B)
  }

  colnames(model_params) <- base::paste(rep(group_id, each = 2),
                                        c("amp", "phase"), sep = "_")

  return(model_params)
}

categorize <- function(a, b, dr) {
  if (a && !b && dr) {
    category <- "loss"
  } else if (!a && b && dr) {
    category <- "gain"
  } else if (a && b && dr) {
    category <- "change"
  } else {
    category <- "same"
  }
  return(category)
}
