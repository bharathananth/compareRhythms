compute_circ_params <- function(y, t, period) {

  inphase <- cos(2 * pi * t / period)
  outphase <- sin(2 * pi * t / period)

  X <- stats::model.matrix(~inphase + outphase)

  fit <- stats::lm.fit(X, t(y))

  amps <- 2 * sqrt(base::colSums(fit$coefficients[-1, ]^2))

  phases <- atan2(fit$coefficients[3, ], fit$coefficients[2, ])

  return(base::cbind(amps = amps, phases = phases))

}
