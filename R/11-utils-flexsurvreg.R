# Log-likelihood extractor for flexsurvreg
logLik.flexsurvreg <- function(object, ...) {
  structure(object$loglik, df = object$npars, class = "logLik")
}
# Empirical estimating function for flexsurvreg (Gompertz only)
estfun.flexsurvreg <- function(object, ...) {
  if (object$dlist$name == "gompertz") {
    Y <- object$data$Y[, 1:3]
    Li <- Y[, 1]
    Ci <- Y[, 2]
    di <- Y[, 3]
    a <- object$coefficients[1]
    log_b <- object$coefficients[2]
    b <- exp(log_b)
    est_shape <- di * Ci - b / (a^2) * (exp(a * Ci) * (a * Ci - 1) - exp(a * Li) * (a * Li - 1))
    est_rate <- di / b - (exp(a * Ci) - exp(a * Li)) / a
    return(cbind(est_shape, est_rate))
  }
}
