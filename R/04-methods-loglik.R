#' @include 00-generics.R 01-classes.R 02-model-fns.R
NULL

#' Log-likelihood for PAC (vectorized over observations)
#' @param model PAC
#' @param par numeric parameter vector (defaults to model@params if omitted)
#' @param x birth age
#' @param r censored/event age
#' @param i event indicator (1=death, 0=censored)
#' @export
setMethod("logllh", "PAC", function(model, par, x, r, i) {
  if (missing(par) || is.null(par)) par <- model@params
  if (any(!is.finite(par)) || length(par) < 2L || abs(par[2]) <= 1e-6) {
    return(rep(-Inf, length(x)))
  }
  H   <- model@cumhaz_(par, x, r - x)
  muR <- model@mu_(par, r)
  mu0 <- model@mu_(par, 0)

  bad <- !is.finite(H) | H <= 0 | !is.finite(muR) | muR <= 0 | !is.finite(mu0) | mu0 <= 0
  ll  <- i * log(muR) - H
  ll[bad] <- -Inf
  ll
})

#' Dataset log-likelihood extractor for PAC
#' @export
setMethod("get_logllh", "PAC", function(model, data_train, par = NULL) {
  if (is.null(par)) par <- model@params
  logllh(model, par = par,
         x = data_train$birthage,
         r = data_train$censoredage,
         i = data_train$DthIndicator)
})
