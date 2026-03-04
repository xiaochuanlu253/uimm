#' Gompertz mortality law
#'
#' The Gompertz model assumes a hazard that grows exponentially with age.
#'
#' @details Parameterization: \eqn{\mu(x) = \exp(a + b x)}.
#' The survival is \eqn{S(x) = \exp\{-\exp(a) (e^{b x} - 1)/b\}}.
#'
#' @param par numeric length-2 vector \code{c(a, b)}.
#' @param x numeric vector of ages or durations since origin.
#' @return For \code{survival_gompertz}: numeric vector of survival probabilities.
#' @family mortality-laws
#' @examples
#' survival_gompertz(c(-10, 0.09), x = 80:100)
#' @export
survival_gompertz <- function(par, x) {
  a <- par[1]; b <- par[2]
  exp(-exp(a) * (exp(b * x) - 1) / b)
}

#' @rdname survival_gompertz
#' @param t numeric vector of time intervals; used by \code{cumulative_hazard_gompertz}.
#' @return For \code{cumulative_hazard_gompertz}: numeric vector of cumulative hazards over \code{t}.
#' @export
cumulative_hazard_gompertz <- function(par, x, t) {
  a <- par[1]; b <- par[2]
  (exp(b * t) - 1) / b * exp(a + b * x)
}

#' @rdname survival_gompertz
#' @return For \code{forceofmortality_gompertz}: numeric vector of hazards \eqn{\mu(x)}.
#' @export
forceofmortality_gompertz <- function(par, x) {
  a <- par[1]; b <- par[2]
  exp(a + b * x)
}


#' Makeham mortality law
#'
#' Gompertz with an age-independent background hazard \eqn{\lambda}.
#'
#' @details \eqn{\mu(x) = \exp(a + b x) + \lambda}.
#' @param par numeric length-3 vector \code{c(a, b, lambda)}.
#' @inheritParams survival_gompertz
#' @family mortality-laws
#' @examples
#' survival_makeham(c(-10, 0.09, 0.0005), x = 80:100)
#' @export
survival_makeham <- function(par, x) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  exp(-exp(a) * (exp(b * x) - 1) / b - x * lambda)
}

#' @rdname survival_makeham
#' @export
cumulative_hazard_makeham <- function(par, x, t) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  (exp(b * t) - 1) / b * exp(a + b * x) + lambda * t
}

#' @rdname survival_makeham
#' @export
forceofmortality_makeham <- function(par, x) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  exp(a + b * x) + lambda
}


#' Perks mortality law
#'
#' Logistic hazard with late-age deceleration (Perks). Optional background \eqn{\lambda}.
#'
#' @details \eqn{\mu(x) = \frac{\exp(a + b x)}{1 + \exp(a + \gamma + b x)} + \lambda}.
#' @param par numeric length-4 vector \code{c(a, b, gam, lambda)} where \code{gam = \gamma}.
#' @inheritParams survival_gompertz
#' @family mortality-laws
#' @examples
#' survival_perks(c(-10, 0.09, -3, 0), x = 80:100)
#' @export
survival_perks <- function(par, x) {
  a <- par[1]; b <- par[2]; gam <- par[3]; lambda <- par[4]
  ((1 + exp(a + gam)) / (1 + exp(a + gam + b * x)))^(exp(-gam) / b) / exp(lambda * x)
}

#' @rdname survival_perks
#' @export
cumulative_hazard_perks <- function(par, x, t) {
  a <- par[1]; b <- par[2]; gam <- par[3]; lambda <- par[4]
  tolog <- if (b > 2) {
    exp(b * t) + exp(-(a + gam + b * (x - t)))
  } else if (b > 1) {
    exp(b * t) + (exp(-b * t) - 1) / (exp(-b * t) + exp(a + gam + b * (x - t)))
  } else {
    (1 + exp(a + gam + b * (x + t))) / (1 + exp(a + gam + b * x))
  }
  exp(-gam) / b * log(tolog) + lambda * t
}

#' @rdname survival_perks
#' @export
forceofmortality_perks <- function(par, x) {
  a <- par[1]; b <- par[2]; gam <- par[3]; lambda <- par[4]
  if (b > 1) 1 / (exp(-(a + b * x)) + exp(gam)) + lambda else exp(a + b * x) / (1 + exp(a + gam + b * x)) + lambda
}


#' Beard mortality law
#'
#' Logistic hazard with late-age leveling (Beard).
#'
#' @details \eqn{\mu(x) = \frac{\exp(a + b x)}{1 + \exp(a + \gamma + b x)}}.
#' @param par numeric length-3 vector \code{c(a, b, gam)} where \code{gam = \gamma}.
#' @inheritParams survival_gompertz
#' @family mortality-laws
#' @examples
#' survival_beard(c(-10, 0.09, -3), x = 80:100)
#' @export
survival_beard <- function(par, x) {
  a <- par[1]; b <- par[2]; gam <- par[3]
  ((1 + exp(a + gam)) / (1 + exp(a + gam + b * x)))^(exp(-gam) / b)
}

#' @rdname survival_beard
#' @export
cumulative_hazard_beard <- function(par, x, t) {
  a <- par[1]; b <- par[2]; gam <- par[3]
  tolog <- if (b > 2) {
    exp(b * t) + exp(-(a + gam + b * (x - t)))
  } else if (b > 1) {
    exp(b * t) + (exp(-b * t) - 1) / (exp(-b * t) + exp(a + gam + b * (x - t)))
  } else {
    (1 + exp(a + gam + b * (x + t))) / (1 + exp(a + gam + b * x))
  }
  exp(-gam) / b * log(tolog)
}

#' @rdname survival_beard
#' @export
forceofmortality_beard <- function(par, x) {
  a <- par[1]; b <- par[2]; gam <- par[3]
  if (b > 1) 1 / (exp(-(a + b * x)) + exp(gam)) else exp(a + b * x) / (1 + exp(a + gam + b * x))
}


#' Kannisto mortality law
#'
#' Two-parameter logistic mortality with late-age deceleration.
#'
#' @details \eqn{\mu(x) = \frac{\exp(a + b x)}{1 + \exp(a + b x)}}.
#' @param par numeric length-2 vector \code{c(a, b)}.
#' @inheritParams survival_gompertz
#' @family mortality-laws
#' @examples
#' survival_kannisto(c(-10, 0.09), x = 80:100)
#' @export
survival_kannisto <- function(par, x) {
  a <- par[1]; b <- par[2]
  ((1 + exp(a)) / (1 + exp(a + b * x)))^(1 / b)
}

#' @rdname survival_kannisto
#' @export
cumulative_hazard_kannisto <- function(par, x, t) {
  a <- par[1]; b <- par[2]
  tolog <- if (b > 2) {
    exp(b * t) + exp(-(a + b * (x - t)))
  } else if (b > 1) {
    exp(b * t) + (exp(-b * t) - 1) / (exp(-b * t) + exp(a + b * (x - t)))
  } else {
    (1 + exp(a + b * (x + t))) / (1 + exp(a + b * x))
  }
  (1 / b) * log(tolog)
}

#' @rdname survival_kannisto
#' @export
forceofmortality_kannisto <- function(par, x) {
  a <- par[1]; b <- par[2]
  if (b > 1) 1 / (exp(-(a + b * x)) + 1) else exp(a + b * x) / (1 + exp(a + b * x))
}


#' Thatcher mortality law
#'
#' Kannisto with an age-independent background hazard \eqn{\lambda}.
#'
#' @details \eqn{\mu(x) = \frac{\exp(a + b x)}{1 + \exp(a + b x)} + \lambda}.
#' @param par numeric length-3 vector \code{c(a, b, lambda)}.
#' @inheritParams survival_gompertz
#' @family mortality-laws
#' @examples
#' survival_thatcher(c(-10, 0.09, 0.0005), x = 80:100)
#' @export
survival_thatcher <- function(par, x) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  ((1 + exp(a)) / (1 + exp(a + b * x)))^(1 / b) / exp(lambda * x)
}

#' @rdname survival_thatcher
#' @export
cumulative_hazard_thatcher <- function(par, x, t) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  tolog <- if (b > 2) {
    exp(b * t) + exp(-(a + b * (x - t)))
  } else if (b > 1) {
    exp(b * t) + (exp(-b * t) - 1) / (exp(-b * t) + exp(a + b * (x - t)))
  } else {
    (1 + exp(a + b * (x + t))) / (1 + exp(a + b * x))
  }
  (1 / b) * log(tolog) + lambda * t
}

#' @rdname survival_thatcher
#' @export
forceofmortality_thatcher <- function(par, x) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  if (b > 1) 1 / (exp(-(a + b * x)) + 1) + lambda else exp(a + b * x) / (1 + exp(a + b * x)) + lambda
}


# method bindings
#' @export
setMethod("mu", "PAC", function(model, x) model@mu_(model@params, x))
#' @export
setMethod("survival", "PAC", function(model, x) model@survival_(model@params, x))
#' @export
setMethod("cumhaz", "PAC", function(model, x, t) model@cumhaz_(model@params, x, t))


# stats form distribution functions

#' The Gompertz distribution
#'
#' Density, distribution function, quantile function and random generation for the Gompertz distribution.
#' @details The Gompertz distribution is defined by the hazard function \eqn{\mu(x) = \exp(a + b x)}.
#' The survival function is \eqn{S(x) = \exp\{-\exp(a) (e^{b x} - 1)/b\}}.
#' @param x,q numeric vector of quantiles.
#' @param p numeric vector of probabilities.
#' @param n number of observations.
#' @param par numeric length-2 vector of parameters \code{c(a, b)}.
#' @param lower.tail logical; if \code{FALSE}, probabilities are \eqn{P(X > x)}.
#' @return
#' * `dgompertz()`: density.
#' * `pgompertz()`: distribution function.
#' * `qgompertz()`: quantiles.
#' * `rgompertz()`: random values.
#'
#' @seealso [`survival_gompertz()`] given the same parameterization.
#' @seealso [flexsurv::rgompertz()] for an implementation in \pkg{flexsurv}.
#'
#' @examples
#' dgompertz(80:100, c(-10, 0.09))
#' pgompertz(80:100, c(-10, 0.09))
#' qgompertz(seq(0, 1, by = 0.1), c(-10, 0.09))
#' rgompertz(10, c(-10, 0.09))
#'
#' @name Gompertz
#' @aliases dgompertz pgompertz qgompertz rgompertz
#'
#' @rdname Gompertz
#' @export
dgompertz <- function(x, par) {
  forceofmortality_gompertz(par, x) * survival_gompertz(par, x)
}
#' @rdname Gompertz
#' @export
pgompertz <- function(x, par) {
  1 - survival_gompertz(par, x)
}
#' @rdname Gompertz
#' @export
qgompertz <- function(p, par, lower.tail = TRUE) {
  if (!lower.tail) p <- 1 - p
  uniroot(function(x) survival_gompertz(par, x) - (1 - p), lower = 0, upper = 150)$root
}
#' @rdname Gompertz
#' @export
rgompertz <- function(n, par) {
  u <- runif(n)
  qgompertz(u, par)
}

#' The Makeham distribution
#'
#' Density, distribution function, quantile function and random generation for the Makeham distribution.
#' @details The Makeham distribution is defined by the hazard function \eqn{\mu(x) = \exp(a + b x) + \lambda}.
#' The survival function is \eqn{S(x) = \exp\{-\exp(a) (e^{b x} - 1)/b - x \lambda\}}.
#' @param x,q numeric vector of quantiles.
#' @param p numeric vector of probabilities.
#' @param n number of observations.
#' @param par numeric length-3 vector of parameters \code{c(a, b, lambda)}.
#' @param lower.tail logical; if \code{FALSE}, probabilities are \eqn{P(X > x)}.
#' @return
#' * `dmakeham()`: density.
#' * `pmakeham()`: distribution function.
#' * `qmakeham()`: quantiles.
#' * `rmakeham()`: random values.
#' @examples
#' dmakeham(80:100, c(-10, 0.09, 0.0005))
#' pmakeham(80:100, c(-10, 0.09, 0.0005))
#' qmakeham(seq(0, 1, by = 0.1), c(-10, 0.09, 0.0005))
#' rmakeham(10, c(-10, 0.09, 0.0005))
#' @name Makeham
#' @aliases dmakeham pmakeham qmakeham rmakeham
#' @rdname Makeham
#' @export
dmakeham <- function(x, par) {
  forceofmortality_makeham(par, x) * survival_makeham(par, x)
}
#' @rdname Makeham
#' @export
pmakeham <- function(x, par) {
  1 - survival_makeham(par, x)
}
#' @rdname Makeham
#' @export
qmakeham <- function(p, par, lower.tail = TRUE) {
  if (!lower.tail) p <- 1 - p
  uniroot(function(x) survival_makeham(par, x) - (1 - p), lower = 0, upper = 150)$root
}
#' @rdname Makeham
#' @export
rmakeham <- function(n, par) {
  u <- runif(n)
  qmakeham(u, par)
}

#' The Perks distribution
#'
#' Density, distribution function, quantile function and random generation for the Perks distribution.
#' @details The Perks distribution is defined by the hazard function \eqn{\mu(x) = \frac{\exp(a + b x)}{1 + \exp(a + \gamma + b x)} + \lambda}.
#' The survival function is \eqn{S(x) = ((1 + \exp(a + \gamma)) / (1 + \exp(a + \gamma + b x)))^{\exp(-\gamma) / b} / \exp(\lambda x)}.
#' @param x,q numeric vector of quantiles.
#' @param p numeric vector of probabilities.
#' @param n number of observations.
#' @param par numeric length-4 vector of parameters \code{c(a, b, gam, lambda)}.
#' @param lower.tail logical; if \code{FALSE}, probabilities are \eqn{P(X > x)}.
#' @return
#' * `dperks()`: density.
#' * `pperks()`: distribution function.
#' * `qperks()`: quantiles.
#' * `rperks()`: random values.
#' @examples
#' dperks(80:100, c(-10, 0.09, -3, 0))
#' pperks(80:100, c(-10, 0.09, -3, 0))
#' qperks(seq(0, 1, by = 0.1), c(-10, 0.09, -3, 0))
#' rperks(10, c(-10, 0.09, -3, 0))
#' @name Perks
#' @aliases dperks pperks qperks rperks
#' @rdname Perks
#' @export
dperks <- function(x, par) {
  forceofmortality_perks(par, x) * survival_perks(par, x)
}
#' @rdname Perks
#' @export
pperks <- function(x, par) {
  1 - survival_perks(par, x)
}
#' @rdname Perks
#' @export
qperks <- function(p, par, lower.tail = TRUE) {
  if (!lower.tail) p <- 1 - p
  uniroot(function(x) survival_perks(par, x) - (1 - p), lower = 0, upper = 150)$root
}
#' @rdname Perks
#' @export
rperks <- function(n, par) {
  u <- runif(n)
  qperks(u, par)
}

#' The Beard distribution
#' Density, distribution function, quantile function and random generation for the Beard distribution.
#' @details The Beard distribution is defined by the hazard function \eqn{\mu(x) = \frac{\exp(a + b x)}{1 + \exp(a + \gamma + b x)}}.
#' The survival function is \eqn{S(x) = ((1 + \exp(a + \gamma)) / (1 + \exp(a + \gamma + b x)))^{\exp(-\gamma) / b}}.
#' @param x,q numeric vector of quantiles.
#' @param p numeric vector of probabilities.
#' @param n number of observations.
#' @param par numeric length-3 vector of parameters \code{c(a, b, gam)}.
#' @param lower.tail logical; if \code{FALSE}, probabilities are \eqn{P(X > x)}.
#' @return
#' * `dbeard()`: density.
#' * `pbeard()`: distribution function.
#' * `qbeard()`: quantiles.
#' * `rbeard()`: random values.
#' @examples
#' dbeard(80:100, c(-10, 0.09, -3))
#' pbeard(80:100, c(-10, 0.09, -3))
#' qbeard(seq(0, 1, by = 0.1), c(-10, 0.09, -3))
#' rbeard(10, c(-10, 0.09, -3))
#' @name Beard
#' @aliases dbeard pbeard qbeard rbeard
#' @rdname Beard
#' @export
dbeard <- function(x, par) {
  forceofmortality_beard(par, x) * survival_beard(par, x)
}
#' @rdname Beard
#' @export
pbeard <- function(x, par) {
  1 - survival_beard(par, x)
}
#' @rdname Beard
#' @export
qbeard <- function(p, par, lower.tail = TRUE) {
  if (!lower.tail) p <- 1 - p
  uniroot(function(x) survival_beard(par, x) - (1 - p), lower = 0, upper = 150)$root
}
#' @rdname Beard
#' @export
rbeard <- function(n, par) {
  u <- runif(n)
  qbeard(u, par)
}

#' The Kannisto distribution
#'
#' Density, distribution function, quantile function and random generation for the Kannisto distribution.
#' @details The Kannisto distribution is defined by the hazard function \eqn{\mu(x) = \frac{\exp(a + b x)}{1 + \exp(a + b x)}}.
#' The survival function is \eqn{S(x) = ((1 + \exp(a)) / (1 + \exp(a + b x)))^{1 / b}}.
#' @param x,q numeric vector of quantiles.
#' @param p numeric vector of probabilities.
#' @param n number of observations.
#' @param par numeric length-2 vector of parameters \code{c(a, b)}.
#' @param lower.tail logical; if \code{FALSE}, probabilities are \eqn{P(X > x)}.
#' @return
#' * `dkannisto()`: density.
#' * `pkannisto()`: distribution function.
#' * `qkannisto()`: quantiles.
#' * `rkannisto()`: random values.
#' @examples
#' dkannisto(80:100, c(-10, 0.09))
#' pkannisto(80:100, c(-10, 0.09))
#' qkannisto(seq(0, 1, by = 0.1), c(-10, 0.09))
#' rkannisto(10, c(-10, 0.09))
#' @name Kannisto
#' @aliases dkannisto pkannisto qkannisto rkannisto
#' @rdname Kannisto
#' @export
dkannisto <- function(x, par) {
  forceofmortality_kannisto(par, x) * survival_kannisto(par, x)
}
#' @rdname Kannisto
#' @export
pkannisto <- function(x, par) {
  1 - survival_kannisto(par, x)
}
#' @rdname Kannisto
#' @export
qkannisto <- function(p, par, lower.tail = TRUE) {
  if (!lower.tail) p <- 1 - p
  uniroot(function(x) survival_kannisto(par, x) - (1 - p), lower = 0, upper = 150)$root
}
#' @rdname Kannisto
#' @export
rkannisto <- function(n, par) {
  u <- runif(n)
  qkannisto(u, par)
}

#' The Thatcher distribution
#'
#' Density, distribution function, quantile function and random generation for the Thatcher distribution.
#' @details The Thatcher distribution is defined by the hazard function \eqn{\mu(x) = \frac{\exp(a + b x)}{1 + \exp(a + b x)} + \lambda}.
#' The survival function is \eqn{S(x) = ((1 + \exp(a)) / (1 + \exp(a + b x)))^{1 / b} / \exp(\lambda x)}.
#' @param x,q numeric vector of quantiles.
#' @param p numeric vector of probabilities.
#' @param n number of observations.
#' @param par numeric length-3 vector of parameters \code{c(a, b, lambda)}.
#' @param lower.tail logical; if \code{FALSE}, probabilities are \eqn{P(X > x)}.
#' @return
#' * `dthatcher()`: density.
#' * `pthatcher()`: distribution function.
#' * `qthatcher()`: quantiles.
#' * `rthatcher()`: random values.
#' @examples
#' dthatcher(80:100, c(-10, 0.09, 0.0005))
#' pthatcher(80:100, c(-10, 0.09, 0.0005))
#' qthatcher(seq(0, 1, by = 0.1), c(-10, 0.09, 0.0005))
#' rthatcher(10, c(-10, 0.09, 0.0005))
#' @name Thatcher
#' @aliases dthatcher pthatcher qthatcher rthatcher
#' @rdname Thatcher
#' @export
dthatcher <- function(x, par) {
  forceofmortality_thatcher(par, x) * survival_thatcher(par, x)
}
#' @rdname Thatcher
#' @export
pthatcher <- function(x, par) {
  1 - survival_thatcher(par, x)
}
#' @rdname Thatcher
#' @export
qthatcher <- function(p, par, lower.tail = TRUE) {
  if (!lower.tail) p <- 1 - p
  uniroot(function(x) survival_thatcher(par, x) - (1 - p), lower = 0, upper = 150)$root
}
#' @rdname Thatcher
#' @export
rthatcher <- function(n, par) {
  u <- runif(n)
  qthatcher(u, par)
}
