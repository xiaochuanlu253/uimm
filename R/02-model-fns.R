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
