# Gompertz
survival_gompertz <- function(par, x) {
  a <- par[1]; b <- par[2]
  exp(-exp(a) * (exp(b * x) - 1) / b)
}
cumulative_hazard_gompertz <- function(par, x, t) {
  a <- par[1]; b <- par[2]
  (exp(b * t) - 1) / b * exp(a + b * x)
}
forceofmortality_gompertz <- function(par, x) {
  a <- par[1]; b <- par[2]
  exp(a + b * x)
}

# Makeham
survival_makeham <- function(par, x) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  exp(-exp(a) * (exp(b * x) - 1) / b - x * lambda)
}
cumulative_hazard_makeham <- function(par, x, t) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  (exp(b * t) - 1) / b * exp(a + b * x) + lambda * t
}
forceofmortality_makeham <- function(par, x) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  exp(a + b * x) + lambda
}

# Perks
survival_perks <- function(par, x) {
  a <- par[1]; b <- par[2]; gam <- par[3]; lambda <- par[4]
  ((1 + exp(a + gam)) / (1 + exp(a + gam + b * x)))^(exp(-gam) / b) / exp(lambda * x)
}
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
forceofmortality_perks <- function(par, x) {
  a <- par[1]; b <- par[2]; gam <- par[3]; lambda <- par[4]
  if (b > 1) 1 / (exp(-(a + b * x)) + exp(gam)) + lambda else exp(a + b * x) / (1 + exp(a + gam + b * x)) + lambda
}

# Beard
survival_beard <- function(par, x) {
  a <- par[1]; b <- par[2]; gam <- par[3]
  ((1 + exp(a + gam)) / (1 + exp(a + gam + b * x)))^(exp(-gam) / b)
}
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
forceofmortality_beard <- function(par, x) {
  a <- par[1]; b <- par[2]; gam <- par[3]
  if (b > 1) 1 / (exp(-(a + b * x)) + exp(gam)) else exp(a + b * x) / (1 + exp(a + gam + b * x))
}

# Kannisto
survival_kannisto <- function(par, x) {
  a <- par[1]; b <- par[2]
  ((1 + exp(a)) / (1 + exp(a + b * x)))^(1 / b)
}
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
forceofmortality_kannisto <- function(par, x) {
  a <- par[1]; b <- par[2]
  if (b > 1) 1 / (exp(-(a + b * x)) + 1) else exp(a + b * x) / (1 + exp(a + b * x))
}

# Thatcher
survival_thatcher <- function(par, x) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  ((1 + exp(a)) / (1 + exp(a + b * x)))^(1 / b) / exp(lambda * x)
}
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
