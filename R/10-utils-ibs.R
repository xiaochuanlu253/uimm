#' Empirical survival cache retriever
#' @keywords internal
.empirical_from_cache <- function(key) {
  if (!exists(key, envir = .pacsurv_cache, inherits = FALSE)) return(NULL)
  get(key, envir = .pacsurv_cache, inherits = FALSE)
}

#' Store in cache
#' @keywords internal
.empirical_to_cache <- function(key, val) assign(key, val, envir = .pacsurv_cache)

#' Compute empirical survival info for IBS weighting
#' @param data data.frame
#' @param key character cache key
#' @param x_min,x_max integer age range (inclusive)
#' @return list with survfit objects and step values
#' @export
empirical_survival <- function(data, key, x_min = 80, x_max = 115) {
  stopifnot(is.data.frame(data))
  stopifnot(is.character(key), length(key) == 1L)
  cached <- .empirical_from_cache(key)
  if (!is.null(cached)) return(cached)

  BA <- survival::Surv(time = data$birthage, event = rep(1, nrow(data)))
  sf_BA <- survival::survfit(BA ~ 1)
  BA_vals <- summary(sf_BA, times = x_min:x_max)$surv

  CA <- survival::Surv(time = data$censoredage, event = 1 - data$DthIndicator)
  sf_CA <- survival::survfit(CA ~ 1)
  CA_vals <- summary(sf_CA, times = x_min:x_max)$surv

  out <- list(survfit_BA = sf_BA, survfit_CA = sf_CA, sf_BA = BA_vals, sf_CA = CA_vals)
  .empirical_to_cache(key, out)
  out
}

#' Raw Brier Score contribution (single i, single age)
#' @keywords internal
BS_raw <- function(x, birthage, censoredage, DthIndicator, Sfunc, par, Surv_info, x_min = 80) {
  S <- Sfunc(par, x); S0 <- Sfunc(par, birthage); S <- S / S0
  if (x < birthage) return(0)
  if (censoredage <= x && DthIndicator == 1) {
    rawbs <- S^2
    CA_fit <- Surv_info$survfit_CA
    if (max(CA_fit$time) <= censoredage) {
      SCA <- 1 / CA_fit$n
    } else {
      SCA <- summary(CA_fit, times = censoredage)$surv
    }
    return(rawbs / SCA)
  }
  if (censoredage > x) {
    rawbs <- (1 - S)^2
    idx <- as.integer(x - x_min + 1L)
    sfBA <- Surv_info$sf_BA[idx]; sfCA <- Surv_info$sf_CA[idx]
    bs <- rawbs / (1 - sfBA) / sfCA
    if (!is.finite(bs) || is.nan(bs) || is.na(bs)) {
      if (max(Surv_info$survfit_BA$time) <= x) sfBA <- 0
      if (max(Surv_info$survfit_CA$time) <= x) sfCA <- 1 / Surv_info$survfit_CA$n
      bs <- rawbs / (1 - sfBA) / sfCA
    }
    return(bs)
  }
  0
}

#' @export
setMethod("BS", "PACFit", function(model_fit, data_valid, Surv_info, x_min = 80, x_max = 115) {
  if (!isTRUE(model_fit@fitted)) stop("Model has not been fitted yet.")
  ages <- x_min:x_max
  n <- nrow(data_valid); m <- length(ages)
  acc <- 0
  for (age in ages) {
    tmp <- mapply(function(ba, ca, di) {
      BS_raw(age, ba, ca, di, model_fit@model@survival_, model_fit@model@params, Surv_info, x_min)
    }, data_valid$birthage, data_valid$censoredage, data_valid$DthIndicator)
    acc <- acc + mean(tmp)
  }
  acc / m
})
