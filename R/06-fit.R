#' Fit a PAC model using GenSA (if available)
#' @param data_train data.frame with birthage,censoredage,DthIndicator
#' @param patience integer early-stop counter
#' @param tol numeric tolerance
#' @param verbose logical
#' @export
setMethod("fit", "PACFit", function(model_fit, data_train, patience = 30, tol = 1e-6, verbose = FALSE) {
  if (!requireNamespace("GenSA", quietly = TRUE)) {
    stop("Package 'GenSA' is required for fitting but is not installed.", call. = FALSE)
  }
  PAC <- model_fit@model
  raw_fn <- function(par) {
    -sum(logllh(PAC, par = par,
                x = data_train$birthage,
                r = data_train$censoredage,
                i = data_train$DthIndicator))
  }
  lower <- PAC@bounds$lower
  upper <- PAC@bounds$upper
  initial <- PAC@params

  # early-stopping wrapper
  best_val <- Inf; counter <- 0L
  patience_fn <- function(par) {
    val <- raw_fn(par)
    if ((best_val - val) > tol) { best_val <<- val; counter <<- 0L } else { counter <<- counter + 1L }
    if (counter >= patience) return(val + 1e9) else return(val)
  }

  result <- GenSA::GenSA(
    par = initial, fn = patience_fn, lower = lower, upper = upper,
    control = list(max.call = 2000L, verbose = if (isTRUE(verbose)) 1L else 0L)
  )
  model_fit@model@params <- result$par
  model_fit@fitted <- TRUE
  model_fit@fit_results <- result
  model_fit
})

#' @export
setMethod("AIC", "PACFit", function(object, ...) {
  if (!isTRUE(object@fitted)) stop("Model has not been fitted yet.")
  k <- length(object@model@params)
  loglik <- -object@fit_results$value
  -2 * loglik + 2 * k
})
