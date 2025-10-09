#' @title PAC and PACFit S4 Classes
#' @description S4 classes for parametric age-cohort (PAC) survival models and fits.
#' @slot name model name
#' @slot params numeric vector of parameters
#' @slot survival_ function(par, x) -> S(x)
#' @slot mu_ function(par, x) -> hazard at x
#' @slot cumhaz_ function(par, x, t) -> cumulative hazard from x to x+t
#' @slot bounds list(lower=, upper=)
#' @exportClass PAC
setClass(
  "PAC",
  slots = list(
    name = "character",
    params = "numeric",
    survival_ = "function",
    mu_ = "function",
    cumhaz_ = "function",
    bounds = "list"
  )
)

#' @slot model a PAC
#' @slot data list(train=, test=)
#' @slot data_train data.frame or NULL
#' @slot data_test  data.frame or NULL
#' @slot fit_results list
#' @slot fitted logical
#' @exportClass PACFit
setClassUnion("data.frameOrNULL", c("data.frame", "NULL"))
setClass(
  "PACFit",
  slots = list(
    model = "PAC",
    data = "list",
    data_train = "data.frameOrNULL",
    data_test  = "data.frameOrNULL",
    fit_results = "list",
    fitted = "logical"
  )
)

#' @export
setMethod("show", "PAC", function(object) {
  cat("PAC Model:", object@name, "\n")
})

#' @export
setMethod("show", "PACFit", function(object) {
  cat("PACFit Model:", object@model@name, "\n")
  if (isTRUE(object@fitted)) {
    if (!is.null(object@fit_results$value)) cat("NegLogLik:", object@fit_results$value, "\n")
    cat("Parameters:", paste(signif(object@model@params, 6), collapse = " "), "\n")
  }
})

# ---- generics ----

#' @export
setGeneric("mu", function(model, ...) standardGeneric("mu"))
#' @export
setGeneric("survival", function(model, ...) standardGeneric("survival"))
#' @export
setGeneric("cumhaz", function(model, ...) standardGeneric("cumhaz"))
#' @export
setGeneric("logllh", function(model, ...) standardGeneric("logllh"))
#' @export
setGeneric("get_logllh", function(model, ...) standardGeneric("get_logllh"))
#' @export
setGeneric("fit", function(model_fit, ...) standardGeneric("fit"))
#' @export
setGeneric("BS", function(model_fit, ...) standardGeneric("BS"))
