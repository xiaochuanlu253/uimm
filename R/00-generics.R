# Keep all S4 generics here so they load before any setMethod()

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

# (package-wide evaluation helpers)
#' @export
setGeneric("evaluate", function(object, ...) standardGeneric("evaluate"))

# ---- The two missing ones causing your error ----
#' @export
setGeneric("getIBS", function(object, ...) standardGeneric("getIBS"))

#' @export
setGeneric("getlogScore", function(object, ...) standardGeneric("getlogScore"))
