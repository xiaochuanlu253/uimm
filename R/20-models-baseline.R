#' @exportClass MortalityModel
setClass(
  "MortalityModel",
  slots = list(
    data = "list",
    data_train = "data.frameOrNULL",
    data_test = "data.frameOrNULL",
    Surv_info = "list"
  ),
  prototype = list(Surv_info = list())
)

#' @export
MortalityModel <- function(data = NULL, data_train = NULL, data_test = NULL) {
  if (is.data.frame(data)) { data_train <- data; data_test <- NULL; data <- list(train = data_train, test = data_test) }
  if (is.null(data)) data <- list(train = data_train, test = data_test)
  methods::new("MortalityModel",
               data = data,
               data_train = data$train,
               data_test = data$test,
               Surv_info = list()
  )
}

#' @export
setMethod("show", "MortalityModel", function(object) {
  cat("Training data size:", if (is.null(object@data_train)) 0L else nrow(object@data_train), "\n")
  cat("Test data size:",     if (is.null(object@data_test))  0L else nrow(object@data_test),  "\n")
})

#' @exportClass BaselineModel
setClass(
  "BaselineModel",
  slots = list(
    IMM = "MortalityModel",
    model_fit = "PACFit"
  )
)

#' @export
BaselineModel <- function(data = NULL, data_train = NULL, data_test = NULL, AC = BeardPAC()) {
  IMM <- MortalityModel(data = data, data_train = data_train, data_test = data_test)
  set.seed(5346615)
  pac_fit_obj <- methods::new("PACFit", model = AC, data = IMM@data)
  pac_fit_obj <- fit(pac_fit_obj, IMM@data_train, patience = 500)
  methods::new("BaselineModel", IMM = IMM, model_fit = pac_fit_obj)
}

#' @export
setMethod("show", "BaselineModel", function(object) {
  cat("Training data size:", nrow(object@IMM@data_train), "\n")
  cat("Test data size:",     if (is.null(object@IMM@data_test)) 0L else nrow(object@IMM@data_test), "\n")
  cat("Model params:", paste(signif(object@model_fit@model@params,6), collapse=" "), "\n")
})

#' @export
setMethod("predict", "BaselineModel",
          function(object, newdata = NULL,
                   type = c("survival","hazard","loglik"),
                   age_range = 80:115) {
            type <- match.arg(type)
            if (is.null(newdata)) newdata <- object@IMM@data_test
            predict.PACFit(object@model_fit, newdata = newdata, type = type, age_range = age_range)
          })

#' @export
setMethod("evaluate", "BaselineModel",
          function(object, newdata = NULL,
                   metrics = c("AIC","loglik_sum","loglik_vec","IBS"),
                   ibs_method = c("predict","ipcw"),
                   age_range = 80:115,
                   normalize_IBS = TRUE,
                   ...) {
            if (is.null(newdata)) newdata <- object@IMM@data_test
            ibs_method <- match.arg(ibs_method)
            # reuse evaluate.PACFit
            evaluate.PACFit(object@model_fit, newdata = newdata, metrics = metrics,
                            ibs_method = ibs_method, age_range = age_range,
                            normalize_IBS = normalize_IBS, ...)
          })

#' @export
setMethod("getIBS", "BaselineModel", function(object) {
  stop("Use evaluate() with ibs_method='ipcw' or 'predict'.")
})
#' @export
setMethod("getlogScore", "BaselineModel", function(object) {
  sum(get_logllh(object@model_fit@model, object@IMM@data_test))
})
