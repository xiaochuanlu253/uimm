methods::setOldClass("flexsurvreg")
#' @exportClass GPHModel
setClass(
  "GPHModel",
  slots = list(
    IMM = "MortalityModel",
    formula = "formula",
    gpfit = "flexsurvreg",
    gpfit_pred = "numeric",
    gpfit_shape = "numeric"
  )
)

#' @export
GPHModel <- function(formula, data = NULL, data_train = NULL, data_test = NULL) {
  IMM <- MortalityModel(data = data, data_train = data_train, data_test = data_test)
  set.seed(5346615)
  gpfit <- flexsurv::flexsurvreg(formula, data = .fill_na_cats(IMM@data_train), dist = "gompertz")
  lp <- predict(gpfit, newdata = .fill_na_cats(IMM@data_test), type = "lp")
  gpfit_pred <- log(lp$.pred_link)
  gpfit_shape <- gpfit$coefficients[1]
  methods::new("GPHModel", IMM = IMM, formula = formula, gpfit = gpfit,
               gpfit_pred = gpfit_pred, gpfit_shape = gpfit_shape)
}

# S3 registration helpers
#' @export
logLik.flexsurvreg <- function(object, ...) structure(object$loglik, df = object$npars, class = "logLik")
#' @export
estfun.flexsurvreg <- function(object, ...) {
  if (object$dlist$name != "gompertz") return(NULL)
  Y <- object$data$Y[, 1:3]
  Li <- Y[, 1]; Ci <- Y[, 2]; di <- Y[, 3]
  a <- object$coefficients[1]; log_b <- object$coefficients[2]; b <- exp(log_b)
  est_shape <- di * Ci - b / (a^2) * (exp(a * Ci) * (a * Ci - 1) - exp(a * Li) * (a * Li - 1))
  est_rate  <- di / b - (exp(a * Ci) - exp(a * Li)) / a
  cbind(est_shape, est_rate)
}

#' @export
setMethod("predict", "GPHModel",
          function(object, newdata = NULL,
                   type = c("survival","hazard","loglik", "par"),
                   age_range = 80:115) {
            type <- match.arg(type)
            ages <- sort(as.numeric(age_range))
            if (is.null(newdata)) newdata <- object@IMM_filled@test_fillna

            # shape (b) common; alpha_i varies via gpfit_pred (stored on object)
            b <- as.numeric(object@gpfit_shape)
            # If predicting on newdata different from test_fillna, recompute lp quickly:
            if (!identical(newdata, object@IMM_filled@test_fillna)) {
              # best-effort: use flexsurv predict with type = "lp"
              lp <- predict(object@gpfit, newdata = newdata, type = "lp")$.pred_link
              alpha_vec <- log(lp)
            } else {
              alpha_vec <- object@gpfit_pred
            }

            if (type %in% c("survival","hazard")) {
              M <- matrix(NA_real_, nrow(newdata), length(ages))
              for (i in seq_len(nrow(newdata))) {
                pac <- GompertzPAC()
                pac@params <- c(alpha_vec[i], b)
                if (type == "survival") {
                  S_age <- pac@survival_(pac@params, ages)
                  S_ba  <- pac@survival_(pac@params, newdata$birthage[i])
                  v     <- S_age / S_ba
                  v[ages < newdata$birthage[i]] <- 1
                  M[i, ] <- v
                } else {
                  M[i, ] <- pac@mu_(pac@params, ages)
                }
              }
              colnames(M) <- as.character(ages)
              return(M)
            }

            # type == "loglik" (per-individual)
            ll <- rep(NA_real_, nrow(newdata))
            for (i in seq_len(nrow(newdata))) {
              pac <- GompertzPAC()
              pac@params <- c(alpha_vec[i], b)
              ll[i] <- logllh(pac,
                              par = pac@params,
                              x = newdata$birthage[i],
                              r = newdata$censoredage[i],
                              i = newdata$DthIndicator[i])
            }

            if (type == "par") return( cbind(alpha = alpha_vec, shape = rep(b, nrow(newdata))) )

            names(ll) <- rownames(newdata)
            ll
          })

#' @export
setMethod("evaluate", "GPHModel",
          function(object, newdata = NULL,
                   metrics = c("AIC","loglik_sum","loglik_vec","IBS"),
                   age_range = 80:115, normalize_IBS = TRUE, ...) {
            if (is.null(newdata)) newdata <- object@IMM_filled@test_fillna
            ages <- sort(as.numeric(age_range))
            out <- list()

            # AIC available via flexsurvreg
            if ("AIC" %in% metrics) out$AIC <- tryCatch(stats::AIC(object@gpfit), error = function(e) NA_real_)

            if ("loglik_vec" %in% metrics || "loglik_sum" %in% metrics) {
              ll_vec <- predict(object, newdata, type = "loglik")
              if ("loglik_vec" %in% metrics) out$loglik_vec <- ll_vec
              if ("loglik_sum" %in% metrics) out$loglik_sum <- sum(ll_vec, na.rm = TRUE)
            }


            if ("IBS" %in% metrics) {
              S_pred <- predict(object, newdata, type = "survival", age_range = ages)
              ibs <- .ibs_from_Spred(newdata, S_pred, ages, normalize = normalize_IBS)
              out$IBS <- ibs$IBS
              out$IBS_by_age <- ibs$IBS_by_age
            }

            out
          })
