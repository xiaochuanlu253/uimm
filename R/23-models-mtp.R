methods::setOldClass("modelparty")
#' @exportClass MTPModel
setClass(
  "MTPModel",
  slots = list(
    IMM_filled = "MortalityModel",
    formula = "formula",
    mobfit = "modelparty",
    mob_coef = "matrix",
    mob_tnodes_test = "numeric"
  )
)

#' @export
MTPModel <- function(formula, data = NULL, data_train = NULL, data_test = NULL) {
  IMM <- MortalityModel(data = data, data_train = data_train, data_test = data_test)

  gpreg <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
    flexsurv::flexsurvreg(y ~ 0 + x, weights = weights, dist = "gompertz", ...)
  }
  mobfit <- partykit::mob(
    formula,
    data = .fill_na_cats(IMM@data_train),
    fit = gpreg,
    control = partykit::mob_control(minsize = nrow(IMM@data_train) / 50, maxdepth = 10)
  )
  mob_coef <- coef(mobfit)
  test_filled <- .fill_na_cats(IMM@data_test)
  mob_nodes <- if (is.null(test_filled)) numeric(0) else predict(mobfit, test_filled, type = "node")

  if (length(unique(mob_nodes)) <= 1L) {
    # fallback baseline
    GompertzFit <- methods::new("PACFit", model = GompertzPAC(), data = IMM@data)
    GompertzFit@model@params <- c(mob_coef[2], mob_coef[1])
    GompertzFit@fitted <- TRUE
    return(methods::new("BaselineModel", IMM = IMM, model_fit = GompertzFit))
  }

  methods::new("MTPModel",
               IMM_filled = IMM, formula = formula,
               mobfit = mobfit, mob_coef = mob_coef, mob_tnodes_test = mob_nodes)
}

.fill_na_cats <- function(df) {
  if (is.null(df)) return(NULL)
  df2 <- df
  if (ncol(df2) >= 5) {
    for (i in 5:ncol(df2)) {
      xi <- df2[[i]]
      if (is.factor(xi) || is.character(xi)) {
        xi[is.na(xi)] <- "Missing"; df2[[i]] <- factor(xi)
      }
    }
  }
  df2
}

#' @export
setMethod("predict", "MTPModel",
          function(object, newdata = NULL,
                   type = c("survival","hazard","loglik", "group", "par"),
                   age_range = 80:115) {
            type <- match.arg(type)
            ages <- sort(as.numeric(age_range))
            tnodes <- rownames(object@mob_coef)
            if (is.null(newdata)) {
              newdata <- object@IMM_filled@test_fillna
              newdata$group <- object@mob_tnodes_test
            } else {
              newdata$group <- predict(object@mobfit, newdata, type = "node")
            }
            mob_coef <- object@mob_coef[, 2:1]  # reorder to (a, b)

            if (type %in% c("survival","hazard")) {
              M <- matrix(NA_real_, nrow(newdata), length(ages))
              for (nodei in tnodes) {
                idx <- which(newdata$group == nodei)
                if (!length(idx)) next
                pac <- GompertzPAC()
                pac@params <- mob_coef[match(nodei, tnodes), ]
                if (type == "survival") {
                  S_age <- pac@survival_(pac@params, ages)
                  S_ba  <- pac@survival_(pac@params, newdata$birthage[idx])
                  B     <- matrix(rep(S_age, each = length(idx)), nrow = length(idx))
                  B     <- sweep(B, 1, S_ba, "/")
                  for (k in seq_along(idx)) {
                    younger <- which(ages < newdata$birthage[idx[k]])
                    if (length(younger)) B[k, younger] <- 1
                  }
                  M[idx, ] <- B
                } else {
                  mu_age <- pac@mu_(pac@params, ages)
                  M[idx, ] <- matrix(rep(mu_age, each = length(idx)), nrow = length(idx))
                }
              }
              colnames(M) <- as.character(ages)
              return(M)
            }

            # type == "loglik"
            ll <- rep(NA_real_, nrow(newdata))
            for (nodei in tnodes) {
              idx <- which(newdata$group == nodei)
              if (!length(idx)) next
              pac <- GompertzPAC()
              pac@params <- mob_coef[match(nodei, tnodes), ]
              ll[idx] <- logllh(pac,
                                par = pac@params,
                                x = newdata$birthage[idx],
                                r = newdata$censoredage[idx],
                                i = newdata$DthIndicator[idx])
            }

            if (type == "group") return( newdata$group)
            if (type == "par")   return( mob_coef[match(newdata$group, tnodes), , drop = FALSE])

            names(ll) <- rownames(newdata)
            ll
          })

#' @export
setMethod("evaluate", "MTPModel",
          function(object, newdata = NULL,
                   metrics = c("AIC","loglik_sum","loglik_vec","IBS"),
                   age_range = 80:115, normalize_IBS = TRUE, ...) {
            if (is.null(newdata)) {
              newdata <- object@IMM_filled@test_fillna
              newdata$group <- object@mob_tnodes_test
            }
            ages <- sort(as.numeric(age_range))
            out <- list()

            # AIC not defined for the piecewise-node model; return NA
            if ("AIC" %in% metrics) out$AIC <- NA_real_

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
