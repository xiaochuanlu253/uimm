#' @exportClass STPModel
setClass(
  "STPModel",
  slots = list(
    STModel = "STModel",
    AC = "PAC",
    par_train = "matrix"
  )
)

#' @export
STPModel <- function(data = NULL, data_train = NULL, data_test = NULL, AC = BeardPAC(), ST = NULL) {
  if (is.null(ST)) ST <- STModel(formula = survival::Surv(birthage, censoredage, DthIndicator) ~ 1,
                                 data = data, data_train = data_train, data_test = data_test)
  tnodes <- ST@tnodes
  par_train <- NULL
  for (nodei in tnodes) {
    set.seed(5346615)
    pac_fit_obj <- methods::new("PACFit", model = AC, data = ST@IMM@data)
    idx <- ST@IMM@data_train$node == nodei
    pac_fit_obj <- fit(pac_fit_obj, ST@IMM@data_train[idx, , drop = FALSE], patience = 500)
    par_train <- rbind(par_train, pac_fit_obj@model@params)
  }
  rownames(par_train) <- tnodes
  methods::new("STPModel", STModel = ST, AC = AC, par_train = par_train)
}

# show/getIBS/getlogScore similar to your originals but with base subsetting

#' @export
setMethod("predict", "STPModel",
          function(object, newdata = NULL,
                   type = c("survival","hazard","loglik", "group", "par"),
                   age_range = 80:115) {
            type <- match.arg(type)
            ages <- sort(as.numeric(age_range))
            if (is.null(newdata)) {
              newdata <- object@STModel@IMM@data_test
              if (is.null(newdata$group)) {
                newdata$group <- predict(object@STModel@model %>% partykit::as.party(),
                                         newdata = newdata, type = "node")
              }
            } else {
              newdata$group <- predict(object@STModel@model %>% partykit::as.party(),
                                       newdata = newdata, type = "node")
            }

            tnodes <- object@STModel@tnodes
            # Precompute per-node functions at 'ages'
            if (type %in% c("survival","hazard")) {
              S_or_mu <- matrix(NA_real_, nrow(newdata), length(ages))
              for (nodei in tnodes) {
                idx <- which(newdata$group == nodei)
                if (!length(idx)) next
                pac <- object@AC
                pac@params <- object@par_train[match(nodei, tnodes), ]
                if (type == "survival") {
                  S_age <- pac@survival_(pac@params, ages)             # |ages|
                  S_ba  <- pac@survival_(pac@params, newdata$birthage[idx])
                  M     <- matrix(rep(S_age, each = length(idx)), nrow = length(idx))
                  M     <- sweep(M, 1, S_ba, "/")
                  # ages < birthage -> 1
                  for (k in seq_along(idx)) {
                    younger <- which(ages < newdata$birthage[idx[k]])
                    if (length(younger)) M[k, younger] <- 1
                  }
                  S_or_mu[idx, ] <- M
                } else {
                  mu_age <- pac@mu_(pac@params, ages)
                  S_or_mu[idx, ] <- matrix(rep(mu_age, each = length(idx)), nrow = length(idx))
                }
              }
              colnames(S_or_mu) <- as.character(ages)
              return(S_or_mu)
            }

            # type == "loglik"
            # Sum per-row log-likelihood using node-specific PAC
            ll <- rep(NA_real_, nrow(newdata))
            for (nodei in tnodes) {
              idx <- which(newdata$group == nodei)
              if (!length(idx)) next
              pac <- object@AC
              pac@params <- object@par_train[match(nodei, tnodes), ]
              ll[idx] <- logllh(pac,
                                par = pac@params,
                                x = newdata$birthage[idx],
                                r = newdata$censoredage[idx],
                                i = newdata$DthIndicator[idx])
            }

            if (type == "group") return( newdata$group)
            if (type == "par")   return( object@par_train[match(newdata$group, tnodes), , drop = FALSE])

            names(ll) <- rownames(newdata)
            ll
          })


#' @export
setMethod("evaluate", "STPModel",
          function(object, newdata = NULL,
                   metrics = c("AIC","loglik_sum","loglik_vec","IBS"),
                   age_range = 80:115, normalize_IBS = TRUE, ...) {
            if (is.null(newdata)) {
              newdata <- object@STModel@IMM@data_test
              if (is.null(newdata$group)) {
                newdata$group <- predict(object@STModel@model %>% partykit::as.party(),
                                         newdata = newdata, type = "node")
              }
            }
            ages <- sort(as.numeric(age_range))
            out <- list()

            # AIC not defined for STP (no single parametric likelihood)
            if ("AIC" %in% metrics) out$AIC <- NA_real_

            # loglik
            if ("loglik_vec" %in% metrics || "loglik_sum" %in% metrics) {
              ll_vec <- predict(object, newdata, type = "loglik")
              if ("loglik_vec" %in% metrics) out$loglik_vec <- ll_vec
              if ("loglik_sum" %in% metrics) out$loglik_sum <- sum(ll_vec, na.rm = TRUE)
            }

            # IBS via predict-based method

            if ("IBS" %in% metrics) {
              S_pred <- predict(object, newdata, type = "survival", age_range = ages)  # n x |ages|
              ibs <- .ibs_from_Spred(newdata, S_pred, ages, normalize = normalize_IBS)
              out$IBS <- ibs$IBS
              out$IBS_by_age <- ibs$IBS_by_age
            }

            out
          })
