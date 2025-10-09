#' Unified formula-first model fitters
#'
#' These helpers provide a consistent front door for fitting baseline
#' parametric age-only models and their tree/partitioning variants.
#'
#' @name mortality_wrappers
NULL

# ---- Internals -------------------------------------------------------------

#' @keywords internal
.AC_registry <- list(
  Gompertz  = GompertzPAC,
  Makeham   = MakehamPAC,
  Perks     = PerksPAC,
  Beard     = BeardPAC,
  Kannisto  = KannistoPAC,
  Thatcher  = ThatcherPAC
)

#' Resolve an AC family specification to a PAC object
#' @param AC a `PAC` object, a family name, or a family constructor
#' @keywords internal
AC_from <- function(AC = "Beard") {
  if (inherits(AC, "PAC")) return(AC)
  if (is.character(AC)) {
    ctor <- .AC_registry[[AC]]
    if (is.null(ctor)) stop("Unknown AC family: ", AC, call. = FALSE)
    return(ctor())
  }
  if (is.function(AC)) return(AC())
  stop("AC must be a PAC object, a family name, or a family constructor.", call. = FALSE)
}

#' Guard: baseline is age-only (~ 1)
#' @keywords internal
.ensure_age_only <- function(formula) {
  tl <- attr(stats::terms(formula), "term.labels")
  if (length(tl) > 0) {
    warning("Baseline models are age-only; RHS covariates will be ignored.")
  }
}

# ---- Public wrappers -------------------------------------------------------

#' Fit a baseline age-only parametric model (PAC)
#'
#' @param formula \code{Surv(birthage, censoredage, DthIndicator) ~ 1}
#' @param data,data_train,data_test Training/test data. If \code{data} is a
#'   data frame it is treated as training data.
#' @param AC AC family (\code{"Beard"}, \code{"Gompertz"}, …) or a \code{PAC} object.
#' @param patience,tol,verbose Passed to \code{fit()}.
#' @return A \code{BaselineModel} object.
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   birthage     = runif(100, 80, 90),
#'   censoredage  = runif(100, 90, 105),
#'   DthIndicator = rbinom(100, 1, 0.4)
#' )
#' \donttest{
#' if (requireNamespace("GenSA", quietly = TRUE)) {
#'   b0 <- baseline_fit(Surv(birthage, censoredage, DthIndicator) ~ 1,
#'                      data = df, AC = "Gompertz", patience = 30)
#'   b0
#' }
#' }
#' @export
baseline_fit <- function(formula,
                         data = NULL,
                         data_train = NULL,
                         data_test = NULL,
                         AC = "Beard",
                         patience = 500, tol = 1e-6, verbose = FALSE) {
  .ensure_age_only(formula)

  # resolve training frame
  train_df <- if (is.data.frame(data)) {
    data
  } else if (is.list(data) && !is.null(data$train)) {
    data$train
  } else if (is.data.frame(data_train)) {
    data_train
  } else {
    stop("Provide training data via `data` (data.frame) or `data_train`.", call. = FALSE)
  }

  # pac_fit expects an AC name; resolve if a PAC or constructor is supplied
  AC_name <- if (inherits(AC, "PAC")) {
    AC@name
  } else if (is.function(AC)) {
    AC()@name
  } else {
    as.character(AC)
  }

  pf <- pac_fit(formula, data = train_df,
                AC = AC_name, patience = patience, tol = tol, verbose = verbose)

  IMM <- MortalityModel(data_train = train_df, data_test = data_test)
  methods::new("BaselineModel", IMM = IMM, model_fit = pf)
}

#' Fit an LTRC survival tree (structure only)
#'
#' @param formula A survival formula.
#' @param data,data_train,data_test See \code{baseline_fit()}.
#' @return An \code{STModel} object.
#' @export
st_fit <- function(formula,
                   data = NULL,
                   data_train = NULL,
                   data_test = NULL) {
  STModel(formula = formula, data = if (is.null(data)) data_train else data)
}

#' Fit a survival tree with per-node PAC models (STP)
#'
#' @param formula A survival formula.
#' @param data,data_train,data_test See \code{baseline_fit()}.
#' @param AC AC family or PAC object for node-wise fits.
#' @return An \code{STPModel} object.
#' @export
stp_fit <- function(formula,
                    data = NULL,
                    data_train = NULL,
                    data_test = NULL,
                    AC = "Beard") {
  st <- STModel(formula = formula, data = if (is.null(data)) data_train else data)
  STPModel(ST = st, AC = AC_from(AC))
}

#' Fit a model-based recursive partitioning model (mob + Gompertz)
#'
#' @param formula A survival formula with partitioning variables on the RHS via \code{|}.
#' @param data,data_train,data_test See \code{baseline_fit()}.
#' @return An \code{MTPModel} object, or a \code{BaselineModel} fallback if only one node results.
#' @export
mtp_fit <- function(formula,
                    data = NULL,
                    data_train = NULL,
                    data_test = NULL) {
  MTPModel(formula = formula, data = if (is.null(data)) data_train else data)
}

#' Fit a Gompertz proportional hazards model via \pkg{flexsurv}
#'
#' @param formula A standard survival formula with covariates.
#' @param data,data_train,data_test See \code{baseline_fit()}.
#' @return A \code{GPHModel} object.
#' @export
gph_fit <- function(formula,
                    data = NULL,
                    data_train = NULL,
                    data_test = NULL) {
  GPHModel(formula = formula, data = if (is.null(data)) data_train else data)
}

#' Unified front door for all supported model types
#'
#' @param model One of \code{"baseline"}, \code{"st"}, \code{"stp"}, \code{"mtp"}, \code{"gph"}.
#' @param formula A survival formula.
#' @param data,data_train,data_test See \code{baseline_fit()}.
#' @param ... Passed through to the respective \code{*_fit()} function.
#' @return An object of class \code{BaselineModel}, \code{STModel}, \code{STPModel},
#'   \code{MTPModel}, or \code{GPHModel}.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   birthage     = runif(80, 80, 90),
#'   censoredage  = runif(80, 90, 105),
#'   DthIndicator = rbinom(80, 1, 0.4)
#' )
#' \donttest{
#' if (requireNamespace("GenSA", quietly = TRUE)) {
#'   m0 <- mort_fit("baseline",
#'     formula = Surv(birthage, censoredage, DthIndicator) ~ 1,
#'     data = df, AC = "Gompertz", patience = 25)
#'   m0
#' }
#' }
#' @export
mort_fit <- function(model = c("baseline","st","stp","mtp","gph"),
                     formula,
                     data = NULL,
                     data_train = NULL,
                     data_test = NULL,
                     ...) {
  # If user passed a plain data.frame via `data`, treat as training data
  if (is.data.frame(data)) {
    data_train <- data
    data_test  <- if (is.list(data) && !is.null(data$test)) data$test else data_test
  }

  model <- match.arg(model)
  switch(model,
         baseline = baseline_fit(formula, data = data, data_train = data_train, data_test = data_test, ...),
         st       = st_fit      (formula, data = data, data_train = data_train, data_test = data_test),
         stp      = stp_fit     (formula, data = data, data_train = data_train, data_test = data_test, ...),
         mtp      = mtp_fit     (formula, data = data, data_train = data_train, data_test = data_test),
         gph      = gph_fit     (formula, data = data, data_train = data_train, data_test = data_test),
         stop("Unknown model type.", call. = FALSE)
  )
}
