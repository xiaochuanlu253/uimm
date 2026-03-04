# Registry
.pac_registry <- list(
  Gompertz  = GompertzPAC,
  Makeham   = MakehamPAC,
  Perks     = PerksPAC,
  Beard     = BeardPAC,
  Kannisto  = KannistoPAC,
  Thatcher  = ThatcherPAC
)

#' Fit a PAC model via a formula
#' @param formula Surv(birthage, censoredage, DthIndicator) ~ 1
#' @param data data.frame
#' @param AC character in names of .pac_registry
#' @param start,bounds optional overrides
#' @param patience,tol,verbose passed to fit()
#' @return A \code{PACFit} object
#'
#' @examplesIf requireNamespace("GenSA", quietly = TRUE)
#' # Minimal, fast example (runs in <1s)
#' set.seed(1)        # \dontshow{set.seed(1)}
#' n <- 120
#' df <- data.frame(
#'   birthage     = runif(n, 80, 90),
#'   censoredage  = runif(n, 90, 100),
#'   DthIndicator = rbinom(n, 1, 0.4)
#' )
#'
#' fit <- pac_fit(
#'   Surv(birthage, censoredage, DthIndicator) ~ 1,
#'   data = df, AC = "Gompertz", patience = 30, tol = 1e-6
#' )
#' fit
#' AIC(fit)
#'
#' @export
pac_fit <- function(formula, data,
                    AC = c("Gompertz","Makeham","Perks","Beard","Kannisto","Thatcher"),
                    start = NULL, bounds = NULL,
                    patience = 500, tol = 1e-6, verbose = FALSE) {
  AC <- match.arg(AC)
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  y  <- stats::model.response(mf)
  sType <- attr(y, "type")
  if (!(sType %in% c("counting"))) stop("Use counting Surv: Surv(birthage, censoredage, DthIndicator) ~ 1")

  if (ncol(y) == 3){
    df_train <- data.frame(
      birthage = as.numeric(y[,1]),
      censoredage = as.numeric(y[,2]),
      DthIndicator = as.integer(y[,3])
    )
  } else if (ncol(y) == 2){
    df_train <- data.frame(
      birthage = rep(0, nrow(y)),
      censoredage = as.numeric(y[,1]),
      DthIndicator = as.integer(y[,2])
    )
  } else {
    stop("Unexpected Surv response format with ", ncol(y), " columns.")
  }

  ctor <- .pac_registry[[AC]]
  if (is.null(ctor)) stop("Unknown AC family: ", AC)
  pac_model <- ctor()
  if (!is.null(start))  pac_model@params <- start
  if (!is.null(bounds)) pac_model@bounds <- bounds

  fit_obj <- methods::new("PACFit", model = pac_model, data = list(train = df_train, test = NULL))
  fit(fit_obj, data_train = df_train, patience = patience, tol = tol, verbose = verbose)
}
