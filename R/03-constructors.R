#' @export
InitPAC <- function(name, params, survival_, mu_, cumhaz_, bounds) {
  methods::new("PAC",
               name = name, params = params,
               survival_ = survival_, mu_ = mu_, cumhaz_ = cumhaz_, bounds = bounds
  )
}

# Defaults & family constructors
DEFAULT_PARAMS_GOMPERTZ <- c(-10, 0.1)
DEFAULT_BOUNDS_GOMPERTZ <- list(lower = c(-20, 0.0001), upper = c(20, 1))
#' @export
GompertzPAC <- function(params = DEFAULT_PARAMS_GOMPERTZ,
                        bounds = DEFAULT_BOUNDS_GOMPERTZ) {
  InitPAC("Gompertz", params, survival_gompertz, forceofmortality_gompertz, cumulative_hazard_gompertz, bounds)
}

DEFAULT_PARAMS_MAKEHAM <- c(-10, 0.1, 0)
DEFAULT_BOUNDS_MAKEHAM <- list(lower = c(-20, 0.0001, 0), upper = c(20, 1, 1))
#' @export
MakehamPAC <- function(params = DEFAULT_PARAMS_MAKEHAM,
                       bounds = DEFAULT_BOUNDS_MAKEHAM) {
  InitPAC("Gompertz-Makeham", params, survival_makeham, forceofmortality_makeham, cumulative_hazard_makeham, bounds)
}

DEFAULT_PARAMS_PERKS <- c(-15, 0.1, 0.7, 0)
DEFAULT_BOUNDS_PERKS <- list(lower = c(-20, 0.0001, -1, 0), upper = c(20, 1, 5, 1))
#' @export
PerksPAC <- function(params = DEFAULT_PARAMS_PERKS,
                     bounds = DEFAULT_BOUNDS_PERKS) {
  InitPAC("Perks", params, survival_perks, forceofmortality_perks, cumulative_hazard_perks, bounds)
}

DEFAULT_PARAMS_BEARD <- c(-15, 0.1, 0.7)
DEFAULT_BOUNDS_BEARD <- list(lower = c(-20, 0.0001, -1), upper = c(20, 1, 5))
#' @export
BeardPAC <- function(params = DEFAULT_PARAMS_BEARD,
                     bounds = DEFAULT_BOUNDS_BEARD) {
  InitPAC("Beard", params, survival_beard, forceofmortality_beard, cumulative_hazard_beard, bounds)
}

DEFAULT_PARAMS_KANNISTO <- c(-10, 0.1)
DEFAULT_BOUNDS_KANNISTO <- list(lower = c(-20, 0.0001), upper = c(20, 1))
#' @export
KannistoPAC <- function(params = DEFAULT_PARAMS_KANNISTO,
                        bounds = DEFAULT_BOUNDS_KANNISTO) {
  InitPAC("Kannisto", params, survival_kannisto, forceofmortality_kannisto, cumulative_hazard_kannisto, bounds)
}

DEFAULT_PARAMS_THATCHER <- c(-10, 0.1, 0)
DEFAULT_BOUNDS_THATCHER <- list(lower = c(-20, 0.0001, 0), upper = c(20, 1, 1))
#' @export
ThatcherPAC <- function(params = DEFAULT_PARAMS_THATCHER,
                        bounds = DEFAULT_BOUNDS_THATCHER) {
  InitPAC("Thatcher", params, survival_thatcher, forceofmortality_thatcher, cumulative_hazard_thatcher, bounds)
}
