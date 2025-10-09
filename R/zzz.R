#' @keywords internal
"_PACKAGE"

# Param cache env (no global <<- or disk writes)
.pacsurv_cache <- new.env(parent = emptyenv())

.onUnload <- function(libpath) {
  rm(list = ls(.pacsurv_cache), envir = .pacsurv_cache)
}
