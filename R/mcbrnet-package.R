#' @useDynLib mcbrnet, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"

.onUnload = function(libpath) {
  library.dynam.unload("mcbrnet", libpath)
}

## usethis namespace: start
## usethis namespace: end
NULL
