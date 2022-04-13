#' Internal function: resample
#'
#' @param x Vector
#' @param ... Additional arguments passed to \code{sample}
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

resample <- function(x, ...) x[sample.int(length(x), ...)]
