#' Internal function: get fcl
#'
#' @param x Community matrix
#' @param delta Preference to basal over ig-prey
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

fun_get_fcl <- function(x, delta) {

  if (ncol(x) != length(delta)) stop("invalid dimension: error in x or delta")

  # raw fcl - 1, 2, 3
  v_fcl_raw <- apply(X = x,
                     MARGIN = 2,
                     FUN = function(y) {
                         fcl <- ifelse(y[1] > 0,
                                       yes = max(which(y > 0)),
                                       no = 0)
                     })
  # omnivory
  omn <- ifelse(v_fcl_raw == 3,
                yes = 1,
                no = 0)

  v_fcl <- (1 - omn) * v_fcl_raw + omn * (1 * delta + 2 * (1 - delta) + 1)

  return(v_fcl)
}
