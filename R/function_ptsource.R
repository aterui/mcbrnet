#' Create point source disturbance
#'
#' @param x `brnet()` object
#' @param n_source Number of point sources
#' @param p Reduction factor for downstream propagation
#' @param q Reduction factor for upstream propagation
#' @param pattern Spatial pattern of point sources. "random", "cluster", "upstream", "downstream"
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @section Reference: see \href{https://aterui.github.io/mcbrnet/}{package webpage} for instruction
#'
#' @importFrom dplyr %>% filter
#' @importFrom rlang .data
#'
#' @export
#'

ptsource <- function(x,
                     n_source,
                     p,
                     q,
                     pattern = "random") {


  # parameter check ---------------------------------------------------------

  if (p < 0 | p > 1 | q < 0 | q > 1) stop("p & q must be 0 - 1")
  if (nrow(x$distance_matrix) < n_source) stop("n_source must be < n_patch")
  if (!inherits(x, what = "brnet")) stop("x must be a 'brnet' object")


  # prepare basic objects ---------------------------------------------------

  ## basic objects
  v_patch_order <- x$distance_matrix[1, ]
  n_patch <- nrow(x$adjacency_matrix)

  v_d_to_root <- sort(x$distance_matrix[1, ])
  m_adj <- x$adjacency_matrix[order(v_patch_order), order(v_patch_order)]
  m_dist <- x$distance_matrix[order(v_patch_order), order(v_patch_order)]
  v_wa <- x$df_patch$n_patch_upstream[order(v_patch_order)]

  # extract upstream adjacency
  m_adj_up <- m_adj
  m_adj_up[lower.tri(m_adj_up)] <- 0

  # upstream tributary proportion
  m_wa_prop <- t(apply(X = m_adj_up,
                       MARGIN = 1,
                       FUN = function(k) {
                         (k * v_wa) / ifelse(sum(k) == 0,
                                             yes = 1,
                                             no = sum(k * v_wa))
                       }))


  # impact ------------------------------------------------------------------

  if (pattern == "random") {
    pts <- resample(seq_len(n_patch), size = n_source)
  }

  if (pattern == "cluster") {
    s <- resample(seq_len(n_patch), size = 1)
    distance <- m_dist[s, ]
    distance[distance == 0] <- 2
    prob <- 1 / dis
    pts <- resample(seq_len(n_patch),
                    size = n_source,
                    prob = prob)
  }

  if (pattern == "downstream") {
    prob <- v_wa
    pts <- resample(seq_len(n_patch),
                    size = n_source,
                    prob = prob)
  }

  if (pattern == "upstream") {
    prob <- 1 / v_wa
    pts <- resample(seq_len(n_patch),
                    size = n_source,
                    prob = prob)
  }

  v_d <- v_u <- v_d_attr <- v_u_attr <- rep(0, n_patch)
  v_d[pts] <- v_u[pts] <- 1
  v_d_attr[pts] <- v_u_attr[pts] <- 1

  for (i in seq_len(max(v_d_to_root))) {
    ## upstream propagation
    v_u <- t(m_wa_prop) %*% (q * v_u)
    v_u_attr <- v_u + v_u_attr

    ## downstream propagation
    v_d <- m_wa_prop %*% (p * v_d)
    v_d_attr <- v_d + v_d_attr
  }

  v_u_attr[pts] <- v_u_attr[pts] - 1
  v_attr <- v_u_attr + v_d_attr


  # export ------------------------------------------------------------------

  v_source <- rep(0, n_patch)
  v_source[pts] <- 1

  df0 <- dplyr::tibble(patch_id = x$df_patch$patch_id[order(v_patch_order)],
                       impact = c(v_attr),
                       point_source = v_source) %>%
    dplyr::arrange(.data$patch_id) %>%
    dplyr::right_join(y = x$df_patch,
                      by = "patch_id")

  x$df_patch <- df0

  return(x)

}
