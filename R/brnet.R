#' Generate a random branching network
#'
#' @param n_patch Number of patches in a network.
#' @param p_branch Branching probability (success probability of a geometric distribution).
#' @param mean_env_source Mean value of environmental condition at upstream terminals.
#' @param sd_env_source SD of environmental condition at upstream terminals.
#' @param rho Strength of spatial autocorrelation in environmental condition.
#' @param sd_env_lon SD of longitudinal environmental noise.
#' @param mean_disturb_source Mean disturbance strength at headwaters. The value is assumed to represent the proportional mortality (0 - 1.0) at the patch level.
#' @param sd_disturb_source SD of disturbance strength at headwaters. The SD is defined in a logit scale with a normal distribution.
#' @param sd_disturb_lon SD of longitudinal noise of disturbance strength. The SD is defined in a logit scale with a normal distribution.
#' @param randomize_patch Whether randomize patches or not. If \code{FALSE}, the function may generate a biased network with ordered patches. Default \code{TRUE}.
#' @param plot Whether a plot should be shown or not. If \code{FALSE}, a plot of the generated network will not be shown. Default \code{TRUE}.
#' @param patch_color Type of patch (vertex) label (either \code{"env"}, \code{"disturbance"} or any color code). Default \code{"env"}.
#' @param patch_label Type of patch (vertex) label (either \code{"patch", "branch", "n_upstream"}). \code{"patch"} shows patch ID, \code{"branch"} branch ID, and \code{"n_upstream"} the number of upstream contributing patches. If \code{"none"}, no label will be shown on patches in the plot. Default \code{"none"}.
#' @param patch_size Patch (vertex) size in the plot.
#' @param n_patch_free Whether imposing a constraint on \code{n_patch}. If \code{TRUE}, the number of patches a random variable following a negative binomial distribution.
#' @param ... Arguments passed to \code{ggraph::geom_node_label()}.
#'
#' @importFrom dplyr %>%
#' @importFrom grDevices grey
#' @importFrom graphics par text
#' @importFrom stats complete.cases rbinom rgeom rnorm runif
#'
#' @section Reference: see \href{https://aterui.github.io/mcbrnet/}{package webpage} for instruction
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @examples
#' \dontrun{
#' # not run
#' brnet(n_patch = 10, p_branch = 0.5)
#' }
#'
#' @export
#'

brnet <- function(n_patch = 50,
                  p_branch = 0.5,
                  mean_env_source = 0,
                  sd_env_source = 1,
                  rho = 1,
                  sd_env_lon = 0.1,
                  mean_disturb_source = 0.9,
                  sd_disturb_source = 1,
                  sd_disturb_lon = 0.1,
                  randomize_patch = TRUE,
                  plot = FALSE,
                  patch_color = "env",
                  patch_label = "none",
                  patch_size = 3,
                  n_patch_free = FALSE,
                  ...) {


  # define variables --------------------------------------------------------

  ## internal function: see "fun_get_n_branch.R"
  n_branch <- fun_get_n_branch(n_patch = n_patch,
                               p_branch = p_branch)


  # adjacency matrix --------------------------------------------------------

  list_adj <- fun_m_adj(n_patch = n_patch,
                        p_branch = p_branch,
                        n_branch = n_branch,
                        n_patch_free = n_patch_free)

  v_n_patch_branch <- list_adj$v_n_patch_branch
  m_adj <- list_adj$m_adj


  # distance matrix ---------------------------------------------------------

  ## internal function: see "function_adjtodist.R"
  m_distance <- m_adj %>%
    igraph::graph_from_adjacency_matrix("undirected") %>%
    igraph::distances()


  # upstream watershed area ------------------------------------------------

  ## internal function: see "fun_wa.R"
  m_wa <- fun_wa(x = m_adj)
  v_wa <- rowSums(m_wa)


  # environmental condition -------------------------------------------------

  ## internal function: see "fun_patch_attr.R"
  v_env <- fun_patch_attr(x = m_adj,
                          n_branch = n_branch,
                          mean_source = mean_env_source,
                          sd_source = sd_env_source,
                          sd_lon = sd_env_lon,
                          m_distance = m_distance,
                          rho = rho,
                          v_wa = v_wa,
                          logit = FALSE)


  # disturbance -------------------------------------------------------------

  if (!(mean_disturb_source <= 1 && mean_disturb_source >= 0)) {

    stop("mean_disturb_source must be between 0 and 1")

  }

  ## internal function: see "fun_patch_attr.R"
  v_disturb <- fun_patch_attr(x = m_adj,
                              n_branch = n_branch,
                              mean_source = mean_disturb_source,
                              sd_source = sd_disturb_source,
                              sd_lon = sd_disturb_lon,
                              m_distance = m_distance,
                              rho = 1,
                              v_wa = v_wa,
                              logit = TRUE)


  # randomize nodes ---------------------------------------------------------

  branch <- unlist(lapply(seq_len(n_branch),
                          function(x) rep(x, each = v_n_patch_branch[x])))

  patch <- seq_len(n_patch)

  if (randomize_patch == TRUE) {

    if (n_branch > 1) {

      df_id <- dplyr::tibble(branch = as.character(c(1, resample(2:n_branch)))) %>%
        dplyr::left_join(data.frame(patch,
                                    branch = as.character(branch)),
                         by = "branch")
      v_wa <- v_wa[df_id$patch]
      v_env <- v_env[df_id$patch]
      v_disturb <- v_disturb[df_id$patch]
      m_adj <- m_adj[df_id$patch, df_id$patch]
      m_distance <- m_distance[df_id$patch, df_id$patch]

    } else {

      df_id <- dplyr::tibble(branch = 1,
                             patch = seq_len(n_patch))

    }
  } else {

    df_id <- dplyr::tibble(branch = branch,
                           patch = patch)

  }

  df_patch <- dplyr::tibble(patch_id = seq_len(n_patch),
                            branch_id = as.numeric(df_id$branch),
                            environment = c(v_env),
                            disturbance = c(boot::inv.logit(v_disturb)),
                            n_patch_upstream = c(v_wa))


  # visualization -----------------------------------------------------------

  if (plot == TRUE) {

    g <- ggbrnet(x = list(adjacency_matrix = m_adj,
                          distance_matrix = m_distance,
                          df_patch = df_patch),
                 patch_color = patch_color,
                 patch_label = patch_label,
                 patch_size = patch_size,
                 ...)

    print(g)

  }

  # return ------------------------------------------------------------------

  list_brnet <- list(adjacency_matrix = m_adj,
                     distance_matrix = m_distance,
                     df_patch = df_patch)

  class(list_brnet) <- "brnet"

  return(list_brnet)
}
