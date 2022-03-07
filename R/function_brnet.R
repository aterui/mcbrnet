#' Generate a random branching network
#'
#' @param n_patch numeric value indicating the number of patches in a network.
#' @param p_branch numeric value indicating the branching probability (success probability of a geometric distribution).
#' @param mean_env_source numeric value indicating the mean value of environmental condition at upstream terminals.
#' @param sd_env_source numeric value indicating the SD of environmental condition at upstream terminals.
#' @param rho numeric value indicating the strength of spatial autocorrelation in environmental condition. The environmental condition at patch x \eqn{z}\out{<sub>x</sub>} is determined as \eqn{z}\out{<sub>x</sub>}\eqn{ = \rho}z\out{<sub>x-1</sub>}\eqn{ + \epsilon}\out{<sub>x</sub>}, where \eqn{\epsilon}\out{<sub>x</sub>} is the random variable drawn from a normal distribution with mean 0 and SD \eqn{\sigma}\out{<sub>env</sub>}. See \href{https://github.com/aterui/mcbrnet}{github page} for further details.
#' @param sd_env_lon numeric value indicating the SD of longitudinal environmental noise.
#' @param mean_disturb_source numeric value indicating the mean of disturbance strength at headwaters. The value is assumed to represent the proportional mortality (0 - 1.0) at the patch level.
#' @param sd_disturb_source numeric value indicating the SD of disturbance strength at headwaters. The SD is defined in a logit scale with a normal distribution.
#' @param sd_disturb_lon numeric value indicating the SD of longitudinal noise of disturbance strength. The SD is defined in a logit scale with a normal distribution.
#' @param asymmetry_factor numeric value rescaling upstream distance. If \code{asymmetry_factor = 1}, distance from a downstream patch x to another upstream patch y would be multiplied by the factor of \code{asymmetry_factor}. This argument does not affect separation distance from an upstream patch to another downstream patch. Default \code{asymmetry_factor = 1} (no asymmetry).
#' @param randomize_patch logical indicating whether randomize patches or not. If \code{FALSE}, the function may generate a biased network with ordered patches. Default \code{TRUE}.
#' @param plot logical indicating if a plot should be shown or not. If \code{FALSE}, a plot of the generated network will not be shown. Default \code{TRUE}.
#' @param patch_label character string indicating a type of patch (vertex) label (either \code{"patch", "branch", "n_upstream"}). \code{"patch"} shows patch ID, \code{"branch"} branch ID, and \code{"n_upstream"} the number of upstream contributing patches. If \code{NULL}, no label will be shown on patches in the plot. Default \code{NULL}.
#' @param patch_size patch (vertex) size in the plot. Default 6.
#' @param patch_scaling logical. If \code{TRUE}, patch (vertex) size will be proportional to the number of upstream contributing patches. The patch (vertex) size will be equal to \code{0.3 * scale_factor} at the upstream terminals and \code{1.3 * scale_factor} at the root. Overrides \code{patch_size}.
#' @param scale_factor numeric value scaling patch (vertex) size. Enabled if \code{patch_scaling = TRUE}.
#' @param n_patch_free logical value indicating whether imposing a constraint on \code{n_patch}. If \code{TRUE}, the number of patches a random variable following a negative binomial distribution.
#'
#' @return \code{adjacency_matrix} adjacency matrix for the generated network.
#' @return \code{distance_matrix} distance matrix for the generated network.
#' @return \code{weighted_distance_matrix} weighted distance matrix for the generated network.
#' @return \code{df_patch} data frame containing patch attributes.
#'
#' @importFrom dplyr %>%
#' @importFrom grDevices grey
#' @importFrom graphics par text
#' @importFrom stats complete.cases rbinom rgeom rnorm runif
#'
#' @section Reference: see \href{https://github.com/aterui/mcbrnet}{github page} for instruction
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @examples
#' brnet(n_patch = 100, p_branch = 0.5)
#'
#' @export
#'
brnet <- function(n_patch,
                  p_branch,
                  mean_env_source = 0,
                  sd_env_source = 1,
                  rho = 1,
                  sd_env_lon = 0.1,
                  mean_disturb_source = 0.9,
                  sd_disturb_source = 1,
                  sd_disturb_lon = 0.1,
                  asymmetry_factor = 1,
                  randomize_patch = TRUE,
                  plot = TRUE,
                  patch_label = NULL,
                  patch_size = 6,
                  patch_scaling = TRUE,
                  scale_factor = 8,
                  n_patch_free = FALSE) {

  # define functions and variables ------------------------------------------

  resample <- function(x, ...) x[sample.int(length(x), ...)]

  if (p_branch > 0 & p_branch < 1) {

    repeat {

      n_branch <- rbinom(n = 1, size = n_patch, prob = p_branch)
      if (n_branch %% 2 == 1) break

    }

  } else {

    if (p_branch == 0) n_branch <- 1

    if (p_branch == 1) {

      if (n_patch %% 2 == 0) stop("n_patch must be an odd number when p_branch = 1")

      n_branch <- n_patch

    }

  }


  # adjacency matrix: linear network ----------------------------------------

  if (p_branch == 0 | n_branch == 1) {

    v_n_patch_branch <- n_patch
    m_adj <- matrix(0, nrow = n_patch, ncol = n_patch)
    x <- seq_len(n_patch - 1)
    y <- 2:n_patch
    m_adj[cbind(x, y)] <- 1
    m_adj[cbind(y, x)] <- 1

  }


  # adjacency matrix: branched network --------------------------------------

  if (p_branch > 0 & n_branch > 1) {

    # vector of the number of patches in each branch
    if (n_patch_free == FALSE) {
      repeat{

        repeat{

          v_n_patch_branch <- rgeom(n = n_branch, prob = p_branch) + 1
          if (sum(v_n_patch_branch) >= n_patch) break

        }

        if (sum(v_n_patch_branch) == n_patch) break

      }
    } else {

      v_n_patch_branch <- rgeom(n = n_branch, prob = p_branch) + 1
      n_patch <- sum(v_n_patch_branch)

    }

    # start_id, end_id, and neighbor list for each branch
    v_end_id <- cumsum(v_n_patch_branch)
    v_start_id <- v_end_id - (v_n_patch_branch - 1)

    list_neighbor_inbranch <- lapply(seq_len(n_branch),
                                     function(i) {
                                       fun_adj(n = v_n_patch_branch[i],
                                               start_id = v_start_id[i])
                                       }
                                     )

    m_neighbor_inbranch <- do.call(rbind, list_neighbor_inbranch)

    # combine parent and offspring branches at confluences
    if (n_branch == 3) {

      parent <- c(1, 1)
      offspg <- c(2, 3)
      m_po <- cbind(parent, offspg)

      list_confluence <- lapply(seq_len(nrow(m_po)),
                                function(x) cbind(v_end_id[m_po[x, 1]],
                                                  v_start_id[m_po[x, 2]]))

      m_confluence <- do.call(rbind, list_confluence)
      m_confluence <- rbind(m_confluence, m_confluence[, c(2, 1)])

    } else {

      n_confluence <- 0.5 * (n_branch - 1)
      v_parent_branch <- seq_len(n_confluence)
      v_offspg_branch <- 2:n_branch

      m_offspg <- matrix(NA,
                         nrow = 2,
                         ncol = n_confluence)

      for (i in n_confluence:1) {

        v_y <- resample(v_offspg_branch[v_offspg_branch > v_parent_branch[i]],
                        size = 2)
        v_offspg_branch <- setdiff(v_offspg_branch, v_y)
        m_offspg[, i] <- v_y

      }

      parent <- rep(v_parent_branch, each = 2)
      offspg <- c(m_offspg)
      m_po <- cbind(parent, offspg)

      list_confluence <- lapply(seq_len(nrow(m_po)),
                                function(x) cbind(v_end_id[m_po[x, 1]],
                                                  v_start_id[m_po[x, 2]]))

      m_confluence <- do.call(rbind, list_confluence)
      m_confluence <- rbind(m_confluence, m_confluence[, c(2, 1)])

    }

    m_neighbor_patch <- rbind(m_neighbor_inbranch, m_confluence)
    m_neighbor_patch <- m_neighbor_patch[complete.cases(m_neighbor_patch), ]

    m_adj <- matrix(0, nrow = n_patch, ncol = n_patch)
    m_adj[m_neighbor_patch] <- 1

  }


  # distance matrix ---------------------------------------------------------

  ## exported function: see "function_adjtodist.R"
  m_distance <- adjtodist(m_adj)


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

  if (!(mean_disturb_source <= 1 & mean_disturb_source >= 0)) {

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


  # asymmetry ---------------------------------------------------------------

  v_distance_to_root <- m_distance[1, ]

  df_asymmetry <- expand.grid(patch1 = seq_len(n_patch),
                              patch2 = seq_len(n_patch)) %>%
    dplyr::mutate(d1 = v_distance_to_root[.data$patch1],
                  d2 = v_distance_to_root[.data$patch2]) %>%
    dplyr::mutate(delta = .data$d2 - .data$d1,
                  distance = m_distance[cbind(.data$patch1, .data$patch2)]) %>%
    dplyr::mutate(x = 0.5 * (.data$distance - .data$delta)) %>%
    dplyr::mutate(weighted_distance = (1 + asymmetry_factor) * .data$x +
                                       asymmetry_factor * .data$delta) %>%
    dplyr::filter(.data$patch1 != .data$patch2)

  m_weighted_distance <- m_distance
  m_weighted_distance[cbind(df_asymmetry$patch1, df_asymmetry$patch2)] <- df_asymmetry$weighted_distance


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
      m_weighted_distance <- m_weighted_distance[df_id$patch, df_id$patch]

    } else {

      df_id <- dplyr::tibble(branch = 1,
                             patch = seq_len(n_patch))

    }
  } else {

    df_id <- dplyr::tibble(branch = branch,
                           patch = patch)

  }


  # visualization -----------------------------------------------------------

  if (plot == TRUE) {

    adj <- igraph::graph.adjacency(adjmatrix = m_adj,
                                   mode = "undirected")

    colvalue <- data.frame(color = viridis::viridis(n_patch, alpha = 0.6),
                           value = sort(v_env))

    layout_tree <- igraph::layout_as_tree(adj,
                                          root = 1,
                                          flip.y = F)

    if (is.null(patch_label)) {

      vertex_label <- NA

    } else {

      if (patch_label == "patch") vertex_label <- seq_len(n_patch)
      if (patch_label == "branch") vertex_label <- df_id$branch
      if (patch_label == "n_upstream") vertex_label <- v_wa

      if (!(patch_label %in% c("patch", "branch", "n_upstream"))) {

        stop("patch_label must be either patch, branch, or n_upstrem")

      }

    }

    if (patch_scaling == TRUE) {

      vertex_size <- I(scale(v_wa,
                             center = min(v_wa),
                             scale = max(v_wa) - min(v_wa)) + 0.3) * scale_factor

    } else {

      vertex_size <- patch_size

    }

    par(mar = c(5.1, 8, 4.1, 2.1))
    igraph::plot.igraph(adj, layout = layout_tree,
                        vertex.size = vertex_size,
                        vertex.label = vertex_label,
                        vertex.label.cex = 0.8,
                        vertex.label.dist = 0.8,
                        vertex.label.degree = pi,
                        vertex.label.color = grey(0.3),
                        vertex.frame.color = grey(0.5),
                        vertex.color = colvalue$color[match(v_env, colvalue$value)],
                        edge.width = 1.8,
                        edge.color = "steelblue")

    plotfunctions::gradientLegend(valRange = range(v_env),
                                  color = viridis::viridis(n_patch),
                                  pos = 0.8,
                                  side = 2,
                                  dec = 2)

    pc <- c(plotfunctions::getCoords(0, side = 1),
            plotfunctions::getCoords(1, side = 2))

    text(x = pc[1],
         y = pc[2],
         labels = "Environmental value",
         adj = 1)

  }


  # return ------------------------------------------------------------------

  return(list(adjacency_matrix = m_adj,
              distance_matrix = m_distance,
              weighted_distance_matrix = m_weighted_distance,
              df_patch = dplyr::tibble(patch_id = seq_len(n_patch),
                                       branch_id = as.numeric(df_id$branch),
                                       environment = c(v_env),
                                       disturbance = c(boot::inv.logit(v_disturb)),
                                       n_patch_upstream = c(v_wa))))
}
