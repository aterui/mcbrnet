#' Generate a random branching network
#'
#' @param n_patch numeric value indicating the number of patches in a network.
#' @param p_branch numeric value indicating the branching probability (success probability of a geometric distribution).
#' @param min_env numeric value indicating minimum value of environmental condition at upstream terminals (minimum of a uniform distribution).
#' @param max_env numeric value indicating maximum value of environmental condition at upstream terminals (maximum of a uniform distribution).
#' @param rho numeric value indicating the strength of spatial autocorrelation in environmental condition. The environmental condition at patch x \eqn{z}\out{<sub>x</sub>} is determined as \eqn{z}\out{<sub>x</sub>}\eqn{ = \rho}z\out{<sub>x-1</sub>}\eqn{ + \epsilon}\out{<sub>x</sub>}, where \eqn{\epsilon}\out{<sub>x</sub>} is the random variable drawn from a normal distribution with mean 0 and SD \eqn{\sigma}\out{<sub>env</sub>}. See \href{https://github.com/aterui/mcbrnet}{github page} for further details.
#' @param sd_env numeric value indicating the SD of spatial environmental noise (\eqn{\sigma}\out{<sub>env</sub>}).
#' @param randomize_patch logical indicating whether randomize patches or not. If \code{FALSE}, the function may generate a biased network with ordered patches. Default \code{TRUE}.
#' @param plot logical indicating if a plot should be shown or not. If \code{FALSE}, a plot of the generated network will not be shown. Default \code{TRUE}.
#' @param patch_label character string indicating a type of patch (vertex) label (either \code{"patch", "branch", "n_upstream"}). \code{"patch"} shows patch ID, \code{"branch"} branch ID, and \code{"n_upstream"} the number of upstream contributing patches. If \code{NULL}, no label will be shown on patches in the plot. Default \code{NULL}.
#' @param patch_size patch (vertex) size in the plot. Default 3.
#' @param patch_scaling logical. If \code{TRUE}, patch (vertex) size will be proportional to the number of upstream contributing patches. The patch (vertex) size will be equal to \code{0.3 * scale_factor} at the upstream terminals and \code{1.3 * scale_factor} at the root. Overrides \code{patch_size}.
#' @param scale_factor numeric value scaling patch (vertex) size. Enabled if \code{patch_scaling = TRUE}.
#'
#' @return \code{adjacency_matrix} adjacency matrix for the generated network.
#' @return \code{distance_matrix} distance matrix for the generated network.
#' @return \code{patch_df} a data frame containing patch attributes.
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
                  min_env = -1,
                  max_env = 1,
                  rho = 1,
                  sd_env = 0.1,
                  randomize_patch = TRUE,
                  plot = TRUE,
                  patch_label = NULL,
                  patch_size = 3,
                  patch_scaling = FALSE,
                  scale_factor = 8) {

  # define functions and variables ------------------------------------------

  resample <- function(x, ...) x[sample.int(length(x), ...)]

  fun_adj <- function(n, start_id = 1) {

    if (n == 1) {
      m_y <- cbind(1, NA)
    }

    if (n == 2) {
      m_y <- cbind(c(1, 2), c(2, 1))
    }

    if (n > 2) {
      y1 <- c(1, rep(2:(n - 1), each = 2), n)
      y2 <- c(2, sapply(2:(n - 1), function(i) c(i - 1, i + 1)), n - 1)
      m_y <- cbind(y1, y2)
    }

    colnames(m_y) <- c("patch_1", "patch_2")
    m_y <- m_y + (start_id - 1)
    return(m_y)
  }


  if (p_branch > 0 & p_branch < 1) {
    repeat {
      n_branch <- rbinom(n = 1, size = n_patch, prob = p_branch)
      if (n_branch %% 2 == 1) break
    }
  } else {
    if (p_branch == 0) {
      n_branch <- 1
    }

    if (p_branch == 1) {
      if (n_patch %% 2 == 0) stop("n_patch must be an odd number when p_branch = 1")
      n_branch <- n_patch
    }
  }


  # adjacency matrix: linear network ----------------------------------------

  if (p_branch == 0 | n_branch == 1) {
    v_n_patch_branch <- n_patch
    m_adj <- matrix(0, nrow = n_patch, ncol = n_patch)
    x <- 1:(n_patch - 1)
    y <- 2:n_patch
    m_adj[cbind(x, y)] <- 1
    m_adj[cbind(y, x)] <- 1
  }


  # adjacency matrix: branched network --------------------------------------

  if (p_branch > 0 & n_branch > 1) {

    # vector of the number of patches in each branch
    repeat{
      repeat{
        v_n_patch_branch <- rgeom(n = n_branch, prob = p_branch) + 1
        if (sum(v_n_patch_branch) >= n_patch) break
      }
      if (sum(v_n_patch_branch) == n_patch) break
    }

    # start_id, end_id, and neighbor list for each branch
    v_end_id <- cumsum(v_n_patch_branch)
    v_start_id <- v_end_id - (v_n_patch_branch - 1)
    list_neighbor_inbranch <- lapply(1:n_branch, function(i) fun_adj(n = v_n_patch_branch[i], start_id = v_start_id[i]))
    m_neighbor_inbranch <- do.call(rbind, list_neighbor_inbranch)

    # combine parent and offspring branches at confluences
    if (n_branch == 3) {
      parent <- c(1, 1)
      offspg <- c(2, 3)
      m_po <- cbind(parent, offspg)
      list_confluence <- lapply(seq_len(nrow(m_po)), function(x) cbind(v_end_id[m_po[x, 1]], v_start_id[m_po[x, 2]]))
      m_confluence <- do.call(rbind, list_confluence)
      m_confluence <- rbind(m_confluence, m_confluence[, c(2, 1)])
    } else {
      n_confluence <- 0.5 * (n_branch - 1)
      v_parent_branch <- 1:n_confluence
      v_offspg_branch <- 2:n_branch

      m_offspg <- matrix(NA, nrow = 2, ncol = n_confluence)
      for (i in n_confluence:1) {
        v_y <- resample(v_offspg_branch[v_offspg_branch > v_parent_branch[i]], size = 2)
        v_offspg_branch <- setdiff(v_offspg_branch, v_y)
        m_offspg[, i] <- v_y
      }

      parent <- rep(v_parent_branch, each = 2)
      offspg <- c(m_offspg)
      m_po <- cbind(parent, offspg)

      list_confluence <- lapply(seq_len(nrow(m_po)), function(x) cbind(v_end_id[m_po[x, 1]], v_start_id[m_po[x, 2]]))
      m_confluence <- do.call(rbind, list_confluence)
      m_confluence <- rbind(m_confluence, m_confluence[, c(2, 1)])
    }

    m_neighbor_patch <- rbind(m_neighbor_inbranch, m_confluence)
    m_neighbor_patch <- m_neighbor_patch[complete.cases(m_neighbor_patch), ]

    m_adj <- matrix(0, nrow = n_patch, ncol = n_patch)
    m_adj[m_neighbor_patch] <- 1
  }


  # distance matrix ---------------------------------------------------------

  m_distance <- matrix(0, ncol = n_patch, nrow = n_patch)
  m_identity <- diag(x = 1, nrow = n_patch, ncol = n_patch)

  for (i in 1:n_patch) {
    m_identity <- m_identity %*% m_adj
    m_distance[m_identity != 0 & m_distance == 0] <- i
    if (length(which(m_distance == 0)) == 0) break
  }

  diag(m_distance) <- 0
  m_identity <- NULL


  # upstream watershed area ------------------------------------------------

  m_adj_up <- m_adj
  m_adj_up[lower.tri(m_adj_up)] <- 0
  m_wa <- matrix(0, ncol = n_patch, nrow = n_patch)
  m_identity <- diag(1, nrow = n_patch, ncol = n_patch)

  for (i in 1:n_patch) {
    m_identity <- m_identity %*% m_adj_up
    m_wa[m_identity != 0 & m_wa == 0] <- 1
    if (length(which(m_wa[upper.tri(m_wa)] == 0)) == 0) break
  }

  diag(m_wa) <- 1
  v_wa <- rowSums(m_wa)


  # environmental condition -------------------------------------------------

  m_wa_prop <- t(apply(X = m_adj_up, MARGIN = 1, function(x) (x * v_wa) / ifelse(sum(x) == 0, 1, sum(x * v_wa))))
  n_source <- 0.5 * (n_branch + 1)
  source <- which(rowSums(m_adj_up) == 0)
  v_z_dummy <- v_z <- v_env <- rep(0, n_patch)
  v_z[source] <- v_env[source] <- runif(n_source, min = min_env, max = max_env)
  v_z_dummy[source] <- 1

  if (!(rho <= 1 & rho >= 0)) stop("rho must be between 0 and 1")
  for (i in 1:max(m_distance[1, ])) {
    v_eps <- rep(0, n_patch)
    v_eps[v_z_dummy != 0] <- rnorm(n = length(v_eps[v_z_dummy != 0]), mean = 0, sd = sd_env)
    v_z <- m_wa_prop %*% ((rho * v_z) + v_eps)
    v_z_dummy <- m_wa_prop %*% v_z_dummy
    v_env <- v_z + v_env
  }


  # randomize nodes ---------------------------------------------------------

  branch <- unlist(lapply(1:n_branch, function(x) rep(x, each = v_n_patch_branch[x])))
  patch <- 1:n_patch

  if (randomize_patch == TRUE) {
    if (n_branch > 1) {
      df_id <- dplyr::tibble(branch = as.character(c(1, resample(2:n_branch)))) %>%
        dplyr::left_join(data.frame(patch, branch = as.character(branch)), by = "branch")
      v_wa <- v_wa[df_id$patch]
      v_env <- v_env[df_id$patch]
      m_adj <- m_adj[df_id$patch, df_id$patch]
      m_distance <- m_distance[df_id$patch, df_id$patch]
    } else {
      df_id <- data.frame(branch = 1, patch = 1:n_patch)
    }
  } else {
    df_id <- data.frame(branch = branch, patch = patch)
  }


  # visualization -----------------------------------------------------------

  if (plot == TRUE) {
    adj <- igraph::graph.adjacency(adjmatrix = m_adj, mode = "undirected")
    colvalue <- data.frame(color = viridis::viridis(n_patch, alpha = 0.6), value = sort(v_env))
    layout_tree <- igraph::layout_as_tree(adj, root = 1, flip.y = F)

    if (is.null(patch_label)) {
      vertex_label <- NA
    } else {
      if (patch_label == "patch") vertex_label <- 1:n_patch
      if (patch_label == "branch") vertex_label <- df_id$branch
      if (patch_label == "n_upstream") vertex_label <- v_wa
      if (!(patch_label %in% c("patch", "branch", "n_upstream"))) stop("patch_label must be either patch, branch, or n_upstrem")
    }

    if (patch_scaling == TRUE) {
      vertex_size <- I(scale(v_wa, center = min(v_wa), scale = max(v_wa) - min(v_wa)) + 0.3) * scale_factor
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
    plotfunctions::gradientLegend(valRange = range(v_env), color = viridis::viridis(n_patch),
                                  pos = 0.8, side = 2, dec = 2)
    pc <- c(plotfunctions::getCoords(0, side = 1), plotfunctions::getCoords(1, side = 2))
    text(x = pc[1], y = pc[2], labels = "Environmental value", adj = 1)
  }


  # return ------------------------------------------------------------------

  rownames(m_adj) <- colnames(m_adj) <- sapply(seq_len(n_patch), function(x) paste0("patch", x))
  rownames(m_distance) <- colnames(m_distance) <- sapply(seq_len(n_patch), function(x) paste0("patch", x))

  return(list(adjacency_matrix = m_adj,
              distance_matrix = m_distance,
              df_patch = dplyr::tibble(patch_id = 1:n_patch,
                                       branch_id = as.numeric(df_id$branch),
                                       environment = c(v_env),
                                       n_patch_upstream = c(v_wa))))
}
