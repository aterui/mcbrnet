#' Return adjacency in a linear habitat
#'
#' @param n Number of nodes
#' @param start_id Start node ID
#'
#' @export

fun_adj <- function(n, start_id = 1) {

  if (n == 1) m_y <- cbind(1, NA)

  if (n == 2) m_y <- cbind(c(1, 2), c(2, 1))

  if (n > 2) {

    y1 <- c(1, rep(2:(n - 1), each = 2), n)
    y2 <- c(2, sapply(2:(n - 1), function(i) c(i - 1, i + 1)), n - 1)
    m_y <- cbind(y1, y2)

  }

  colnames(m_y) <- c("patch_1", "patch_2")

  m_y <- m_y + (start_id - 1)

  return(m_y)
}

#' To dispersal matrix
#'
#' @inheritParams mcsim
#'
#' @export

fun_disp_mat <- function(n_patch,
                         landscape_size,
                         theta,
                         xy_coord = NULL,
                         distance_matrix = NULL,
                         dispersal_matrix = NULL) {

  # with no landscape information
  if (is.null(xy_coord) &&
      is.null(distance_matrix) &&
      is.null(dispersal_matrix)) {

    message("neither xy_coord nor distance_matrix is given: generate a square landscape with landscape_size (default: 10) on a side")

    v_x_coord <- runif(n = n_patch,
                       min = 0,
                       max = landscape_size)

    v_y_coord <- runif(n = n_patch,
                       min = 0,
                       max = landscape_size)

    df_xy_coord <- dplyr::tibble(x_coord = v_x_coord,
                                 y_coord = v_y_coord)

    m_distance <- data.matrix(dist(df_xy_coord,
                                   diag = TRUE,
                                   upper = TRUE))

    m_dispersal <- data.matrix(exp(-theta * m_distance))

  }

  # with xy coordinates
  if (!is.null(xy_coord) &&
      is.null(distance_matrix) &&
      is.null(dispersal_matrix)) {

    if (nrow(xy_coord) != n_patch) stop("row numbers must match n_patch")
    if (ncol(xy_coord) != 2) stop("the number of columns must be two, describing x- and y-cooridnates")

    colnames(xy_coord) <- c("x_coord", "y_coord")

    df_xy_coord <- dplyr::as_tibble(xy_coord)

    m_distance <- data.matrix(dist(df_xy_coord,
                                   diag = TRUE,
                                   upper = TRUE))

    m_dispersal <- data.matrix(exp(-theta * m_distance))

  }

  # with distance matrix
  if (!is.null(distance_matrix) &&
      is.null(dispersal_matrix)) {

    message("distance_matrix is provided: distance_matrix is used to simulate dispersal process")
    if (!is.matrix(distance_matrix)) stop("distance matrix must be provided as matrix")
    if (nrow(distance_matrix) != n_patch) stop("invalid dimension: distance matrix must have a dimension of n_patch * n_patch")
    if (any(diag(distance_matrix) != 0)) stop("invalid distance matrix: diagonal elements must be zero")

    df_xy_coord <- NULL

    m_distance <- distance_matrix

    m_dispersal <- data.matrix(exp(-theta * m_distance))
  }

  # with dispersal matrix
  if (!is.null(dispersal_matrix)) {

    message("dispersal_matrix is provided: dispersal_matrix is used to simulate dispersal process")
    if (!is.matrix(dispersal_matrix)) stop("dispersal_matrix must be provided as matrix")
    if (nrow(dispersal_matrix) != n_patch) stop("invalid dimension: dispersal_matrix must have a dimension of n_patch * n_patch")
    if (any(diag(dispersal_matrix) != 0)) stop("invalid dispersal_matrix: diagonal elements must be zero")

    df_xy_coord <- NULL

    m_distance <- NULL

    m_dispersal <- data.matrix(dispersal_matrix)

  }

  diag(m_dispersal) <- 0

  return(list(m_distance = m_distance,
              m_dispersal = m_dispersal,
              df_xy_coord = df_xy_coord))

}

#' Dispersal
#'
#' @param x Population size matrix
#' @param v_p_dispersal Vector of dispersal probability
#' @param m_dispersal Dispersal matrix
#'
#' @export

fun_dispersal <- function(x,
                          v_p_dispersal,
                          m_dispersal) {

  if (any(diag(m_dispersal) != 0)) stop("error in m_dispersal")

  # x is n_species (row) x n_patch (column) matrix
  # m_e_hat: expected number of emigrants from each habitat patch
  m_e_hat <- x * v_p_dispersal

  # m_e_sum: summed across patches
  v_e_sum <- rowSums(m_e_hat)

  # immigration potential for each patch = m_e_hat x m_dispersal (unit: individuals
  m_i_raw <- m_e_hat %*% m_dispersal

  # v_i_sum: summed across patches
  v_i_sum <- rowSums(m_i_raw)

  # insert 1 if v_i_sum = 0 (to avoid NaN for division)
  v_i_sum[v_i_sum == 0] <- 1

  # immigration probability = patch-specific potential / summed potential across patches
  m_i_prob <- m_i_raw / v_i_sum

  # expected immigrants: immigration prob. x total emigrants across patches
  m_i_hat <- m_i_prob * v_e_sum
  m_n_prime <- x + m_i_hat - m_e_hat

  return(m_n_prime)
}

#' Get fcl
#'
#' @param x Community matrix
#' @param delta Preference to basal over ig-prey
#'
#' @export

fun_get_fcl <- function(x, delta) {

  if (ncol(x) != length(delta)) stop("invalid dimension: error in x or delta")

  # raw fcl - 1, 2, 3
  v_fcl_raw <- apply(X = x,
                     MARGIN = 2,
                     FUN = function(y) {
                       ifelse(y[1] > 0,
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

#' Get n_branch
#'
#' @inheritParams brnet
#'
#' @export

fun_get_n_branch <- function(n_patch,
                             p_branch) {

  if (p_branch > 0 && p_branch < 1) {

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

  return(n_branch)

}

#' IGP function
#'
#' @param x Community matrix for basal, intraguild prey (ig-prey), and intraguild predator (ig-predator)
#' @param r_b Maximum reproductive rate of basal species
#' @param e Energetic conversion efficiency; must be given by the order of basal to ig-prey, basal to ig-predator, and ig-prey to ig-predator
#' @param k Carrying capacity of basal species
#' @param a Attack rate. Must be given by the order of basal to ig_prey, basal to ig-predator, and ig-prey to ig-predator
#' @param h Handling time. Must be given by the order of basal to ig_prey, basal to ig-predator, and ig-prey to ig-predator
#' @param s Strength of prey switching from ig-prey to basal
#'
#' @export

fun_igp <- function(x,
                    r_b = 5,
                    k = 100,
                    e = c(0.8, 0.8, 0.8),
                    a = c(2, 2, 2),
                    h = c(0.1, 0.1, 0.1),
                    s = 0
) {

  # check inputs ------------------------------------------------------------

  if (!is.matrix(x)) stop("error in x; x must be matrix")
  if (dim(x)[1] != 3) stop("error in x; dim(x)[1] must be 3")
  if (length(s) != 1) stop("error in s; s must be a scalar")
  if (s < 0 || s > 1) stop("error in s; s must be 0 - 1")

  if (any(x < 0) ||
      any(r_b < 0) ||
      any(k < 0) ||
      any(a < 0) ||
      any(e < 0) ||
      any(h < 0)
  ) stop("negative values detected in parameters or community matrix")

  # define parameters -------------------------------------------------------

  ## basal species, intraguild prey, intraguild predator
  v_n_b <- x[1, ]
  v_n_c <- x[2, ]
  v_n_p <- x[3, ]

  ## growth rate / conversion efficiency
  e_bc <- e[1]
  e_bp <- e[2]
  e_cp <- e[3]

  ## attack rate
  a_bc <- a[1]
  a_bp <- a[2]
  a_cp <- a[3]

  ## attack rate x handling time
  h_bc <- a[1] * h[1]
  h_bp <- a[2] * h[2]


  # trophic dynamics: B to C ------------------------------------------------

  ## Basal species; Beverton-Holt growth
  nu <- (r_b - 1) / k
  v_n_b_plus <- v_n_b * (r_b / (1 + nu * v_n_b))

  ## v_p_bc: fraction of prey survived after predation by consumer
  v_p_bc <- exp(-((a_bc * v_n_c) / (1 + h_bc * v_n_b_plus)))

  ## predation C on B & conversion
  v_n_b_minus_c <- v_n_b_plus * v_p_bc
  v_n_c_plus <- e_bc * v_n_b_plus * (1 - v_p_bc)


  # trophic dynamics: B to P, C to P ----------------------------------------

  ## preference function P on B & C
  ## phi: switching function
  ## s: strength of switching
  phi <- s * (v_n_b_minus_c - v_n_c_plus) / (v_n_b_minus_c + v_n_c_plus)
  phi <- ifelse(is.nan(phi), 0, phi) # NaN produced when n_b = n_c = 0
  om_b <- 1 + phi # preference to basal over ig-prey
  om_c <- 1 - phi # preference to ig-prey over basal

  ## v_p_xp: fraction of prey survived after predation by predator
  v_p_bp <- exp(-((om_b * a_bp * v_n_p) / (1 + h_bp * v_n_b_minus_c)))
  v_p_cp <- exp(-((om_c * a_cp * v_n_p) / (1 + h_bp * v_n_c_plus)))

  ## Basal
  ### predation P on B
  v_n_b_hat <- v_n_b_minus_c * v_p_bp

  ## IG prey
  ### conversion x basal n x faction captured x predation P on C
  v_n_c_hat <- v_n_c_plus * v_p_cp

  ## IG predator
  ### conversion x (basal n x fraction remained) x fraction captured
  ### conversion x consumer n x fraction captured
  v_n_p_hat <- e_bp * v_n_b_minus_c * (1 - v_p_bp) + e_cp * v_n_c_plus * (1 - v_p_cp)


  # omnivory level ----------------------------------------------------------

  v_delta <- (e_bp * v_n_b_minus_c * (1 - v_p_bp)) / v_n_p_hat
  v_delta[is.nan(v_delta)] <- 0


  # export ------------------------------------------------------------------

  m_n_hat <- rbind(v_n_b_hat,
                   v_n_c_hat,
                   v_n_p_hat)

  dimnames(m_n_hat) <- NULL

  return(list(m_n_hat = m_n_hat,
              delta = v_delta))

}

#' Interaction matrix
#'
#' @inheritParams mcsim
#'
#' @export

fun_int_mat <- function(n_species,
                        alpha,
                        min_alpha,
                        max_alpha,
                        interaction_type) {

  if (interaction_type == "constant") {

    if (alpha < 0 || length(alpha) != 1) stop("invalid value of alpha - the value must be a positive scalar")

    m_interaction <- matrix(alpha,
                            nrow = n_species,
                            ncol = n_species)

  } else {

    if (interaction_type != "random") stop("invalid interaction_type")
    if (is.null(min_alpha) || is.null(max_alpha)) stop("provide min_alpha and max_alpha")
    if (min_alpha < 0 || max_alpha < 0) stop("invalid values of min_alpha and/or max_alpha - values must be positive values")
    if (min_alpha > max_alpha) stop("max_alpha must exceed min_alpha")

    alpha <- runif(n = n_species * n_species,
                   min = min_alpha,
                   max = max_alpha)

    m_interaction <- matrix(alpha,
                            nrow = n_species,
                            ncol = n_species)

  }

  diag(m_interaction) <- 1

  return(m_interaction)

}

#' Adjacency matrix in branching networks
#'
#' @inheritParams brnet
#' @param n_branch Number of branches
#'
#' @export

fun_m_adj <- function(n_patch,
                      p_branch,
                      n_branch,
                      n_patch_free) {

  # adjacency matrix: linear network ----------------------------------------

  if (p_branch == 0 || n_branch == 1) {

    v_n_patch_branch <- n_patch
    m_adj <- matrix(0, nrow = n_patch, ncol = n_patch)
    x <- seq_len(n_patch - 1)
    y <- 2:n_patch
    m_adj[cbind(x, y)] <- 1
    m_adj[cbind(y, x)] <- 1

  }


  # adjacency matrix: branched network --------------------------------------

  if (p_branch > 0 && n_branch > 1) {

    # vector of the number of patches in each branch
    if (n_patch_free == FALSE) {
      repeat {

        repeat {

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
                                function(x) cbind(v_end_id[m_po[x, 1]], v_start_id[m_po[x, 2]]))

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
                                function(x) cbind(v_end_id[m_po[x, 1]], v_start_id[m_po[x, 2]]))

      m_confluence <- do.call(rbind, list_confluence)
      m_confluence <- rbind(m_confluence, m_confluence[, c(2, 1)])

    }

    m_neighbor_patch <- rbind(m_neighbor_inbranch, m_confluence)
    m_neighbor_patch <- m_neighbor_patch[complete.cases(m_neighbor_patch), ]

    m_adj <- matrix(0, nrow = n_patch, ncol = n_patch)
    m_adj[m_neighbor_patch] <- 1

  }

  return(list(v_n_patch_branch = v_n_patch_branch,
              m_adj = m_adj))
}

#' Weighted-mean patch attributes
#'
#' @param n_branch Number of branches
#' @param x Adjacency matrix to be converted
#' @param mean_source Mean value of source attribute
#' @param sd_source SD of source attribute
#' @param sd_lon SD of longitudinal change
#' @param m_distance Distance matrix
#' @param rho Longitudinal autocorrelation
#' @param v_wa Vector of watershed area
#' @param logit Logit transformation of mean source attribute
#'
#' @export

fun_patch_attr <- function(x,
                           n_branch,
                           mean_source,
                           sd_source,
                           sd_lon,
                           m_distance,
                           rho = 1,
                           v_wa,
                           logit = FALSE) {

  # extract upstream adjacency
  m_adj_up <- x
  m_adj_up[lower.tri(m_adj_up)] <- 0
  n_patch <- nrow(x)

  # upstream tributary proportion
  m_wa_prop <- t(apply(X = m_adj_up,
                       MARGIN = 1,
                       function(x) {
                         (x * v_wa) / ifelse(sum(x) == 0,
                                             yes = 1,
                                             no = sum(x * v_wa))
                       }
  )
  )

  n_source <- 0.5 * (n_branch + 1)
  source <- which(rowSums(m_adj_up) == 0)

  v_z_dummy <- v_z <- v_attr <- rep(0, n_patch)
  v_z[source] <- v_attr[source] <- rnorm(n = n_source,
                                         mean = ifelse(logit == TRUE,
                                                       yes = boot::logit(mean_source),
                                                       no = mean_source),
                                         sd = sd_source)
  v_z_dummy[source] <- 1

  if (!(rho <= 1 && rho >= 0)) stop("rho must be between 0 and 1")

  for (i in seq_len(max(m_distance[1, ]))) {

    v_eps <- rep(0, n_patch)
    v_eps[v_z_dummy != 0] <- rnorm(n = length(v_eps[v_z_dummy != 0]),
                                   mean = 0,
                                   sd = sd_lon)

    v_z <- m_wa_prop %*% ((rho * v_z) + v_eps)
    v_z_dummy <- m_wa_prop %*% v_z_dummy
    v_attr <- v_z + v_attr

  }

  return(c(v_attr))
}

#' Resample function
#'
#' @param x Vector
#' @param ... Additional arguments passed to \code{sample}
#'
#' @export

resample <- function(x, ...) x[sample.int(length(x), ...)]

#' To matrix
#'
#' @param x Scalar or vector
#' @param n_species Number of species
#' @param n_patch Number of patches
#' @param param_attr Indicate species or patch attribute
#'
#' @export

fun_to_m <- function(x,
                     n_species,
                     n_patch,
                     param_attr) {

  # check input -------------------------------------------------------------

  if (!is.vector(x)) stop("x must be a vector")


  # conversion to matrix ----------------------------------------------------

  ## patch-wise
  if (param_attr == "patch") {

    if (length(x) == 1) {

      message("single value is given for species or patch attributes: assume identical across species or patches")

      v_x <- rep(x = x, times = n_patch)
      m_x <- matrix(x,
                    nrow = n_species,
                    ncol = n_patch)

    } else {

      if (length(x) != n_patch) stop("x must have a length of one or n_patch")

      v_x <- x
      m_x <- matrix(rep(x = x,
                        each = n_species),
                    nrow = n_species,
                    ncol = n_patch)

    } # ifelse

  } # if (param_attr == "patch")

  ## species-wise
  if (param_attr == "species") {

    if (length(x) == 1) {

      message("single value is given for species or patch attributes: assume identical across species or patches")

      v_x <- rep(x = x, times = n_species)
      m_x <- matrix(x,
                    nrow = n_species,
                    ncol = n_patch)

    } else {

      if (length(x) != n_species) stop("x must have a length of one or n_species")

      v_x <- x
      m_x <- matrix(rep(x = x,
                        times = n_patch),
                    nrow = n_species,
                    ncol = n_patch)

    } # ifelse

  } # if (param_attr == "species")


  # export ------------------------------------------------------------------

  return(list(v_x = v_x,
              m_x = m_x))

}

#' To vector
#'
#' @param x Scalar or vector
#' @param n Replication (n_species or n_patch)
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

fun_to_v <- function(x, n) {

  if (!is.vector(x)) stop("x must be a vector")

  if (length(x) == 1) {

    message("single value is given for species or patch attributes: assume identical across species or patches")
    v_k <- rep(x = x,
               times = n)

  } else {

    if (length(x) != n) stop("x must have a length of one or n_species/n_patch")
    v_k <- x

  }

  return(v_k)

}

#' Watershed area
#'
#' @param x Adjacency matrix to be converted
#'
#' @export

fun_wa <- function(x) {

  m_adj_up <- x
  m_adj_up[lower.tri(m_adj_up)] <- 0
  n_patch <- nrow(x)

  m_wa <- matrix(0,
                 ncol = n_patch,
                 nrow = n_patch)

  m_identity <- diag(1,
                     nrow = n_patch,
                     ncol = n_patch)

  for (i in seq_len(n_patch)) {

    m_identity <- m_identity %*% m_adj_up
    m_wa[m_identity != 0 & m_wa == 0] <- 1

    if (length(which(m_wa[upper.tri(m_wa)] == 0)) == 0) break

  }

  diag(m_wa) <- 1

  return(m_wa)
}

#' Yield partial derivatives
#'
#' @param r Numeric. Intrinsic growth rate of species i.
#' @param a Numeric vector. Interaction coefficients.
#' @param i Integer. Species numeric ID.
#' @inheritParams stability
#'
#' @export

fun_partial <- function(r, a, i, x0, model) {

  # check input -------------------------------------------------------------
  if (length(r) != 1)
    stop(paste("Input 'r' must be a scalar:",
               "'r' has length", length(r)))

  if (length(i) != 1)
    stop(paste("'i' must be a scalar:",
               "'i' has length", length(i)))

  if (length(a) != length(x0))
    stop(paste("Invalid inputs in 'a' or 'x0':",
               "'a' has length =", length(a),
               "while 'x0' has length =", length(x0)))

  if (model %in% c("ricker", "bh", "glv"))
    stop(paste("Invalid model type: 'model' must be either of",
               "'ricker', 'bh', or 'glv'"))

  # model formula -----------------------------------------------------------
  ## declare vectorized parameters and variables
  v_a <- paste0("a[", seq_len(length(a)), "]")
  v_x <- paste0("x[", seq_len(length(x0)), "]")
  arg <- paste(c("r", "x", "a"), collapse = ", ")

  ## linear combination
  lcm <- paste(v_a, "*", v_x)

  ## function
  f <- NULL

  ## Generalized Lotka-Volterra model
  if (model == "glv") {
    ## get a model formula
    fm <- c("r", lcm)
    m <- paste0("x[", i, "]", " * ",
                "(",
                paste(fm, collapse = " + "),
                ")")

    ## function text for evaluation
    f_text <- parse(text = paste0("f <- function(", arg, ") {",
                                  m,
                                  "}"))

    eval(f_text)
  }

  ## Ricker model
  if (model == "ricker") {
    ## get a model formula
    fm <- c("r", lcm)
    m <- paste0("x[", i, "]", " * ",
                "exp(",
                paste(fm, collapse = " + "),
                ")")

    ## function text for evaluation
    f_text <- parse(text = paste0("f <- function(", arg, ") {",
                                  m,
                                  "}"))

    eval(f_text)
  }

  ## Beverton-Holt model
  if (model == "bh") {
    ## get a model formula
    m <- paste0("x[", i, "]", " * ", "exp(r)",
                " * ",
                "(1 + ",
                paste(lcm, collapse = " + "),
                ") ** -1")

    ## function text for evaluation
    f_text <- parse(text = paste0("f <- function(", arg, ") {",
                                  m,
                                  "}"))

    eval(f_text)
  }

  ## return partial derivatives evaluated at x0 (equilibrium)
  return(pracma::jacobian(f, x0 = x0, r = r, a = a))
}
