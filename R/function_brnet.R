#' Generate a random branching network
#'
#' @param n_patch number of patches in a network.
#' @param p_branch branching probability (success probability of a geometric distribution).
#' @param min_env minimum value of environmental condition in source streams (minimum of a uniform distribution).
#' @param max_env maximum value of environmental condition in source streams (maximum of a uniform distribution).
#' @param rho strength of spatial autocorrelation in environmental condition. The environmental condition at patch i \eqn{x}\out{<sub>i</sub>} is determined as \eqn{x}\out{<sub>i</sub>}\eqn{ = \rho}x\out{<sub>i-1</sub>}\eqn{ + \epsilon}\out{<sub>i</sub>}, where \eqn{\epsilon}\out{<sub>i</sub>} is the random variable drawn from a normal distribution with mean 0 and SD \eqn{\sigma}\out{<sub>env</sub>}.
#' @param sd_env SD of environmental noise (\eqn{\sigma}\out{<sub>env</sub>}).
#' @param randomize_patch logical indicating whether randomize or not patches. If \code{FALSE}, the function may generate a biased network with ordered patches. Default \code{TRUE}.
#' @param plot logical. If \code{FALSE}, a plot of the generated network will not be shown. Default \code{TRUE}.
#'
#' @return \code{adjacency_matrix} adjacency matrix for the generated network.
#' @return \code{distance_matrix} distance matrix for the generated network.
#' @return \code{environment} environmental values for patches.
#' @return \code{watershed_area} number of patches upstream. Patch 1 is set to be the root patch (i.e., outlet).
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @examples
#' brnet(n_patch = 100, p_branch = 0.5)
#'
#' @export
#'
brnet <- function(n_patch = 100,
                  p_branch = 0.5,
                  min_env = -1,
                  max_env = 1,
                  rho = 1,
                  sd_env = 0.1,
                  randomize_patch = TRUE,
                  plot = TRUE){

  # define functions and variables ------------------------------------------

  resample <- function(x, ...) x[sample.int(length(x), ...)]

  fun_adj <- function(n, start_id = 1){
    if(n == 1){ Y <- cbind(1, NA) }
    if(n == 2){
      Y <- cbind(c(1, 2), c(2, 1))
    }
    if(n > 2){
      y1 <- c(1, rep(2:(n-1), each = 2), n)
      y2 <- c(2, sapply(2:(n-1), function(i) c(i-1, i+1)), n-1)
      Y <- cbind(y1, y2)
    }
    colnames(Y) <- c("patch_1", "patch_2")
    Y <- Y + (start_id - 1)
    return(Y)
  }# fun_adj


  if(p_branch > 0 & p_branch < 1){
    repeat{
      n_branch <- rbinom(n = 1, size = n_patch, prob = p_branch)
      if(n_branch%%2 == 1) break
    }# repeat
  }else{
    if(p_branch == 0){ n_branch <- 1 }
    if(p_branch == 1){
      if(n_patch%%2 == 0) stop("n_patch must be an odd number when p_branch = 1")
      n_branch <- n_patch
    }
  }# ifelse


  # adjacency matrix: linear network ----------------------------------------

  if(p_branch == 0|n_branch == 1){
    Adj <- matrix(0, nrow = n_patch, ncol = n_patch)
    x <- 1:(n_patch-1)
    y <- 2:n_patch
    Adj[cbind(x,y)] <- 1
    Adj[cbind(y,x)] <- 1
  }# if(p_branch == 0|n_branch == 1)


  # adjacency matrix: branched network --------------------------------------

  if(p_branch > 0 & n_branch > 1){

    # vector of the number of patches in each branch
    repeat{
      repeat{
        v_n_patch_branch <- rgeom(n = n_branch, prob = p_branch) + 1
        if(sum(v_n_patch_branch) >= n_patch) break
      }# repeat
      if(sum(v_n_patch_branch) == n_patch) break
    }# repeat

    # start_id, end_id, and neighbor list for each branch
    v_end_id <- cumsum(v_n_patch_branch)
    v_start_id <- v_end_id - (v_n_patch_branch - 1)
    list_neighbor_branch <- lapply(1:n_branch, function(i) fun_adj(n = v_n_patch_branch[i], start_id = v_start_id[i]) )
    Neighbor_branch <- do.call(rbind, list_neighbor_branch)

    # combine parent and offspring branches at confluences
    if(n_branch == 3){
      parent <- c(1, 1)
      offspg <- c(2, 3)
      Po <- cbind(parent, offspg)
      list_confluence <- lapply(1:nrow(Po), function(x){ cbind(v_end_id[Po[x,1]], v_start_id[Po[x,2]]) })
      Confluence <- do.call(rbind, list_confluence)
      Confluence <- rbind(Confluence, Confluence[,c(2,1)])
    }else{
      n_confluence <- 0.5*(n_branch - 1)
      v_parent_branch <- 1:n_confluence
      v_offspg_branch <- 2:n_branch

      Offspg <- NULL
      for(i in n_confluence:1){
        v_y <- resample(v_offspg_branch[v_offspg_branch > v_parent_branch[i]], size = 2)
        v_offspg_branch <- setdiff(v_offspg_branch, v_y)
        Offspg <- cbind(v_y, Offspg)
      }

      parent <- rep(v_parent_branch, each = 2)
      offspg <- c(Offspg)
      Po <- cbind(parent, offspg)

      list_confluence <- lapply(1:nrow(Po), function(x){ cbind(v_end_id[Po[x,1]], v_start_id[Po[x,2]]) })
      Confluence <- do.call(rbind, list_confluence)
      Confluence <- rbind(Confluence, Confluence[,c(2,1)])
    }# ifelse

    Neighbor <- rbind(Neighbor_branch, Confluence)
    Neighbor <- Neighbor[complete.cases(Neighbor),]

    Adj <- matrix(0, nrow = n_patch, ncol = n_patch)
    Adj[Neighbor] <- 1
  }# if(p_branch > 0)


  # distance matrix ---------------------------------------------------------

  D <- matrix(0, ncol = n_patch, nrow = n_patch)
  A <- diag(x = 1, nrow = n_patch, ncol = n_patch)
  for(i in 1:n_patch){
    A <- A%*%Adj
    D[A!=0 & D==0] <- i
    if(length(which(D==0)) == 0) break
  }# for(i in 1:n_patch)
  diag(D) <- 0
  A <- NULL


  # upstream watershed area ------------------------------------------------

  Adj_up <- Adj
  Adj_up[lower.tri(Adj_up)] <- 0

  Wm <- matrix(0, ncol = n_patch, nrow = n_patch)
  A <- diag(1, nrow = n_patch, ncol = n_patch)

  for(i in 1:n_patch){
    A <- A%*%Adj_up
    Wm[A!=0 & Wm==0] <- 1
    if(length(which(Wm[upper.tri(Wm)]==0)) == 0) break
  }# for(i in 1:n_patch)
  diag(Wm) <- 1

  v_wa <- rowSums(Wm)


  # environmental condition -------------------------------------------------

  Pm <- t(apply(X = Adj_up, MARGIN = 1, function(x) x*v_wa/ifelse(sum(x)==0, 1, sum(x*v_wa) ) ) )

  n_source <- 0.5*(n_branch + 1)
  source <- which(rowSums(Adj_up)==0)
  v_z_dummy <- v_z <- v_env <- rep(0, n_patch)
  v_z[source] <- v_env[source] <- runif(n_source, min = min_env, max = max_env)
  v_z_dummy[source] <- 1

  if(!(rho <= 1&rho >= 0)) stop("rho must be between 0 and 1")
  for(i in 1:max(D[1,]) ){
    v_eps <- rep(0, n_patch)
    v_eps[v_z_dummy != 0] <- rnorm(n = length(v_eps[v_z_dummy != 0]), mean = 0, sd = sd_env)
    v_z <- Pm%*%((rho*v_z) + v_eps)
    v_z_dummy <- Pm%*%v_z_dummy
    v_env <- v_z + v_env
  }# for(i in 1:max(D[1,]) )


  # randomize nodes ---------------------------------------------------------

  branch <- unlist(lapply(1:n_branch, function(x) rep(x, each = v_n_patch_branch[x])))
  patch <- 1:n_patch

  if(randomize_patch == TRUE){
    if(n_branch > 1){
      df_id <- data.frame(branch = as.character(c(1, resample(2:n_branch)))) %>%
        dplyr::left_join(data.frame(patch, branch = as.character(branch)), by = "branch")
      v_wa <- v_wa[df_id$patch]
      v_env <- v_env[df_id$patch]
      Adj <- Adj[df_id$patch, df_id$patch]
      D <- D[df_id$patch, df_id$patch]
    }else{
      df_id <- data.frame(branch = 1, patch = 1:n_patch)
    }# ifelse(n_branch > 1)
  }else{
    df_id <- data.frame(branch = branch, patch = patch)
  }# ifelse(randomize_patch == TRUE)

  # visualization -----------------------------------------------------------

  if(plot == T){
    adj <- igraph::graph.adjacency(adjmatrix = Adj, mode = "undirected")
    colvalue <- data.frame(color = viridis::viridis(n_patch, alpha = 0.6), value = sort(v_env) )
    layout.tree <- igraph::layout_as_tree(adj, root = 1, flip.y = F)
    par(mar = c(5,8,4,3))
    igraph::plot.igraph(adj, layout = layout.tree,
                        vertex.size = I(scale(v_wa, center = min(v_wa), scale = max(v_wa)-min(v_wa)) + 0.3)*10,
                        vertex.label = NA,
                        vertex.frame.color = grey(0.5),
                        vertex.color = colvalue$color[match(v_env, colvalue$value)])
    plotfunctions::gradientLegend(valRange = range(v_env), color = viridis::viridis(n_patch),
                                  pos = 0.8, side = 2, dec = 2)
    pc <- c(plotfunctions::getCoords(0, side = 1), plotfunctions::getCoords(1, side = 2))
    text(x = pc[1], y = pc[2], labels = "Environmental value", adj = 1)
  }# if(plot == T)

  return(list(adjacency_matrix = Adj,
              distance_matrix = D,
              environment = c(v_env),
              watershed_area = c(v_wa),
              branch_id = as.numeric(df_id$branch)) )
}
