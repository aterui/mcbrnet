#' Simulate trophic dynamics in 3 species communities through space and time
#'
#' @param n_patch single numeric value. Total number of patches in simulation
#' @param n_0 `NULL` or numeric vector length either 1 or 3. Starting abundances for species in each patch. If only one value is supplied it is assumed to be the same for all three species. Can designate starting values for each species with a vector of length 3. If `NULL` starting abundances for B, C, and P are 0.8, 0.5, and 0.25 * `k_base`, respectively.
#' @param dist_mat square distance matrix describing distance of all patches. If it is `NULL` simulates a square landscape with `length = landscape_size`, with `n_patch` randomly distributed.
#' @param landscape_size single numeric value. Length of a square landscape on a side.
#' @param adjacency_matrix square matrix describing which patches are adjacent to one another. For branching river networks, use `net$adjacency_matrix` from `mcbrnet::brnet()` as the input argument.
#' @param k_function character string, one of: `c("patches-upstream", "environment", NULL)`. `"patches-upstream"` calculates carrying capacity based on the number of nodes upstream. See `?k_n_upstream` for details.
#' `"environment"` calculates K based on `df_patch$environment` value returned from `brnet()` function. `NULL` assumes that `k = k_base` in all patches.
#' @param n_upstream numeric vector of `length = n_patch`. Used if `k_function = "patches-upstream"`. Use `df_patch$n_patch_upstream` from `mcbrnet::brnet()`.
#' @param environment_value numeric vector of `length = n_patch`. used if `k_function = "environment"`. Use `df_patch$environment` from `mcbrnet::brnet()`.
#' @param k_base single numeric value. Used to calculate carrying capacity, K. see `?k_function_internal` for more.
#' @param k_c single numeric value. Constant used to calculate K when `k_function = "patches-upstream"`. Default = 10
#' @param k_min_exponent single numeric value. Minimum value of a uniform distribution for the exponent used to calculate K when `k_function = "patches_upstream"`
#' @param k_max_exponent single numeric value. Maximum value of a uniform distribution for the exponent used to calculate K when `k_function = "patches_upstream"`
#' @param p_dispersal numeric value describing the probability of dispersal. For example, a value of 0.1 means that 10% of the population abundance emmigrates from each patch in each time step. If length = 1, assumed to be same for all three species. To set probability independently for species, B, C, and P, use a vector of length = 3.
#' @param theta numeric value of length 1 or 3. Dispersal parameter describing dispersal capability of species. If length(theta) = 1 it is the same for all three species. Can be set for each species with vector of length = 3.
#' @param r_max single numeric value. Maximum reproductive number for the basal species, B, of the Beverton-Holt model.
#' @param ebc single Numeric value. Conversion efficiency of turning B biomass into new C biomass.
#' @param ebp single Numeric value. Conversion efficiency of turning B biomass into new P biomass
#' @param ecp single Numeric value. Conversion efficiency of turning C biomass into new P biomass
#' @param alphabc single numeric value. Parameter controlling the predation of resource B by consumer C
#' @param betabc single numeric value. Parameter controlling the predation of resource B by consumer C
#' @param alphap single numeric value. Parameter controlling the predation of both resources by consumer P
#' @param betap single numeric value. Parameter controlling the predation of both resources by consumer P
#' @param P_pref single numeric value between 0-1 or NULL. Predator preference of B over C. If `P_pref = NULL`, default, calculates predator preference independently for each patch and time step, using the `prey_preference()` function. Preference is a function of conversion efficiency rates `(ebp, ecp)` and observed prey abundances of both prey species in each patch and time step. i.e., Holling's "type 3" response.
#' If value is supplied to `P_pref`, it is fixed in all patches and all time steps. i.e., `P_pref = 0.25` means that the predator will always allocate 25% of its search effort to B, and 75% to C.
#' @param s0 numeric value length should be either 1 or 3. Base survival probability for all (`length = 1`) or each (`length = 3`) species
#' @param disturb_type character string, one of: `c("regional", or NULL)`.
#'  `NULL` does not apply any disturbance events throughout the simulation.
#' `"regional"` applies a disturbance to the entire network in a given time step with probability = `disturb_p`. The magnitude of the disturbance varies based on the `disturb_value`.
#' Future releases are intended to include a `point-source` argument, where disturbances occur randomly in a patch, and the magnitude decays downstream. This is not currently implemented.
#' @param disturb_value Numeric vector of length = `n_patch`.
#' This is optimized to work with branching river networks generated from the the `mcbrnet` package by taking the output of `df_patch$disturbance` from `brnet()` as an input. See `?brnet` for more details.
#' For 2D habitats, the input value is ignored. The disturbance value is sampled randomly by converting the `disturb_mag_mean` and `disturb_mag_sd` arguments using the inverse logit function. Currently, the value for each patch is randomly sampled and set for the entire simulation run. i.e., there is no spatial correlation or variation through time. This has not been optimized for 2D habitats
#' @param disturb_p single numeric value between 0-1. The probability that a disturbance occurs at a given time step.
#' @param disturb_mag_mean single numeric value between 0-1. This controls the mean disturbance magnitude in 2D habitat network structures. This argument is ignored in simulations with branching network structure.
#' @param disturb_mag_sd single numeric value. This controls the variation in the disturbance magnitude in 2D habitat network structures. This argument is ignored in simulations with branching network structure.
#' @param n_burnin Single numeric value. The number of time-steps to occur before recording values. default = 200,
#' @param n_timestep Single numeric value. The number of time-steps to be saved. default =  1000
#' @param plot_patch_dynamics logical. If `TRUE` plots population abundances through time for 5 random patches.
#' @param plot_mean_fcl
#' @param plot_fcl_state logical. If `TRUE` plots the proportion of patches with a given food chain length through all time-steps. FCL is a state value describing the number and identity of species present: FCL = 0 when no species are present, 1 = B only, 2 = B + C, 2.5 = B + P, 3 = B + C + P.
#' @return `sp_dynamics` a data frame containing simulated IGP community dynamics.
#' @return `df_patch` a data frame containing information describing each patch
#' @return `fcl_state_prop` a data frame containing the proportion of patches that have a given community composition for each time step.
#' @return `mean_fcl` a data frame containing the mean food chain length across all patches in each time step.
#' @return `sim_params` a data frame which contains the parameter values used in the simulation. This data frame is a single row, unless one of the following argument is of length = 3: `p_dispersal`, `s0`, or `theta`, in which case the data frame will be 3 rows, corresponding to species B, C, and P, respectively.
#'
#' @importFrom dplyr %>% bind_rows filter pull select summarize group_by
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes facet_grid facet_wrap geom_line geom_point label_both labeller labs geom_step theme_bw geom_smooth geom_vline
#' @importFrom boot logit inv.logit
#'
#' @export
#'
#' @examples
#' igp_sim()
#'
igp_sim <- function(n_patch = 20,
                    n_0 = NULL,
                    dist_mat = NULL,
                    landscape_size = 10,
                    adjacency_matrix = NULL,
                    k_function = NULL,
                    #c("patches-upstream", "environment")
                    n_upstream = NULL,
                    environment_value = NULL,
                    k_base = 150,
                    k_c = 10,
                    k_min_exponent = 1.10,
                    k_max_exponent = 1.35,
                    p_dispersal = 0.1,
                    theta = 1,
                    r_max = 2.5,
                    ebc = 2,
                    ebp = 2,
                    ecp = 2,
                    alphabc = 4,
                    alphap = 4,
                    betabc = 20,
                    betap = 20,
                    P_pref = NULL, # preference of B over C
                    s0 = 0.75,
                    disturb_type = NULL, #c(NULL, "point-source", "regional")
                    disturb_value = NULL,
                    disturb_p = 1e-4,
                    disturb_mag_mean = 0.9,
                    disturb_mag_sd = 0.1,
                    n_burnin = 200,
                    n_timestep = 1000,
                    plot_patch_dynamics = FALSE,
                    plot_mean_fcl = FALSE,
                    plot_fcl_state = FALSE){

  param_df <- data.frame(
    k_base = k_base,
    k_function = ifelse(is.null(k_function), NA, k_function),
    r_max = r_max,
    alphabc = alphabc, betabc = betabc, ebc = ebc,
    alphap = alphap, betap = betap, ebp = ebp, ecp =ecp,
    P_pref = ifelse(is.null(P_pref), NA, P_pref),
    p_dispersal = p_dispersal, theta = theta,
    s0 = s0,
    disturb_type = ifelse(is.null(disturb_type), NA, disturb_type),
    disturb_p = disturb_p, disturb_mag_mean = disturb_mag_mean,
    disturb_mag_sd = disturb_mag_sd)
  if(nrow(param_df) == 1){
    row.names(param_df) <- "all_sp"
  }
  if(nrow(param_df) == 3){
    row.names(param_df) <- c("B", "C", "P")
  }

  n_sp = 3
  k_input = k_base

  ## distance matrix ####
  dist_structure <- dist_mat_internal(dist_mat = dist_mat,
                    landscape_size = landscape_size,
                    n_patch = n_patch)
  dist_mat = dist_structure$dist_mat
  river_network_structure = dist_structure$river_network_structure

  if(river_network_structure == TRUE & is.null(disturb_value)){
    stop("Disturbance value from 'mcbrnet' is needed for branching habitats")
  }

  if(river_network_structure == FALSE & !is.null(disturb_value)){
    warning("Disturbance value is provided for 2D habitat, disturbance values overwritten with random sample")
  }

  # disturbance value for 2D habitats
  if(river_network_structure == FALSE){
    environment_value <- rnorm(n_patch, mean = 0, sd = 1)
    disturb_value <- boot::inv.logit(rnorm(n_patch,
                           mean = boot::logit(disturb_mag_mean),
                           sd = disturb_mag_sd))
  }

  if(length(theta) == 1){
    v_theta = rep(theta, 3)
  }
  if(length(theta) == 3){
    v_theta = theta
  }
  # probability of dispersal ####
  if(length(p_dispersal) == 1){
    v_p_dispersal <- rep(x = p_dispersal, times = n_sp)
    message(
      paste("1 value of p_dispersal supplied:",
            p_dispersal,
            "for all three species"))
  }else{
    if(length(p_dispersal)!= n_sp)
      stop("length of p_dispersal should be 1  or n_sp")
    v_p_dispersal <- p_dispersal
    message(
      paste("3 values of p_dispersal supplied: B_dispersal =",
            p_dispersal[1],
            "C_dispersal =", p_dispersal[2],
            "and P_dispersal =", p_dispersal[3]))
  }

  # survival probability s0####
  if (length(s0) == 1){
    v_s0 <- rep(x = s0, times = n_sp)
    message(
      paste("1 value of s0 supplied:",
            s0,
            "for all three species"))
  }else{
    if(length(s0)!= n_sp)
      stop("length of s0 should be 1  or n_sp = 3")
    v_s0 <- s0
    message(
      paste("3 values of s0 supplied: B_s0 =", s0[1],
            "C_s0 =", s0[2],
            "and P_s0 =", s0[3]))
  }

  # dispersal distance ####
  # distance decay of dispersal, theta
  if(length(theta) != 1 & length(theta) !=3){
    stop("in igp_sim(), Theta must be length 1 or 3")
  }
  if(length(theta) == 1){
    message(
      paste("1 value of theta supplied: B_theta = C_theta = P_theta =",
            theta))
  }
  if(length(theta) == 3){
    message(
      paste("3 values of theta supplied: B_theta =", v_theta[1],
            "C_theta =", v_theta[2],
            "and P_theta =", v_theta[3]))
  }

  # dispersal matrices by species ####
  # dispersal matrices per species
  # species B
  m_b_dispersal <- data.matrix(exp(-v_theta[1] * dist_mat))
  diag(m_b_dispersal) <- 0
  # species C
  m_c_dispersal <- data.matrix(exp(-v_theta[2] * dist_mat))
  diag(m_c_dispersal) <- 0
  # species P
  m_p_dispersal <- data.matrix(exp(-v_theta[3] * dist_mat))
  diag(m_p_dispersal) <- 0

  ## carrying capacity wrapper function ####
  b_k_list <- k_function_internal(
    k_function = k_function,
    k_c = k_c,
    k_min_exponent = k_min_exponent,
    k_max_exponent = k_max_exponent,
    k_base = k_base,
    r_max = r_max,
    n_upstream = n_upstream,
    n_patch = n_patch,
    river_network_structure = river_network_structure,
    environment_value = environment_value)
  k = b_k_list$k
  b = b_k_list$b

  # P_preference ####
  # predator preference fixed or variable
  # NULL or value from 0-1
  if(is.null(P_pref)){ # is.null(P_pref)
    fixed_P_pref = FALSE
    message("Predator preference varies with resource abundance")
  }
  if(is.numeric(P_pref) & length(P_pref == 1)){
    fixed_P_pref = TRUE
    P_pref = rep(P_pref, n_patch)
    message(paste("Predator preference is set at", P_pref))
  }

  ## initial community ####
  if(is.null(n_0)){
    N <- matrix(rpois(n = n_sp * n_patch,
                      lambda = c(k_input*0.8, k_input*0.5, k_input*0.25)),
                nrow = n_sp,
                ncol = n_patch)
    message(
      paste(
        "initial community abundance for B, C and P = 0.8, 0.5 and 0.25 *",
        k_base,
        "respectively"))
  } else{
    if(length(n_0) != 1 & length(n_0) != 3)
      stop("argument n_0 should be length 1, 3, or NULL")
    if(length(n_0) == 1){
      N <- matrix(rpois(n = n_sp*n_patch,
                        lambda = c(n_0, n_0, n_0)),
                  nrow = n_sp, ncol = n_patch)
      message(
        paste(
          "initial community abundance in all patches for all 3 species is",
          n_0))
    }
    if(length(n_0) == 3){
      N <- matrix(rpois(n = n_sp*n_patch,
                        lambda = c(n_0[1], n_0[2], n_0[3])),
                  nrow = n_sp, ncol = n_patch)
      message(paste("initial community abundance in all
                  patches for B =", n_0[1],
                  ", C =", n_0[2],
                  "and P =", n_0[3]),
              sep = "")
    }}

  # if basal species is 0 in a patch
  # make abundance of C and P = 0
  N[,N[1,] == 0] <- 0

  # result output ####
  output <- list()
  fcl_list <- list()


  counter = 1 # for recording output at end of simulation loop

  n_sim = n_burnin + n_timestep
  ## simulation ####
  for (i in seq_len(n_sim)){
    # dispersal at beginning of time step
    m_n_prime <- dispersal_n(N = N,
                             v_p_dispersal = v_p_dispersal,
                             v_theta = v_theta,
                             dist_mat = dist_mat,
                             m_b_dispersal = m_b_dispersal,
                             m_c_dispersal = m_c_dispersal,
                             m_p_dispersal = m_p_dispersal)

    # turn species abundances into integers
    N <- matrix(rpois(n = n_sp * n_patch,
                      lambda = m_n_prime),
                nrow = n_sp,
                ncol = n_patch)

    # pop_sim() ####
    pop_sim_out = pop_sim(N = N,
            P_pref = P_pref,
            fixed_P_pref = fixed_P_pref,
            alphabc = alphabc,
            betabc = betabc,
            ebc = ebc,
            alphap = alphap,
            betap = betap,
            ebp = ebp,
            ecp = ecp,
            v_s0 = v_s0,
            b = b,
            k = k,
            r_max = r_max)

    N = pop_sim_out$N
    lambda_b = pop_sim_out$lambda_b

    # disturbance ####
    disturb_result = disturb_internal(
      N = N,
      disturb_value = disturb_value,
      disturb_type = disturb_type,
      disturb_p = disturb_p,
      river_network_structure = river_network_structure)

    N = disturb_result$N
    patch_extinction = disturb_result$patch_extinction

    # end of cycle ####
    # number of individuals  in patch x at time t
    N <-  matrix(rpois(n = n_sp * n_patch,
                       lambda = N),
                 nrow = n_sp,
                 ncol = n_patch)

    # if basal species is extinct in a patch
    # make abundance of C and P = 0
    N[,N[1,] == 0] <- 0

    if(i > n_burnin){

      # food chain length state
      fcl = get_fcl(N = N,
                    lambda_b = lambda_b)

      # path_dynamics output
      out = cbind(1:n_patch,
                  t(N),
                  counter,
                  k,
                  patch_extinction,
                  fcl,
                  lambda_b)
      colnames(out) <- c("patch",
                         "B",
                         "C",
                         "P",
                         "time",
                         "basal_k",
                         "disturbance",
                         "fcl",
                         "lambda_b")
      output[[counter]] <- as.data.frame(out)
      fcl_list[[counter]] <- fcl_prop(get_fcl_state(N))
      counter = counter + 1
      #print(paste("loop iteration ", i))
    }
  }

  # make sure that these functions are properly imported
  dat <- dplyr::bind_rows(output)
  dat <- tidyr::pivot_longer(dat, 2:4, names_to = "species", values_to = "pop_density")

  mean_fcl_dat <- dat %>%
    group_by(time) %>%
    summarize(mean_fcl = mean(fcl))

  fcl_df <- data.frame(dplyr::bind_rows(fcl_list),
                       time = 1:n_timestep)
  fcl_df <- tidyr::pivot_longer(fcl_df,
                                1:5,
                                names_to = "community",
                                values_to = "proportion")

  # Plots -------------------------------------------------------------------
# plot FCL ####
  # plot of food chain length
  if(plot_fcl_state == TRUE){
    fcl_df$community = factor(
      fcl_df$community,
      levels = c("No_spp", "B_only", "B_C", "B_P", "B_C_P"))

    fcl_state_plot <-
      ggplot(fcl_df,
             aes(y = proportion, x = time, color = community)) +
      geom_step() +
      facet_wrap(.~community) +
      theme_bw() +
      labs(y = "proportion of patches",
           title = "Community composition") +
      NULL
    print(fcl_state_plot)
  }

  if(plot_mean_fcl == TRUE){
    message("mean FCL plot using 'loess' smoothing function")
    dist_dat <- dplyr::filter(dat, disturbance == 1)


    if(nrow(dist_dat) >0){
      mean_fcl_plot <- mean_fcl_dat %>%
        ggplot(aes(x = time, y = mean_fcl)) +
        geom_point() +
        geom_smooth(method = "loess") +
        geom_vline(#inherit.aes = FALSE,
                   data = dist_dat,
                   mapping = aes(xintercept = time),
                   linetype = 2,
                   size = 1,
                   alpha = 0.25) +
        theme_bw() +
        labs(y = "FCL",
             title = "Mean FCL through time") +
        NULL
    } else {
      mean_fcl_plot <- mean_fcl_dat %>%
        ggplot(aes(x = time, y = mean_fcl)) +
        geom_point() +
        geom_smooth(method = "loess") +
        theme_bw() +
        labs(y = "FCL",
             title = "Mean FCL through time") +
        NULL}
    print(mean_fcl_plot)
  }


  # plot patch dynamics ####
  # plot of patch dynamics
  if(plot_patch_dynamics == TRUE){

    plot_dat <- dat %>% filter(patch %in% sample(1:n_patch, 6))
    dist_dat <- dplyr::filter(plot_dat, disturbance == 1)

    if(nrow(dist_dat) >0){
      patch_plot <- ggplot(plot_dat,
                           aes(y = pop_density,
                               x = time,
                               color = species))+
        geom_line() +
        geom_point(inherit.aes = FALSE,
                   data = dist_dat,
                   mapping = aes(x = time,
                                 y = -5),
                   shape = 2,
                   size = 2) +
        facet_wrap(.~patch, labeller = label_both) +
        labs(y = "Abundance",
             title = "Dynamics for 6 random patches") +
        theme_bw() +
        NULL
      print(patch_plot)
    } else {
      patch_plot <- ggplot(plot_dat,
                           aes(y = pop_density,
                               x = time,
                               color = species))+
        geom_line() +
        facet_wrap(.~patch, labeller = label_both) +
        labs(y = "Abundance",
             title = "Dynamics for 6 random patches") +
        theme_bw() +
        NULL
      print(patch_plot)
    }
  }


  # return ####
  return(list(sp_dynamics = dat,
              df_patch = dplyr::tibble(patch_id = seq_len(n_patch),
                                       disturb_value = c(disturb_value),
                                       basal_k = c(k),
                                       dens_dep_reg = c(b)),
              fcl_state_prop = fcl_df,
              mean_fcl = mean_fcl_dat,
              sim_params = param_df))
}
