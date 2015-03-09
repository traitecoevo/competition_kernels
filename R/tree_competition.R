tree_competition <- function(x_resident, obj, n=101L) {
  x_resident <- check_point(trait_matrix(x_resident, obj$trait),
                            obj$bounds)
  x_mutant <- seq_log_range(obj$bounds, n)
  trait <- rownames(obj$bounds)
  N_resident <- obj$K(x_resident)

  message("*** Setting up resident population")
  p <- obj$p0
  p$strategies <- strategy_list(x_resident, p)
  p$seed_rain <- N_resident
  p <- build_schedule(p)

  message("*** Computing invasion fitness")
  obj$w <- fitness_landscape(trait_matrix(x_mutant, trait), p)

  r <- obj$r(x_mutant)
  K <- obj$K(x_mutant)

  ## Then competition:
  obj$x_mutant   <- x_mutant
  obj$x_resident <- x_resident
  obj$N_resident <- N_resident
  obj$alpha <- compute_alpha(obj$w, r, K, N_resident)
  obj
}

tree_competition_prepare <- function(trait, bounds=NULL, p0=NULL,
                                     check_viable=TRUE,
                                     n=20L, parallel=FALSE) {
  if (is.null(p0)) {
    p0 <- tree_base_parameters()
  }
  if (check_viable || is.null(bounds)) {
    if (is.null(bounds)) {
      bounds <- bounds_infinite(trait)
    }
    bounds <- viable_fitness(bounds, p0)
  }
  bounds <- check_bounds(bounds)

  x <- seq_log_range(bounds, n)
  message("*** Computing r")
  r <- max_growth_rate(trait_matrix(x, trait), p0)

  message("*** Computing K")
  K <- carrying_capacity(trait_matrix(x, trait), p0, parallel=parallel)

  list(p0=p0,
       trait=trait,
       bounds=bounds,
       r=splinefun_log(x, r),
       K=splinefun_loglog(x, K))
}

tree_base_parameters <- function() {
  p <- ebt_base_parameters()
  p$control$equilibrium_nsteps <- 30
  p$control$equilibrium_solver_name <- "hybrid"
  p$disturbance_mean_interval <- 30.0

  ## Will be overriden in specific models:
  p$hmat <- 15.0
  p$lma  <- 0.2
  
  p
}
