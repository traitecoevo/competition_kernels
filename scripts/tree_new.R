tree_competition <- function(x_resident, obj, n=101L) {
  x_resident <- check_point(x_resident, obj$bounds)
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
                                     n=20L, parallel=FALSE) {
  if (is.null(p0)) {
    p0 <- tree_base_parameters()
  }
  if (is.null(bounds)) {
    bounds <- viable_fitness(bounds_infinite(trait), p0)
  }
  bounds <- check_bounds(bounds)

  x <- seq_log_range(bounds, n)
  message("*** Computing r")
  r <- max_growth_rate(trait_matrix(x, trait), p0)

  message("*** Computing K")
  K <- carrying_capacity(trait_matrix(x, trait), p0, parallel=parallel)

  list(p0=p0,
       bounds=bounds,
       r=splinefun_log(x, r),
       K=splinefun_loglog(x, K))
}


## Compute competition coefficient given invader fitness, max growth
## rate, carrying capacity and resident density.
compute_alpha <- function(w_i, r_i, K_i, N_r) {
  (1 - 1 / r_i * w_i) * K_i / N_r
}

make_invasion_fitness <- function(trait_matrix, seed_rain, p) {
  p$strategies <- strategy_list(trait_matrix, p)
  p$seed_rain <- seed_rain
  p <- build_schedule(p)
  function(x) {

  }
}

tree_base_parameters <- function() {
  p <- ebt_base_parameters()
  p$control$equilibrium_nsteps <- 30
  p$control$equilibrium_solver_name <- "hybrid"
  ## I don't think that we actually used these before, but they *are*
  ## used in the successional_diversity models, and I've set up the
  ## python interface to work this way.  See tree #127
  ##   p$strategy_default$c_r1 <- 0.5
  ##   p$strategy_default$c_r2 <- 0
  p$disturbance_mean_interval <- 30.0
  p
}
