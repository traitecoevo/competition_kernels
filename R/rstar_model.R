## TODO: I'm handling the 'x' parameters in the opposite direction to
## plant; should be each column represents traits and each row
## represents species.
rstar_competition <- function(x_resident, p, N_resident=NULL, n=301L) {
  x_invade <- rbind(seq(0, 1, length.out=n))

  if (is.null(N_resident)) {
    eq <- rstar_single_equilibrium(x_resident, p)
  } else {
    eq <- rstar_equilibrium_R(x_resident, N_resident, p)
    eq$N <- N_resident
  }
  K_invade <- rstar_carrying_capacity(x_invade, p)
  r_invade <- rstar_max_growth_rate(x_invade, p)
  w_invade <- rstar_fitness_given_R(x_invade, eq$R, p)
  
  list(p=p,
       x_resident=x_resident,
       x_invade=x_invade,
       N_resident=eq$N,
       R_eq=eq$R,
       r_invade=r_invade,
       K_invade=K_invade,
       w_invade=w_invade,
       alpha=compute_alpha(w_invade, r_invade, K_invade, eq$N))
}

rstar_parameters <- function(K, C, S, r=1, m=0.25, D=0.25) {
  mat <- rstar_matrices(K, C)
  k <- mat$k
  ## Very basic checking:
  assert_scalar(r)
  assert_scalar(m)
  assert_scalar(D)
  if (length(S) == 1L) {
    S <- rep(S, k)
  } else if (length(S) != k) {
    stop("Invalid length S: must be 1 or ", k)
  }

  mat$S <- S
  mat$r <- r
  mat$m <- m
  mat$D <- D

  mat
}

rstar_initial_conditions <- function(N, parameters) {
  c(parameters$S, N)
}

## Specific growth rate at a particular resource level R
rstar_p <- function(x, R, parameters) {
  parameters$r * R / (parameters$K(x) + R)
}

## Resource level that would be required for particular growth rate g
rstar_pinv <- function(x, g, parameters) {
  g * parameters$K(x) / (parameters$r - g)
}

## Compute equilibrium level of each resource (i.e. level it will
## reach when limiting.  This is the resource level where the growth
## rate equals the mortality rate.
rstar_Rstar <- function(x, parameters) {
  rstar_pinv(x, parameters$m, parameters)
}

## For fitness we want the *per-capita growth rate*; that is (min_p -
## m), rather than y * (min_p - m).  We can't just divide by y because
## we want it when y = 0 (rather than just being close).  For now, I
## think we can just fake it by repeating a bunch of the stuff.
##
## The issue that we have in terms of interface is that this doesn't
## really depend on {x,y} as such, but on 'R', which is itself
## determined by {x,y}.  So the key determinant of fitness needs to be
## computed elsewhere for now.
rstar_fitness_given_R <- function(x_invade, R, parameters) {
  rstar_min_p(x_invade, drop(R), parameters) - parameters$m
}

rstar_fitness_given_N <- function(x_invade, x_res, N_res, parameters) {
  R <- rstar_equilibrium_R(x_res, N_res, parameters)$R
  rstar_fitness_given_R(x_invade, R, parameters)
}

rstar_max_growth_rate <- function(x, parameters) {
  ## rstar_fitness(x, parameters$S, parameters)
  rstar_min_p(x, parameters$S, parameters) - parameters$m
}

rstar_min_p <- function(x, R, parameters) {
  colMins(rstar_p(x, R, parameters))
}

## Compute the carrying capacity -- the population size that a species
## with trait x will reach when alone
rstar_carrying_capacity <- function(x, parameters) {
  Rstar <- rstar_Rstar(x, parameters)
  D <- parameters$D
  S <- parameters$S
  C <- parameters$C(x)
  m <- parameters$m

  colMins(D * (S - Rstar) / (m * C))
}

## Change in resources with respect to time
rstar_dRdt <- function(x, N, R, parameters) {
  D <- parameters$D
  S <- parameters$S
  C <- parameters$C(x)

  D * (S - R) - drop(C %*% (N * rstar_min_p(x, R, parameters)))
}

rstar_dNdt <- function(x, N, R, parameters) {
  m <- parameters$m

  N * (rstar_min_p(x, R, parameters) - m)
}

## To use derivs, you'll need to put the 'x' values into the
## parameters object.
rstar_derivs <- function(t, ode_y, parameters) {
  x <- parameters$x
  R <- ode_y[parameters$i_R]
  N <- ode_y[parameters$i_N]
  c(rstar_dRdt(x, N, R, parameters),
    rstar_dNdt(x, N, R, parameters))
}

rstar_run <- function(times, x, N, parameters) {
  parameters$x <- x
  ode_y <- rstar_initial_conditions(N, parameters)
  res <- .lsoda(ode_y, times, rstar_derivs, parameters)[, -1, drop=FALSE]
  list(t=times, x=x, R=res[, parameters$i_R], N=res[, parameters$i_N])
}

## Same but for R only, holding N constant.  Put N into the parameters
## object too.
rstar_derivs_R <- function(t, ode_y, parameters) {
  x <- parameters$x
  R <- ode_y
  N <- parameters$N
  rstar_dRdt(x, N, R, parameters)
}

rstar_run_R <- function(times, x, N, parameters) {
  parameters$x <- x
  parameters$N <- N
  ode_y <- parameters$S
  res <- .lsoda(ode_y, times, rstar_derivs_R, parameters)[, -1, drop=FALSE]
  list(t=times, x=x, R=res)
}

## General, numerical, equilibria:

## Here, N is the vector of *initial* population size.  It's needed.
rstar_equilibrium <- function(x, N, parameters, ...) {
  ode_y <- rstar_initial_conditions(N, parameters)
  parameters$x <- x
  res <- equilibrium(rstar_derivs, ode_y, parameters, ...)
  list(x=x, R=res[parameters$i_R], N=res[parameters$i_N])
}

## ## Given traits 'x' and fixed species densities 'N', compute the
## ## equilibrium.  Return the equilibrium resource availability at the
## ## same time
## rstar_equilibrium_R <- function(x, N, parameters, ...) {
##   ode_y <- parameters$S
##   parameters$x <- x
##   equilibrium(rstar_derivs_R, ode_y, parameters, ...)
## }

rstar_equilibrium_R <- function(x, N, parameters, ...) {
  ode_y <- parameters$S
  parameters$x <- x
  parameters$N <- N
  res <- equilibrium(rstar_derivs_R, ode_y, parameters, ...)
  list(x=x, R=res[parameters$i_R])
}

## Special equilibria:

## Given species traits `x`, compute equilibrium system.  If multiple
## species trait values are given, the output will be multiple
## separate equilibruia, rather than a joint community.
rstar_single_equilibrium <- function(x, parameters) {
  S <- parameters$S

  N <- rstar_carrying_capacity(x, parameters)
  R <- rstar_Rstar(x, parameters)

  is_extinct <- N < 0
  N[is_extinct] <- 0

  ## For 1 resource both R and N must be at equilibrium, but for >1
  ## resource, the nonlimiting resource need adjusting.
  if (nrow(R) > 1) {
    C <- parameters$C(x)
    len_min <- colMins((S - R) / C)
    R <- S - C * rep(len_min, each=nrow(R))
  }
  ## Fixup cases where species went extinct:
  R[, is_extinct] <- S

  list(x=x, N=N, R=R)
}

## Solve for the equilibrium level of a resource, given a vector of
## densities (and species traits).  Do do this, we solve dRdt == 0:
##     $D (S - R) - c N p(R) == 0$
##     $D * (S - R) - c N r R / (K + R) == 0$
##  As a quadratic with respect to 'R':
##     $-D * R^2 + (D (S - K) - c r N) * R + D K S == 0$
rstar_single_equilibrium_R <- function(x, N, parameters) {
  # x, N, parameters
  K <- parameters$K(x)
  C <- parameters$C(x)
  S <- parameters$S
  D <- parameters$D
  r <- parameters$r
  ans <- quadratic_roots(-D, (D * (S - K) - C * r * N), D * K * S)
  ans <- ans[ans >= 0 & ans <= S]
  if (length(ans) != 1) {
    stop("Did not find a unique solution")
  }
  ans
}

## Matrices.  This ugliness is required to create functions that turn
## some number of parameters into the K and C matrices.
rstar_matrices <- function(K, C) {
  if (is.matrix(K)) {
    rstar_matrices(make_rstar_mat_constant(K), C)
  } else if (is.matrix(C)) {
    rstar_matrices(K, make_rstar_mat_constant(C))
  } else {
    n_K <- attr(K, "npar", TRUE)
    n_C <- attr(C, "npar", TRUE)
    i_K <- seq_len(n_K)
    i_C <- seq_len(n_C) + n_K

    ## Number of resouces, computed with dummy vector
    k <- nrow(K(matrix(rep(0.5, n_K), n_K, 1)))
    if (nrow(C(matrix(rep(0.5, n_C), n_C, 1))) != k) {
      stop("K and C disagree on k")
    }

    ## Indices to the derivatives vector...
    i_R <- seq_len(k)
    i_N <- -i_R

    list(K=function(x) K(x[i_K, , drop=FALSE]),
         C=function(x) C(x[i_C, , drop=FALSE]),
         k=k, i_K=i_K, i_C=i_C, i_R=i_R, i_N=i_N)
  }
}

## Note that this does not use make_rstar_mat_constant because these
## matrices are constant with respect to both the number of resources
## *and* the number of species.  This exists only for checking against
## the paper.
rstar_matrices_fixed <- function(K, C) {
  if (!is.matrix(K) || !is.matrix(C)) {
    stop("K and C must both be matrices")
  }
  if (!identical(dim(K), dim(C))) {
    stop("K and C must have the same dimensions")
  }
  k <- nrow(K)
  matrices <- list(K=function(x) K, C=function(x) C,
                   n=ncol(K), k=k,
                   i_K=integer(0), i_C=integer(0),
                   i_R=seq_len(k))
}

## Make a constant function out of a column matrix or vector.
make_rstar_mat_constant <- function(M) {
  if (!is.matrix(M) || ncol(M) != 1) {
    stop("M must be a single column matrix")
  }
  ret <- function(x) {
    matrix(M, nrow=nrow(M), ncol=ncol(x))
  }
  attr(ret, "npar") <- 0
  ret
}

make_rstar_identity <- function(k) {
  ret <- function(x) {
    if (!(is.matrix(x) && nrow(x) == k)) {
      stop("x must be a matrix with ", k, " rows")
    }
    x
  }
  attr(ret, "npar") <- k
  ret
}

rstar_mat_1 <- make_rstar_identity(1L)
rstar_mat_2 <- make_rstar_identity(2L)
rstar_mat_2_tradeoff <- function(x) {
  rbind(x, 1 - x, deparse.level=0)
}
attr(rstar_mat_2_tradeoff, "npar") <- 1
