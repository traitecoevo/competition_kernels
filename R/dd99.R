dd99_parameters <- function(r=1.0, K_0=500.0, s2_C=0.16, s2_K=1) {
  list(r=r, K_0=K_0, s2_C=s2_C, s2_K=s2_K)
}

dd99_fitness <- function(x_invade, x_res, N_res, parameters) {
  r <- parameters$r
  a <- dd99_competition(x_invade, x_res, parameters)
  K <- dd99_carrying_capacity(x_invade, parameters)
  r * (1 - colSums(N_res * a) / K)
}

dd99_carrying_capacity <- function(x, parameters) {
  K_0 <- parameters$K_0
  s2_K <- parameters$s2_K
  K_0 * exp(- x^2 / (2 * s2_K))
}

dd99_max_growth_rate <- function(x, parameters) {
  rep(parameters$r, length.out=length(x))
}

dd99_competition <- function(x_invade, x_res, parameters) {
  s2_C <- parameters$s2_C
  exp(- outer(x_res, x_invade, "-")^2 / (2 * s2_C))
}
