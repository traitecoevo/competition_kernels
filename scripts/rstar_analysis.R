## 1: Why is the optimal strategy all the way over at 0.9
## 2: Fix the break point in simulations
## 3: Different measures in main paper or on SM?

source("rstar_new.R")
source("rstar_support.R")
source("rstar_plot.R")
source("support.R")

p <- rstar_parameters(rstar_mat_2_tradeoff,
                      matrix(c(.3, .7), nrow=2),
                      S=0.5)

col1 <- "blue"
col1_tr <- make_transparent(col1, .2)
col_died <- make_transparent("grey", .2)
x1 <- matrix(0.3, nrow=1)
N1 <- 1.0

rstar_plot(x1, p, col=col1)
for (i in 1:100) {
  rstar_trajectory(x1, N1, p, col=col1_tr, col_died=col_died, S=runif(2))
}

## Consider this case here: supply coefficient at [0.5, 0.5] and
## resident at [0.3] (x1 above).  We'll need the equilibrium resource
## and population size.
p$S <- c(0.5, 0.5)
eq <- rstar_single_equilibrium(x1, p)

## Then, we'll consider mutants changing the K coefficient.  The odd
## choice of length here means that the resident trait is included.
x_invade <- rbind(seq(0, 1, length.out=301))

## Compute carrying capacity, maximum growth rate and fitness at the
## equilibrium R for each of these points.
K_invade <- rstar_carrying_capacity(x_invade, p)
r_invade <- rstar_max_growth_rate(x_invade, p)
w_invade <- rstar_fitness_given_R(x_invade, eq$R, p)

plot(x_invade, w_invade, type="l", las=1,
     xlab="Species trait (K)", ylab="Invasion fitness")
abline(v=eq$x[1, ], h=0, lty=2)

plot(x_invade, K_invade, type="l")
plot(x_invade, r_invade, type="l")

## 1: Based on solving the LV equations for alpha:
alpha_LV <- competition_LV(w_invade, r_invade, K_invade, eq$N)
plot(x_invade, alpha_LV, type="l")
abline(h=1, v=eq$x[1], lty=2)

## 2: Competitive effect of the resident on the mutant.  This requires
## two parts.

## Four different gradients needed:
##   invader growth  / resident density
##   invader growth  / invader density
##   resident growth / resident density
##   resident growth / invader density

## For plain finite differencing (i.e. for things based from zero N)
## we'll use this displacement.
dN <- 1e-4

## 1: derivative of invader growth wrt resident density:
f_dgi_dNj <- function(N_res, x_invade, x_res, p) {
  rstar_fitness_given_N(x_invade, x_res, N_res, p)
}
dgi_dNj <- drop(numDeriv::jacobian(f_dgi_dNj, eq$N,
                                   x_invade=x_invade, x_res=x1, p=p))

## 2: derivative of invader growth wrt invader density
f_dgi_dNi <- function(N_invade, N_res, x_invade, x_res, p) {
  f <- function(x) {
    rstar_fitness_given_N(cbind(x),
                          cbind(x_res, x),
                          c(N_res, N_invade), p)
  }
  apply(x_invade, 2, f)
}
dgi_dNi <- (f_dgi_dNi(dN, eq$N, x_invade, eq$x, p) - w_invade) / dN

## 3: derivative of resident growth with respect to resident density:
f_dgj_dNj <- function(N_res, x_res, p) {
  rstar_fitness_given_N(x_res, x_res, N_res, p)
}
dgj_dNj <- numDeriv::grad(f_dgj_dNj, eq$N,
                          x_res=eq$x, p=p)

## 4: derivative of resident growth with respect to invader density:
f_dgj_dNi <- function(N_invade, N_res, x_invade, x_res, p) {
  f <- function(x) {
    rstar_fitness_given_N(x_res,
                          cbind(x_res, x),
                          c(N_res, N_invade), p)
  }
  apply(x_invade, 2, f)
}
dgj_dNi <- (f_dgj_dNi(dN, eq$N, x_invade, eq$x, p) - 0) / dN

## And moar:

## dgi / dNi @ Ni = Ki, Nj = 0
f1 <- function(x, p) {
  x <- cbind(x)
  eq <- rstar_single_equilibrium(x, p)
  ## TODO: These answers don;t line up, so I've got this a bit wrong.
  ## Should not be too much work to get them sorted though.
  numDeriv::grad(function(N) rstar_fitness_given_N(x, x, N, p), eq$N)
  ## f <- function(N) {
  ##   R <- rstar_single_equilibrium_R(x, N, p)
  ##   rstar_fitness_given_R(x, R, p)
  ## }
  ## numDeriv::grad(f, eq$N)
}
tmp1 <- apply(x_invade, 2, f1, p)

## dgi / dNj @ Ni = Ki, Nj = 0
f2 <- function(x, p, x_res, dN) {
  x <- cbind(x)
  eq <- rstar_single_equilibrium(x, p)
  (rstar_fitness_given_N(x, cbind(x, x_res), c(eq$N, dN), p) - 0.0) / dN
}
tmp2 <- apply(x_invade, 2, f2, p, x1, dN)

alpha_Xij <- tmp2 / tmp1
alpha_Yji <- dgi_dNj / dgi_dNi # this is X^*_ij in the paper
alpha_Zij <- dgj_dNi / dgj_dNj # this is something else entirely.

plot(x_invade, alpha_LV, type="l", main="LV")
abline(h=1.0, v=x1, lty=2)
plot(x_invade, alpha_Xij, type="l", main="X")
abline(h=1.0, v=x1, lty=2)
plot(x_invade, alpha_Yji, type="l", main="Y")
abline(h=1.0, v=x1, lty=2)
plot(x_invade, alpha_Zij, type="l", main="Z")
abline(h=1.0, v=x1, lty=2)

## Try with (0.5, 0.5): the winning strategy?

rstar_plot(x2, p, col=col1)
for (i in 1:100) {
  rstar_trajectory(x2, N1, p, col=col1_tr, col_died=col_died, S=runif(2))
}


p$S <- c(0.5, 0.5)
x2 <- matrix(0.9, nrow=1)
eq <- rstar_single_equilibrium(x2, p)

## Then, we'll consider mutants changing the K coefficient.  The odd
## choice of length here means that the resident trait is included.
x_invade <- rbind(seq(0, 1, length.out=301))

## Compute carrying capacity, maximum growth rate and fitness at the
## equilibrium R for each of these points.
K_invade <- rstar_carrying_capacity(x_invade, p)
r_invade <- rstar_max_growth_rate(x_invade, p)
w_invade <- rstar_fitness_given_R(x_invade, eq$R, p)

plot(x_invade, w_invade, type="l", las=1,
     xlab="Species trait (K)", ylab="Invasion fitness")
abline(v=eq$x[1, ], h=0, lty=2)

alpha_LV <- competition_LV(w_invade, r_invade, K_invade, eq$N)
plot(x_invade, alpha_LV, type="l")
abline(h=1, v=eq$x[1], lty=2)


## For plain finite differencing (i.e. for things based from zero N)
## we'll use this displacement.
dN <- 1e-4

## 1: derivative of invader growth wrt resident density:
f_dgi_dNj <- function(N_res, x_invade, x_res, p) {
  rstar_fitness_given_N(x_invade, x_res, N_res, p)
}
dgi_dNj <- drop(numDeriv::jacobian(f_dgi_dNj, eq$N,
                                   x_invade=x_invade, x_res=x2, p=p))

## 2: derivative of invader growth wrt invader density
f_dgi_dNi <- function(N_invade, N_res, x_invade, x_res, p) {
  f <- function(x) {
    rstar_fitness_given_N(cbind(x),
                          cbind(x_res, x),
                          c(N_res, N_invade), p)
  }
  apply(x_invade, 2, f)
}
dgi_dNi <- (f_dgi_dNi(dN, eq$N, x_invade, eq$x, p) - w_invade) / dN

## 3: derivative of resident growth with respect to resident density:
f_dgj_dNj <- function(N_res, x_res, p) {
  rstar_fitness_given_N(x_res, x_res, N_res, p)
}
dgj_dNj <- numDeriv::grad(f_dgj_dNj, eq$N,
                          x_res=eq$x, p=p)

## 4: derivative of resident growth with respect to invader density:
f_dgj_dNi <- function(N_invade, N_res, x_invade, x_res, p) {
  f <- function(x) {
    rstar_fitness_given_N(x_res,
                          cbind(x_res, x),
                          c(N_res, N_invade), p)
  }
  apply(x_invade, 2, f)
}
dgj_dNi <- (f_dgj_dNi(dN, eq$N, x_invade, eq$x, p) - 0) / dN

## And moar:

## dgi / dNi @ Ni = Ki, Nj = 0
f1 <- function(x, p) {
  x <- cbind(x)
  eq <- rstar_single_equilibrium(x, p)
  ## TODO: These answers don;t line up, so I've got this a bit wrong.
  ## Should not be too much work to get them sorted though.
  numDeriv::grad(function(N) rstar_fitness_given_N(x, x, N, p), eq$N)
  ## f <- function(N) {
  ##   R <- rstar_single_equilibrium_R(x, N, p)
  ##   rstar_fitness_given_R(x, R, p)
  ## }
  ## numDeriv::grad(f, eq$N)
}
tmp1 <- apply(x_invade, 2, f1, p)

## dgi / dNj @ Ni = Ki, Nj = 0
f2 <- function(x, p, x_res, dN) {
  x <- cbind(x)
  eq <- rstar_single_equilibrium(x, p)
  (rstar_fitness_given_N(x, cbind(x, x_res), c(eq$N, dN), p) - 0.0) / dN
}
tmp2 <- apply(x_invade, 2, f2, p, x2, dN)

alpha_Xij <- tmp2 / tmp1
alpha_Yji <- dgi_dNj / dgi_dNi # this is X^*_ij in the paper
alpha_Zij <- dgj_dNi / dgj_dNj # this is something else entirely.

plot(x_invade, alpha_LV, type="l", main="LV")
abline(h=1.0, v=x2, lty=2)
plot(x_invade, alpha_Xij, type="l", main="X")
abline(h=1.0, v=x2, lty=2)
plot(x_invade, alpha_Yji, type="l", main="Y")
abline(h=1.0, v=x2, lty=2)
plot(x_invade, alpha_Zij, type="l", main="Z")
abline(h=1.0, v=x2, lty=2)

plot(x_invade, alpha_LV, type="l")
lines(x_invade, alpha_Xij, col="darkgrey")
lines(x_invade, alpha_Yji, col="red")
abline(h=1.0, v=x2, lty=2)
