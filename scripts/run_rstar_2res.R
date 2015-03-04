## Fundamental result that we need: Species can coexist *at
## equilibrium* only if each species becomes limited by the resource
## for which it has, compared to its competitors, the highest
## requirement.
source("rstar_new.R")
source("rstar_support.R")
source("rstar_plot.R")

col_died <- make_transparent("grey", .2)
p <- rstar_parameters(rstar_mat_2_tradeoff, rstar_mat_2_tradeoff, S=1.0)

## x <- matrix(0.5, nrow=2)
## N0 <- 1


## mat <- rstar_matrices(rstar_mat_2_tradeoff, rstar_mat_2_tradeoff)
## m <- rstar(mat, 1)


## One species:
col1 <- "blue"
col1_tr <- make_transparent(col1, .2)
x1 <- matrix(0.2, nrow=2)
N1 <- 1.0

## END OF WEDNESDAY, 2015-03-04: up to here---
rstar_plot(x1, p, col=col1)
for (i in 1:100) {
  rstar_trajectory(x1, N1, p, col=col1_tr, col_died=col_died, S=runif(2))
}

## Two species:
x2 <- cbind(x1, 0.7, deparse.level=0)
N2 <- c(1, 1)
col2 <- c(col1, "red")
col2_tr <- make_transparent(col2, .2)

rstar_plot(x2, p, col=col2)
for (i in 1:100) {
  rstar_trajectory(x2, N2, p, col=col2_tr, col_died=col_died, S=runif(2))
}

## Now, displace the single species system from equilibrium and look
## at the new level of resources:
p$S <- c(1, 1)
eq <- rstar_single_equilibrium(x1, p)
rstar_plot(x1, p, col=col2)
rstar_trajectory(x1, N1, p, S=p$S)

x1 <- x1
x2 <- cbind(x1, .3, deparse.level=0)

R1 <- rstar_single_equilibrium(x2[, 1, drop=FALSE], p)$R # .8166, .2666
R2 <- rstar_single_equilibrium(x2[, 2, drop=FALSE], p)$R # .6714, .2333
testthat::expect_that(rstar_single_equilibrium(x2, p)$R,
                      testthat::equals(cbind(R1, R2, deparse.level=0)))

xx <- matrix(rep(seq(0, 1, by=.1), each=2), 2)
rstar_single_equilibrium(xx, p)

dN <- eq$N * 0.1
N1 <- eq$N - dN
N2 <- eq$N + dN
t <- seq(0, 30, length.out=201)

res1 <- rstar_run_R(t, x1, N1, p)
res2 <- rstar_run_R(t, x1, N2, p)

eq_y1 <- rstar_equilibrium_R(x1, N1, p)
eq_y2 <- rstar_equilibrium_R(x1, N2, p)

matplot(t, cbind(res1$R, res2$R), type="l",
        col=c(1,1,2,2), lty=c(1,2,1,2))
abline(h=eq_y1$R, lty=3, col=1)
abline(h=eq_y2$R, lty=3, col=2)

# Solving this analytically is going to be harder than for the one
# resource case.  However, it's possible that we can actually skip the
# hard work and solve this for each resource separately and then see
# which ones are plausible solutions (i.e., *assume* that the first or
# the second resource is limiting and then see if that makes sense).
# Realistically I'm going to need an pen and paper for this though.
#
# The alternative approach is to solve numerically.  That should be
# pretty easy to do, either with nleqslv or with runsteady, or
# possibly a uniroot approach on the equation itself.
#
# Even if a semianalytic solution is found, we're going to need the
# same trick before to work out where it lands using the supply curve.

# Look at the fitness landscape: how does the instantaneous growth
# rate look with respect to K (& C):
x_invade <- rbind(seq(0, 1, length.out=301), eq$x[2])
w_invade <- rstar_fitness_given_R(x_invade, eq$R, p)

plot(x_invade[1,], w_invade, type="l", las=1)
abline(h=0, col="grey", lty=3)
abline(v=eq$x, lty=2)
