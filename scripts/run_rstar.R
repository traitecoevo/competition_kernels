source("rstar_new.R")
source("rstar_support.R")

## OK, I've driven the species extinct, which is not great.

p <- rstar_parameters(rstar_mat_2_tradeoff, rstar_mat_2_tradeoff, 1.0)

## Species traits:
x <- rbind(c(0.2, 0.7),
           c(0.2, 0.7))
## Initial population densities:
N0 <- c(0.3, 0.7)

tt <- seq(0, 300, length.out=301)

p$C(x)
p$K(x)
rstar_min_p(x, p$S, p) # 0.05555556 0.05882353
rstar_p(x, p$S, p)

rstar_dNdt(x, N0, p$S, p) #  0.09166667  0.23676471
rstar_dRdt(x, N0, p$S, p) # -0.3215686  -0.2568627

obj_eq <- rstar_equilibrium(x, N0, p)
obj_tr <- rstar_run(tt, x, N0, p)

op <- par(mfrow=c(2, 1), mar=c(4.1, 4.1, .5, .5))
matplot(obj_tr$t, obj_tr$N, type="l", xlab="", ylab="Abundance")
points(rep(max(tt), 2), obj_eq$N, col=1:2, cex=c(1, .5))
matplot(obj_tr$t, obj_tr$R, type="l", xlab="Time", ylab="Resource")
points(rep(max(tt), 2), obj_eq$R, col=1:2, cex=c(1, .5))
par(op)

## Now, look at the equilibrium of the single species setup.  For this
## we'll just start the second species at a density of 0 and with an
## arbitrary 0.5 for their state.
N0[2] <- 0
x[, 2] <- 0.5
x[, 1] <- 0.4

obj_eq <- rstar_equilibrium(x, N0, p)
obj_tr <- rstar_run(tt, x, N0, p)

op <- par(mfrow=c(2, 1), mar=c(4.1, 4.1, .5, .5))
matplot(obj_tr$t, obj_tr$N, type="l", xlab="", ylab="Abundance",
        ylim=c(0, 2))
points(rep(max(tt), 2), obj_eq$N, col=1:2)
matplot(obj_tr$t, obj_tr$R, type="l", xlab="Time", ylab="Resource",
        ylim=c(0, 2))
points(rep(max(tt), 2), obj_eq$R, col=1:2)
par(op)
