# look like.
library(Revolve)
library(plyr)
library(numDeriv)

m <- make_huisman_2001(S=c(1,1), matrices=huisman_matrices(huisman_mat_2))
x <- 0.3
y0 <- 1
t <- seq(0, 30, length=201)

res <- m$run(x, y0, t)
matplot(res$t, cbind(res$R, res$y), type="l", lty=c(2,2,1), xlab="Time",
        ylab="Abundance (solid), Resource (dashed)")

eq <- m$single_equilibrium(x)
abline(h=unlist(eq), col=1:3, lty=3)

# Displace the solution from equilibrium and look at the new level of
# resources:
dy <- eq$y * 0.1
res1 <- m$run_fixed_density(x, eq$y - dy, t, eq$R)
res2 <- m$run_fixed_density(x, eq$y + dy, t, eq$R)
eq1 <- m$equilibrium_R(x, eq$y - dy)
eq2 <- m$equilibrium_R(x, eq$y + dy)
matplot(res1$t, cbind(res1$R, res2$R), type="l",
        col=c(2,2,3,3))
abline(h=eq$R, lty=3, col=1)
abline(h=eq1$R, lty=3, col=2)
abline(h=eq2$R, lty=3, col=3)

# Next, we start working towards the instantaneous growth rate of a
# new type at this equilibrium
x2 <- seq(0, 1, length=101)

# Fitness in the empty environment and in the presence of the resident:
fitness.empty    <- m$fitness(x2, x, 0, m$parameters$get()[["S"]])
fitness.resident <- m$fitness(x2, x, eq$y)

plot(x2, fitness.resident, type="l")
abline(h=0, col="grey", lty=3)
abline(v=x[1], lty=2)

plot(x2, fitness.empty, type="l")
abline(h=0, col="grey", lty=3)
abline(v=x[1], lty=2)

plot(x2, fitness.resident - fitness.empty, type="l")
abline(h=0, col="grey", lty=3)
abline(v=x[1], lty=2)

# Next, see what happpens to this fitness as we vary the population
# density.  This involves recomputing the new resource density for the
# species density:
z <- jacobian(function(y) m$fitness(x2, x, y), eq$y)

plot(x2, z, type="l")
abline(h=0, col="grey", lty=3)
abline(v=x[1], lty=2)
abline(v=m$Rstar(x), col="blue")

plot(x2, z / fitness.empty, type="l")
abline(v=x[1], lty=2)
