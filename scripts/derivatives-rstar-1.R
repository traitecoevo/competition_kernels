# Fundamental result that we need: Species can coexist *at
# equilibrium* only if each species becomes limited by the resource
# for which it has, compared to its competitors, the highest
# requirement.

# Several species competing for one resource (p. 2684)
#
#   The species that has the lowest resource requirement for the
#   limiting resouce (i.e., the species with the lowest R^*_{1i}) will
#   displace all other species.
#
# So let's look at the behaviour of the system with a single resource
# and move towards understanding how the competition function might
# look like.
library(Revolve)
library(plyr)
library(numDeriv)

mat <- rstar_matrices(rstar_mat_1, rstar_mat_1)

m <- make_rstar(mat, S=1)
x <- matrix(0.5, nrow=2)
y0 <- 1
t <- seq(0, 30, length=201)

res <- m$run(x, y0, t)
matplot(res$t, cbind(res$R, res$y), type="l", lty=1:2, xlab="Time",
        ylab="Abundance (red), Resource (black)",
        ylim=c(0, max(res$R, res$y)))

eq <- m$single_equilibrium(x)
abline(h=unlist(eq), col=1:2, lty=3)

# Displace the solution from equilibrium population density and look
# at the approach to a new level of resources:
dy <- eq$y * 0.1
res1 <- m$run_fixed_density(x, eq$y - dy, t, eq$R)
res2 <- m$run_fixed_density(x, eq$y + dy, t, eq$R)
matplot(res1$t, cbind(res1$R, res2$R), type="l", col=2:3)
# Determined by running things out a bunch:
eq1 <- m$equilibrium_R(x, eq$y - dy)
eq2 <- m$equilibrium_R(x, eq$y + dy)
abline(h=c(eq$R, eq1$R, eq2$R), col=1:3, lty=3)
# Determined analytically:
points(max(res1$t), m$single_equilibrium_R(x, eq$y - dy), col=2)
points(max(res1$t), m$single_equilibrium_R(x, eq$y + dy), col=3)

# Next, we start working towards the instantaneous growth rate of a
# new type at this equilibrium

# Look at the fitness landscape: how does the instantaneous growth
# rate look with respect to mutant C and K.
xx <- seq(0, 1, length=101)
x.K <- rbind(xx, x[2], deparse.level=0) # K varying
x.C <- rbind(x[2], xx, deparse.level=0) # C varying

# Fitness does not vary with respect to c when rare, because species
# only consume resources when not rare.  Which doesn't mean that there
# will be no competition between species with different c values but
# the same K values of course: a different c value can drive a
# different species extinct!
plot(xx, m$fitness(x.C, x, eq$y, eq$R), type="l")
abline(h=0, col="grey", lty=3)

# Fitness does vary with K.
w.mutant <- m$fitness(x.K, x, eq$y, eq$R)
plot(xx, w.mutant, type="l")
abline(h=0, col="grey", lty=3)
abline(v=x[1], lty=2)

# But then there's a broader view of what competition is here: The
# species with K above the resident species are all strongly competed
# against by the resident species, regardless of what the derivatives
# show -- because they can't grow there.  That's the insight that we'd
# get from the subtraction case:
w.empty <- m$fitness(x.K, x, 0, m$parameters$get()[["S"]])
plot(xx, w.empty, type="l")
abline(h=0, col="grey", lty=3)
abline(v=x[1], lty=2)

plot(xx, w.mutant - w.empty, type="l")
abline(h=0, col="grey", lty=3)
abline(v=x[1], lty=2)

# Next, see what happpens to this fitness as we vary the population
# density.  This involves recomputing the new resource density for the
# species density:
z <- jacobian(function(y) m$fitness(x.K, x, y), eq$y)

plot(xx, z, type="l")
abline(h=0, col="grey", lty=3)
abline(v=x[1], lty=2)

op <- par(mfrow=c(3, 1), mar=c(2.5, 4, .5, .5))
plot(xx, w.mutant - w.empty, type="l")
abline(h=0, col="grey", lty=3)
abline(v=x[1], lty=2)

plot(xx, z, type="l")
abline(h=0, col="grey", lty=3)
abline(v=x[1], lty=2)

plot(xx, z / w.empty, type="l")
abline(v=x[1], lty=2)
par(op)

# To make this all a little clearer, this is what the fitness gradient
# is measuring.  Focus on a single mutant, whose trait is at the
# resident species' R*:
x.mutant <- rbind(m$Rstar(x), x[[2]])

# Then consider a range of densities of the resident species.
y.resident <- eq$y + seq(-1, 1, length=101)

# Fitness of the mutant as a function of these resident densities:
w.mutant <- sapply(y.resident, function(y) m$fitness(x.mutant, x, y))

w.mutant0 <- m$fitness(x.mutant, x, eq$y)

plot(y.resident, w.mutant, type="l")
abline(h=w.mutant0, v=eq$y, lty=3, col="grey")
z.mutant <- jacobian(function(y) m$fitness(x.mutant, x, y), eq$y)
points(eq$y, w.mutant0, col="red")
abcline(eq$y, w.mutant0, z.mutant, lty=2, col="red")
