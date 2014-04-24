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
library(numDeriv)

mat <- rstar_matrices(rstar_mat_1, rstar_mat_1)

m <- make_rstar(mat, S=1)
sys0 <- sys(matrix(0.5, nrow=2), 1)
t <- seq(0, 30, length=201)

res <- m$run(sys0, t)
matplot(res$t, cbind(res$R, res$y), type="l", lty=1:2, xlab="Time",
        ylab="Abundance (red), Resource (black)",
        ylim=c(0, max(res$R, res$y)))
eq <- m$single_equilibrium(sys0$x)
abline(h=unlist(eq[c("R", "y")]), col=1:2, lty=3)

# Displace the solution from equilibrium and look at the new level of
# resources:
dy <- eq$y * 0.1

sys.y1 <- modifyList(eq, list(y=eq$y - dy))
sys.y2 <- modifyList(eq, list(y=eq$y + dy))

res.y1 <- m$run_fixed_density(sys.y1, t)
res.y2 <- m$run_fixed_density(sys.y2, t)

eq.y1 <- m$equilibrium_R(sys.y1)
eq.y2 <- m$equilibrium_R(sys.y2)

matplot(res.y1$t, cbind(res.y1$R, res.y2$R), type="l", col=2:3)
abline(h=c(eq$R, eq.y1$R, eq.y2$R), col=1:3, lty=3)
# Analytically:
points(max(res.y1$t), m$single_equilibrium_R(sys.y1), col=2)
points(max(res.y2$t), m$single_equilibrium_R(sys.y2), col=3)

# Next, we start working towards the instantaneous growth rate of a
# new type at this equilibrium

# Look at the fitness landscape: how does the instantaneous growth
# rate look with respect to K:
x.mutant <- seq(0, 1, length=101)
x.K <- rbind(x.mutant, eq$x[2], deparse.level=0) # K varying
x.C <- rbind(eq$x[2], x.mutant, deparse.level=0) # C varying

# Fitness does not vary with respect to c when rare, because species
# only consume resources when not rare.  Which doesn't mean that there
# will be no competition between species with different c values but
# the same K values of course: a different c value can drive a
# different species extinct!
plot(x.mutant, m$fitness(x.C, eq$x, eq$y, eq$R), type="l",
     xlab="Trait (C)", ylab="Fitness")
abline(h=0, col="grey", lty=3)

# Fitness does vary with K.
w.mutant <- m$fitness(x.K, eq$x, eq$y, eq$R)
plot(x.mutant, w.mutant, type="l",
     xlab="Trait (K)", ylab="Fitness")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

# But then there's a broader view of what competition is here: The
# species with K above the resident species are all strongly competed
# against by the resident species, regardless of what the derivatives
# show -- because they can't grow there.  That's the insight that we'd
# get from the subtraction case:
w.empty <- m$fitness(x.K, eq$x, 0, m$parameters$get()[["S"]])
plot(x.mutant, w.empty, type="l")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

# Here is how much fitness is suppressed by the presence of the
# resident type (dashed lines).
plot(x.mutant, w.mutant - w.empty, type="l")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

# Next, see what happpens to this fitness as we vary the population
# density.  This involves recomputing the new resource density for the
# species density:
z <- jacobian(function(y) m$fitness(x.K, eq$x, y), eq$y)

plot(x.mutant, z, type="l")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

op <- par(mfrow=c(3, 1), mar=c(2.5, 4, .5, .5))
plot(x.mutant, w.mutant - w.empty, type="l")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

plot(x.mutant, z, type="l")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

plot(x.mutant, z / w.empty, type="l")
abline(v=eq$x[1], lty=2)
par(op)

# To make this all a little clearer, this is what the fitness gradient
# is measuring.  Focus on a single mutant, whose trait is 0.3:
x.mutant1 <- rbind(0.3, eq$x[[2]])

# Then consider a range of densities of the resident species.
y.resident <- eq$y + seq(-1, 1, length=101)

# Fitness of the mutant as a function of these resident densities:
w.mutant1 <- sapply(y.resident, function(y) m$fitness(x.mutant1, eq$x, y))

w.mutant0 <- m$fitness(x.mutant1, eq$x, eq$y)

plot(y.resident, w.mutant1, type="l")
abline(h=w.mutant0, v=eq$y, lty=3, col="grey")
z.mutant1 <- jacobian(function(y)
                      m$fitness(x.mutant1, eq$x, eq$y), eq$y)
points(eq$y, w.mutant0, col="red")
abcline(eq$y, w.mutant0, z.mutant1, lty=2, col="red")
