# Doing the derivatives test is a bit harder for tree, because where
# there are strong empty fitness peaks the surface gets distorted
# badly.  This is interesting in itself probably.
library(Revolve)
library(deSolve)
library(tree)
library(numDeriv)
source("util-tree.R")

p <- new(Parameters)
# p$add_strategy(new(Strategy, list(lma=0.0648406)))
# p$add_strategy(new(Strategy, list(lma=0.1977910)))
p$add_strategy(new(Strategy, list(lma=0.0575703671036349)))
p$add_strategy(new(Strategy, list(lma=0.201106071562295)))
p$seed_rain <- c(1.1, 1.1)               # Starting rain.
p$seed_rain <- c(319.256840535597, 4.51466102145741)
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster

t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)
times0 <- cohort.introduction.times(t.max)
schedule0 <- schedule.from.times(times0, 2L)

# Work out what the equilibrium seed rain is.
res <- equilibrium.seed.rain(p, schedule0, 10,
                             build.args=list(verbose=TRUE),
                             progress=TRUE, verbose=TRUE)
# Set the seed rain back into the parameters
p$seed_rain <- res$seed_rain[,"out"]

# And rerun again to get the ODE times.
schedule <- res$schedule
schedule$ode_times <- attr(res, "ebt")$ode_times

# Resident lma values:
lma.res <- sapply(seq_len(p$size), function(i) p[[i]]$parameters$lma)

# Generate 101 log-spaced LMA values and compute mutant fitness over
# these values.
lma.mutant <- sort(c(seq_log(0.03, 0.8, 101), lma.res))
w.mutant <- landscape(lma.mutant, p, schedule)

# Seed production:
plot(w.mutant ~ lma.mutant, log="x", type="l")
abline(v=lma.res, col=1:2, lty=2)
abline(h=0, col="grey")

# On a log scale so that zero is no net increase:
plot(log(w.mutant) ~ lma.mutant, log="x", type="l")
abline(v=lma.res, col=1:2, lty=2)
abline(h=0, col="grey")

# Now, we want to work out what the derivatives of the fitness
# landscape are with respect to resident density, given these mutants.
z <- model_jacobian_density(lma.mutant, p, schedule)

# If what we care about is seed output on a log scale, then we can use
# the chain rule to compute that from this derivative and the fitness
# landscape:
#   d/dx [log(f(x))] -> f'(x) / f(x)
z.log <- z / w.mutant

# Here is the competitive effect.  Larger negative numbers are
# stronger competition.  Positive numbers are facilitation!
matplot(lma.mutant, z, type="l", log="x")
abline(v=lma.res, lty=3, col=1:2)

matplot(lma.mutant, z.log, type="l", log="x")
abline(v=lma.res, lty=3, col=1:2)

# These are quite different landscapes than either what I'd expect, or
# what we see in the toy models case.  First, there are regions of
# apparent facilitation (at least on a log scale, look at the black
# line and see it increase to the left).

# There are two potential sources of distortion here:

# 1. numerical inaccuracies, particularly in the relative calculations
#    of mutant fitness when very low (cohort refining stops around
#    changes of about 1e-3).
#
# 2. massive changes in the carrying capacity.

# Other things to try:
#
# * how do the landscapes change with population size (i.e., out of
#   demographic equilibrium)?
# * how does the seed rain of an *established* species vary with the
#   changing density? (that requires finding the equilibrium
#   population size, but for small enough changes we could probably
#   just take the output seed rain).
