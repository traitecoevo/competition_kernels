library(Revolve)
library(deSolve)
library(tree)
library(numDeriv)

p <- new(Parameters)
# p$add_strategy(new(Strategy, list(lma=0.0648406)))
# p$add_strategy(new(Strategy, list(lma=0.1977910)))
p$add_strategy(new(Strategy, list(lma=0.0575703671036349)))
p$seed_rain <- c(1.1)               # Starting rain.
p$seed_rain <- c(319.256840535597)
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster

t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)
times0 <- cohort.introduction.times(t.max)
schedule0 <- schedule.from.times(times0, 1L)

# Work out what the equilibrium seed rain is.
res <- equilibrium.seed.rain(p, schedule0, 10,
                             build.args=list(verbose=TRUE),
                             progress=TRUE, verbose=TRUE)
# Set the seed rain back into the parameters
p$seed_rain <- res$seed_rain[,"out"]

# And rerun again to get the ODE times.
schedule <- res$schedule
ebt <- run.ebt(p, schedule$copy())
schedule$ode_times <- ebt$ode_times

# Resident lma values:
lma.res <- sapply(seq_len(p$size), function(i) p[[i]]$parameters$lma)

# Generate 101 log-spaced LMA values and compute mutant fitness over
# these values.
lma.mutant <- sort(c(seq_log(0.03, 0.8, 101), lma.res))
w.mutant <- landscape(lma.mutant, p, schedule)

plot(w.mutant ~ lma.mutant, log="x", type="l")
abline(v=lma.res, col=1:2, lty=2)
abline(h=0, col="grey")

plot(w.mutant ~ lma.mutant, log="xy", type="l")
abline(v=lma.res, col=1:2, lty=2)
abline(h=0, col="grey")


# Now, we want to work out what the derivatives of the fitness
# landscape are with respect to resident density, given these mutants.
# That's not too horrid to do, but we need a support function:
target <- function(N, lma, p, schedule) {
  p <- p$copy()
  p$seed_rain <- N
  landscape(lma, p, schedule)
}

# This is quite slow to evaluate.  Each evaluation takes about 4s, so
# that's 400s (a little over 6 minutes).
z <- jacobian(function(N) target(N, lma.mutant, p, schedule),
              res$seed_rain[,"out"], method="simple")
