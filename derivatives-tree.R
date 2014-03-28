library(Revolve)
library(deSolve)
library(tree)
library(numDeriv)

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
ebt <- run.ebt(p, schedule$copy())
schedule$ode_times <- ebt$ode_times

# Resident lma values:
lma.res <- sapply(seq_len(p$size), function(i) p[[i]]$parameters$lma)

# Generate 101 log-spaced LMA values and compute mutant fitness over
# these values.
lma.mutant <- sort(c(seq_log(0.03, 0.8, 101), lma.res))
w.mutant <- landscape(lma.mutant, p, schedule)

# I'm pretty sure that the 
plot(w.mutant ~ lma.mutant, log="x", type="l")
abline(v=lma.res, col=1:2, lty=2)
abline(h=0, col="grey")

pdf("fitness_tree.pdf")
plot(log(w.mutant) ~ lma.mutant, log="x", type="l")
abline(v=lma.res, col=1:2, lty=2)
abline(h=0, col="grey")
dev.off()

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

# If what we care about is seed output on a log scale, then we can use
# the chain rule to compute that from this derivative and the fitness
# landscape:
#   d/dx [log(f(x))] -> f'(x) / f(x)
z.log <- z / w.mutant

# Here is the competitive effect.  Larger negative numbers are
# stronger competition.  Positive numbers are facilitation!
pdf("competition_tree.pdf")
matplot(lma.mutant, z, type="l", log="x")
abline(v=lma.res, lty=3, col=1:2)

matplot(lma.mutant, z.log, type="l", log="x")
abline(v=lma.res, lty=3, col=1:2)
dev.off()


plot(w.mutant ~ lma.mutant, log="x", type="l")
abline(v=lma.res, col=1:2, lty=2)
abline(h=0, col="grey")

matplot(lma.mutant, z, type="l", log="x", ylim=c(-0.001, 0.001))
abline(v=lma.res, lty=3, col=1:2)

matplot(lma.mutant, z, type="o", log="x", ylim=c(-0.0001, 0.0001),
        pch=19, cex=.3)
abline(v=lma.res, lty=3, col=1:2)

matplot(lma.mutant, sign(z), type="s", log="x",
        pch=19, cex=.3)
abline(v=lma.res, lty=3, col=1:2)




z2 <- jacobian(function(N) log(target(N, lma.mutant, p, schedule)),
               res$seed_rain[,"out"], method="simple")


matplot(lma.mutant, z2, type="l", log="x")
abline(v=lma.res, lty=3, col=1:2)

matpoints(lma.mutant, z / w.mutant, pch=1)





















lma.64 <- seq.log(0.03, 0.8, 64)
w.64 <- landscape(lma.64, p, schedule)

## Note that this should actually be w - 1 I think (input seed rain is
## arbitrarily 1 for mutants). Check with Daniel.  Doesn't affect the
## derivative calculations though.
pdf("landscape_tree.pdf")
plot(w.64 - 1 ~ lma.64, type="l", log="x")
abline(v=lma.res, lty=3, col=1:2)
abline(h=0)
dev.off()

## TODO: This is where I am up to!

## Oh dear - something is up here.  The mutant fitness calculations
## are off for one of the residents.  It should come back with one
## seed out, and I thought I'd checked that before.  Something is
## terribly up.

## Now, we want to work out what the derivatives of the fitness
## landscape are with respect to resident density, given these
## mutants.  That's not too horrid to do.
f <- function(N, lma, p, schedule) {
  p <- p$copy()
  p$seed_rain <- N
  landscape(lma, p, schedule)
}
h <- function(N)
  f(N, lma.64, p, schedule)

tmp <- f(res$seed_rain[,"out"], lma.64, p, schedule)

z <- jacobian(h, res$seed_rain[,"out"], method="simple")

pdf("competition_tree.pdf")
matplot(lma.64, z, type="l", log="x")
abline(v=lma.res, lty=3, col=1:2)
dev.off()

## OK, so this looks stupid, because things aren't at the top of their
## respective hills.  The black point should be left a little bit and
## the red point should be right a little bit.  And there possibly
## should be a third type at around lma 0.50.

plot(w.64 - 1 ~ lma.64, type="l", log="x")
abline(v=lma.res, lty=3, col=1:2)
abline(h=0)

## So, work out where the peak shoud be.  To do this, find all the
## "peaks".  Those are the times where the first derivative changes
## sign from negaitve to positive.  So the optimim is bracketed within
## the points here:
i <- which(diff(sign(diff(w.64))) == -2)
j <- i + 1L
plot(w.64 ~ lma.64, log="x", type="l")
points(w.64[i] ~ lma.64[i], col="red", pch=19)
points(w.64[j] ~ lma.64[j], col="blue", pch=19)
abline(v=lma.res, lty=3, col=1:2)

## Though it might be that as we move that way the optimum will have
## moved.  I really don't want to get involved with a full-on assembly
## here.
lma.res2 <- ((lma.64[i] + lma.64[j])/2)[1:2]
p2 <- new(Parameters)
for (i in lma.res2)
  p2$add_strategy(new(Strategy, list(lma=i)))
p2$seed_rain <- c(348.738300488797, 3.07449195531246) # Starting rain.
p2$set_parameters(list(patch_area=1.0))   # See issue #13
p2$set_control_parameters(fast.control()) # A bit faster
schedule0 <- schedule.from.times(times0, 2L)

res <- equilibrium.seed.rain(p2, schedule0, 10,
                             build.args=list(verbose=TRUE),
                             progress=TRUE, verbose=TRUE)
p2$seed_rain <- res$seed_rain[,"out"]

# Resident lma values:
lma.res <- sapply(seq_len(p2$size), function(i) p2[[i]]$parameters$lma)

## Forcing the resident types in here.  This would be much nicer with
## an adaptive refinement such as Conrad's GP code.
lma.mutant <- sort(c(seq.log(0.03, 0.8, 128), lma.res))
w.mutant <- landscape(lma.mutant, p2, schedule)
plot(w.mutant - 1 ~ lma.mutant, log="x", type="l")
abline(v=lma.res, col=1:2, lty=2)


h <- function(N)
  f(N, lma.mutant, p2, schedule)
z <- jacobian(h, res$seed_rain[,"out"], method="simple")

plot(w.mutant - 1 ~ lma.mutant, type="l", log="x")
abline(v=lma.res, lty=3, col=1:2)
abline(h=0, col="grey", lty=2)
## These are quite different magnitudes, but some of this is probably
## explained by the fact that the first species 
matplot(lma.mutant, z, type="l", log="x", lty=1)
abline(v=lma.res, lty=3, col=1:2)

## Looking close there is a weird faciliation effect at lma values
## very slightly smaller than the first species.  There is also a
## wobble around the second species.  This exactly goes through zero
## at the location of the second species, so I think that's a second
## order effect perhaps?  If we drop the red species, this should
## disappear.  But while that species is here -- you increase the
## abundance of the black species, that decreases the abundance of the
## red species allowing more growth here.  Not sure what is going on
## at the point slightly below the red species though.
matplot(lma.mutant, z, type="l", log="x", lty=1, ylim=c(-0.02, 0.001))
abline(v=lma.res, lty=3, col=1:2)
abline(h=0, col="grey")

## Next steps:
##   1. Tidy up
##      * landscape code and default heights, etc
##      * optimisation code
##      * definition of fitness
##      * quick-and-dirty interpolation to skip the boring bits
##   2. Two dimensions
##      * possibly need the gp code here?
