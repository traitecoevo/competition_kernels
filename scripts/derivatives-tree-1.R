library(Revolve)
library(deSolve)
library(tree)
library(numDeriv)
source("util-tree.R")

p <- new(Parameters)
p$add_strategy(new(Strategy, list(lma=0.0578025395245511)))
p$seed_rain <- c(1.1)               # Starting rain.
p$seed_rain <- c(319.256840535597)
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster

t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)
times0 <- cohort.introduction.times(t.max)
schedule0 <- schedule.from.times(times0, 1L)

## Work out what the equilibrium seed rain is.
res <- equilibrium.seed.rain(p, schedule0, 10,
                             build.args=list(verbose=TRUE),
                             progress=TRUE, verbose=TRUE)
## Set the seed rain back into the parameters
p$seed_rain <- res$seed_rain[,"out"]

## Extract the ODE times so that we run at the same points as the
## resident population was run at.
schedule <- res$schedule
schedule$ode_times <- attr(res, "ebt")$ode_times

## Resident lma values:
lma.res <- sapply(seq_len(p$size), function(i) p[[i]]$parameters$lma)

## Generate 101 log-spaced LMA values and compute mutant fitness over
## these values.
lma.mutant <- sort(c(seq_log(0.03, 0.8, 101), lma.res))
w.mutant <- landscape(lma.mutant, p, schedule)

plot(w.mutant ~ lma.mutant, log="x", type="l")
abline(v=lma.res, col=1:2, lty=2)
abline(h=1, col="grey")

plot(w.mutant ~ lma.mutant, log="xy", type="l")
abline(v=lma.res, col=1:2, lty=2)
abline(h=1, col="grey")

## Zoom in on the area with the resident: looks like we need a higher
## resolution to see what is going on here.
plot(w.mutant ~ lma.mutant, log="xy", type="l", xlim=c(0.054, 0.063))
abline(v=lma.res, col=1:2, lty=2)
abline(h=1, col="grey")
abline(v=c(0.055, 0.061), col="red")

lma.mutant.detail <- seq_log(0.055, 0.061, length=51)
w.mutant.detail <- landscape(lma.mutant.detail, p, schedule)

## This looks bad, but it's just linear interpolation:
plot(w.mutant ~ lma.mutant, log="xy", type="l", xlim=c(0.054, 0.063))
lines(w.mutant.detail ~ lma.mutant.detail, col="red")
abline(v=lma.res, col=1:2, lty=2)
abline(h=1, col="grey")

## Stick these all together:
lma.mutant.all <- c(lma.mutant, lma.mutant.detail)
w.mutant.all <- c(w.mutant, w.mutant.detail)[order(lma.mutant.all)]
lma.mutant.all <- sort(lma.mutant.all)

plot(w.mutant.all ~ lma.mutant.all, log="xy", type="l")
abline(v=lma.res, col=1:2, lty=2)
abline(h=1, col="grey")

## One interesting thing here is that you can see the underlying
## fitness landscape looks like two curves superimposed on each other:
## one peak at the resident species, the other a much broader peak.
plot(w.mutant.all ~ lma.mutant.all, log="xy", type="l", xlim=c(0.054, 0.063))
abline(v=lma.res, col=1:2, lty=2)
abline(h=1, col="grey")

## We'll just use the detailed version below here:
lma.mutant <- lma.mutant.all
w.mutant <- w.mutant.all
rm(w.mutant.all, lma.mutant.all)

## See how the current fitness landscape relates to the empty
## landscape (this is going to give ridiculous values, of course).
w.empty <- landscape.empty(lma.mutant, p, schedule)
plot(w.empty ~ lma.mutant, log="xy", type="l")

## We could try and superimpose these on one another by rescaling a
## peak height arbitrarily to 1, and see if they're the same shape,
## perhaps.

## This makes me think of several aspects of competition: there is a
## clear effect over the entire landscape; the number of seeds produced
## is depressed extremely strongly by the presence of a large number of
## individuals.  But is the per-capita rate of seed production in an
## empty environment even a slightly fair comparison?

## This will also imply an asymmetrical competition function because
## the resident species is on one side of the empty landscape.  So for
## symmetry we might want to see how a non-optimal strategy with high
## LMA will do.  Unfortunately, being on a non-flat part of the fitness
## landscape might cause some issues.

## Now, we want to work out what the derivatives of the fitness
## landscape are with respect to resident density, given these mutants.

## This is quite slow to evaluate.
z <- model_jacobian_density(lma.mutant, p, schedule)

## Derivative of the log of fitness, via the chain rule:
z.log <- z / w.mutant

## Plotted on a non-log scale, ther is a small region of competitive
## effect: it's fairly concentrated around the focal species as a
## region of depression.  That is increase the density of the
## single resident species and the seed production of nearby strategies
## is decreased.

## Then, some distance away at the next peak, there is a much larger
## competitive effect.
plot(lma.mutant, z, type="l", log="x")
abline(v=lma.res, lty=3, col=1:2)

plot(lma.mutant, z, type="l", log="x", ylim=c(-0.02, 0.01))
abline(v=lma.res, lty=3, col=1)
abline(h=0, lty=3, col="grey")

## This version is possibly more meaningful?
plot(lma.mutant, z.log, type="l", log="x")
abline(v=lma.res, lty=3, col=1:2)

## Arond the resident:
plot(lma.mutant, z, type="l", log="x", ylim=c(-0.015, 0.01),
     xlim=c(0.054, 0.063))
abline(v=lma.res, lty=3, col=1)
abline(h=0, lty=3, col="grey")

plot(lma.mutant, z.log, type="l", log="x",
     xlim=c(0.054, 0.063))
abline(v=lma.res, lty=3, col=1)
abline(h=0, lty=3, col="grey")

## Start a single strategy over where the peak is in the empty landscape
p2 <- p$copy()
p2$clear()
p2$add_strategy(new(Strategy, list(lma=lma.mutant[which.max(w.empty)])))

res2 <- equilibrium.seed.rain(p2, schedule0, 10,
                              build.args=list(verbose=TRUE),
                              progress=TRUE, verbose=TRUE)
## Set the seed rain back into the parameters
p2$seed_rain <- res2$seed_rain[,"out"]

## Extract the ODE times so that we run at the same points as the
## resident population was run at.
schedule2 <- res2$schedule
schedule2$ode_times <- attr(res2, "ebt")$ode_times

## Resident lma values:
lma.res2 <- sapply(seq_len(p2$size), function(i) p2[[i]]$parameters$lma)

lma.mutant2 <- sort(c(seq_log(0.03, 0.8, 101), lma.res2))
w.mutant2 <- landscape(lma.mutant2, p2, schedule2)

plot(w.mutant2 ~ lma.mutant2, log="x", type="l")
abline(v=lma.res2, lty=2)
abline(h=0, col="grey")

plot(log(w.mutant2) ~ lma.mutant2, log="x", type="l")
abline(v=lma.res2, lty=2)
abline(h=0, col="grey")

z2 <- model_jacobian_density(lma.mutant2, p2, schedule2)
z.log2 <- z2 / w.mutant2

plot(lma.mutant2, z2, type="l", log="x")
abline(v=lma.res2, lty=3, col=1)
abline(h=0, lty=3, col="grey")

## This version is possibly more meaningful?
plot(lma.mutant2, z.log2, type="l", log="x")
abline(v=lma.res2, lty=3, col=1:2)
abline(h=0, lty=3, col="grey")
