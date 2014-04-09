# Daniel suggested looking at how these derivative approaches vary
# with population size.
library(Revolve)
library(deSolve)
library(tree)
library(numDeriv)
library(parallel)
library(RColorBrewer)
source("util-tree.R")

competition_assay <- function(N, lma.mutant, p) {
  t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)
  times0 <- cohort.introduction.times(t.max)
  schedule0 <- schedule.from.times(times0, 1L)

  p <- p$copy()
  p$seed_rain <- N
  schedule <- build.schedule(p, schedule0, 5, 1e-3,
                             progress=FALSE, verbose=FALSE)
  schedule$ode_times <- attr(schedule, "ebt")$ode_times

  w <- landscape(lma.mutant, p, schedule)
  z <- model_jacobian_density(lma.mutant, p, schedule)
  list(lma=lma.mutant, w=w, z=z)
}

p <- new(Parameters)
p$add_strategy(new(Strategy, list(lma=0.0578025395245511)))
p$seed_rain <- c(1.1)               # Starting rain.
p$seed_rain <- c(319.256840535597)
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast.control()) # A bit faster

lma.res <- sapply(seq_len(p$size), function(i) p[[i]]$parameters$lma)
lma.mutant <- sort(c(seq_log(0.03,  0.8,   101),
                     seq_log(0.055, 0.061, 51),
                     lma.res))

# Equilibriun seed rain for this species is about 323.5, depending on
# exactly how the cohort schedule ends up structured (yes, that's
# quite unfortunate).
NN <- 323.5 * 2^seq(-4, 1)

# This is really slow (several CPU hours)
ans <- mclapply(NN, competition_assay, lma.mutant, p,
                mc.preschedule=FALSE)

cols <- names(ans[[1]])
tmp <- lapply(cols, function(v) sapply(ans, "[[", v))
names(tmp) <- cols

cols <- brewer.pal(length(NN), "Blues")

matplot(tmp$lma, tmp$w, type="l", lty=1, col=cols, log="x")
abline(v=lma.res)
abline(h=1, col="grey", lty=3)

matplot(tmp$lma, tmp$w, type="l", lty=1, col=cols, log="xy")
abline(v=lma.res)
abline(h=1, col="grey", lty=3)

matplot(tmp$lma, tmp$z, type="l", lty=1, col=cols, log="x")
abline(v=lma.res)
abline(h=0, col="grey", lty=3)

matplot(tmp$lma, tmp$z / tmp$w, type="l", lty=1, col=cols, log="x")
abline(v=lma.res)
abline(h=0, col="grey", lty=3)
