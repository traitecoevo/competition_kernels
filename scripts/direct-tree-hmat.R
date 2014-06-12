library(Revolve)
library(deSolve)
library(tree)
library(parallel)

source("direct-tree-fun.R")

s0 <- new(Strategy)

## Range of values to work with for now.  Not sure why hmat=1 is
## failing, but does so with odd error message (height must be
## positive) which suggests that something nasty has happened
## somewhere.
hmat.r <- c(1, 50)

## Base set of 20 parameters for which we'll compute r and K.  This
## should be close to enough?
hmat.base <- seq_log(hmat.r[[1]], hmat.r[[2]], 20)
## More complete set of points for everything else:
hmat.full <- seq_log(hmat.r[[1]], hmat.r[[2]], 101)

r.base <- log(max_growth_rate("hmat", hmat.base))
K.base <- carrying_capacity("hmat", hmat.base, parallel=TRUE)

r.i <- loglog_splinefun(hmat.base, r.base, "x")
K.i <- loglog_splinefun(hmat.base, K.base, "xy")

## Maximum growth rate
plot(hmat.base, r.base, log="x")
lines(hmat.full, r.i(hmat.full))

## Carrying capacity
plot(hmat.base, K.base, log="x")
lines(hmat.full, K.i(hmat.full))

## Resident species:
hmat.r <- 16
N.r <- K.i(hmat.r)

## Generating function that will give us invasion fitness 
f.w <- make_invasion_fitness("hmat", hmat.r, N.r)
w.i <- f.w(hmat.full)

alpha.i <- compute_alpha(w.i, r.i(hmat.full), K.i(hmat.full), N.r)

plot(hmat.full, w.i, log="x", type="l")
abline(h=0, lty=2)
abline(v=hmat.r)

res.hmat <- list(trait="hmat",
                 x.resident=hmat.r,
                 N.resident=N.r,
                 base=hmat.base,
                 full=hmat.full,
                 r.base=r.base,
                 K.base=K.base,
                 w=w.i,
                 alpha=alpha.i,
                 r=r.i, K=K.i)
saveRDS(res.hmat, "../talk/res.hmat.rds")
