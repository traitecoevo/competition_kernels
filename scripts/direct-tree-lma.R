library(Revolve)
library(deSolve)
library(tree)
library(parallel)

source("direct-tree-fun.R")

## Range of values to work with for now:
lma.r <- c(0.03, 3.8)

## Base set of 20 parameters for which we'll compute r and K.  This
## should be close to enough?
lma.base <- seq_log(lma.r[[1]], lma.r[[2]], 20)
## More complete set of points for everything else:
lma.full <- seq_log(lma.r[[1]], lma.r[[2]], 101)

r.base <- log(max_growth_rate("lma", lma.base))
K.base <- carrying_capacity("lma", lma.base, parallel=TRUE)

r.i <- loglog_splinefun(lma.base, r.base, "x")
K.i <- loglog_splinefun(lma.base, K.base, "xy")

## Maximum growth rate
plot(lma.base, r.base, log="x")
lines(lma.full, r.i(lma.full))

## Carrying capacity
plot(lma.base, K.base, log="x")
lines(lma.full, K.i(lma.full))

## Resident species:
lma.r <- 0.2
N.r <- K.i(lma.r)

## Generating function that will give us invasion fitness 
f.w <- make_invasion_fitness("lma", lma.r, N.r)
w.i <- f.w(lma.full)

plot(lma.full, w.i, log="x", type="l")
abline(h=0, lty=2)
abline(v=lma.r)

alpha.i <- compute_alpha(w.i, r.i(lma.full), K.i(lma.full), N.r)

res.lma <- list(trait="lma",
                x.resident=lma.r,
                N.resident=N.r,
                base=lma.base,
                full=lma.full,
                r.base=r.base,
                K.base=K.base,
                w=w.i,
                r=r.i, K=K.i,
                alpha=alpha.i)

saveRDS(res.lma, "res.lma.rds")

plot(lma.full, alpha.lma,
     log="x", type="l")
abline(v=lma.r, lty=2)

## Function to stick that all together for us:
f <- function(lma.r, lma.i, r.i, K.i) {
  N.r <- K.i(lma.r)
  ret <- list(lma.r=lma.r, N.r=N.r,
              lma.i=lma.i, r.i=r.i(lma.i), K.i=K.i(lma.i))
  ret$w.i <- make_invasion_fitness(lma.r, N.r)(lma.i)
  ret$alpha <- compute_alpha(ret$w.i, ret$r.i, ret$K.i, N.r)
  ret
}

## Generate a series of potentially interesting values and look at the
## fitness landscape and competition function over it:
lma.r <- seq_log(.06, .2, 10)
ans <- mclapply(lma.r, f, lma.full, r.i, K.i, mc.preschedule=FALSE)

wlim <- range(sapply(ans, function(x) range(x$w.i)))
alim <- range(sapply(ans, function(x) range(x$alpha)))
par(mfcol=c(2, length(lma.r)), mar=rep(0, 4), oma=c(c(2.5, 2.5, .5, .5)))
for (el in ans) {
  plot(w.i ~ lma.i, el, type="l", log="x", axes=FALSE, ylim=wlim)
  box()
  abline(v=el$lma.r, h=0, lty=2)
  if (identical(el, ans[[1]])) {
    mtext("Fitness", 2, line=1)
  }

  plot(alpha ~ lma.i, el, type="l", log="x", axes=FALSE, ylim=alim)
  box()
  abline(v=el$lma.r, h=c(0, 1), lty=2)
  mtext(round(el$lma.r, 2), 1, line=1)
  if (identical(el, ans[[1]])) {
    mtext("Alpha", 2, line=1)
  }
}

## Then check what is going on right around the fitness peak:
## Little function to zoom in on the maximum:
g <- function(lma.r, r.i, K.i, eps=1e-4) {
  N.r <- K.i(lma.r)
  lma.i <- lma.r + c(-1, 1) * eps
  w <- make_invasion_fitness(lma.r, N.r)
  ww <- w(lma.i)
  gr <- diff(ww) / (2 * eps)
  attr(gr, "w") <- ww
  attr(gr, "lma") <- lma.i
  gr
}

x0 <- .07
x1 <- .1
y0 <- g(x0, r.i, K.i)
y1 <- g(x1, r.i, K.i)
opt <- uniroot(g, c(x0, x1), r.i, K.i, f.lower=y0, f.upper=y1)

## Add some detail around the resident:
i <- which.min(abs(lma.full - opt$root))
lma.extra <- seq_log(lma.full[i - 2], lma.full[i + 2], length=31)
res.root <- f(opt$root, sort(c(lma.full, lma.extra)), r.i, K.i)

plot(res.root$lma.i, res.root$w.i, log="x", type="l",
     xlab="LMA", ylab="Fitness")
abline(v=res.root$lma.r, lty=2)

## Here's a zoom on the detailed region:
plot(res.root$lma.i, res.root$w.i, log="x", type="l",
     xlab="LMA", ylab="Fitness",
     xlim=range(lma.extra),
     ylim=range(res.root$w.i[res.root$lma.i >= min(lma.extra) &
       res.root$lma.i <= max(lma.extra)]))
abline(v=res.root$lma.r, lty=2)

plot(res.root$lma.i, res.root$alpha, log="x", type="l",
     xlab="LMA", ylab="Competition (resident on invader)")
abline(v=res.root$lma.r, lty=2)
