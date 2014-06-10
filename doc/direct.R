source("common.R")

## See the `direct.tex` file for what is going on here.

## So we don't need the derivative at all - we can just plug things
## directly back.

## Start with the "full" single species model.  Each species is
## defined by a vector of two elements: {K, C}.
mat_r1 <- rstar_matrices(rstar_mat_1, rstar_mat_1)
m_r1 <- rstar(mat_r1, S=1)

## First consider an environment with a single species in it, with
## parameters {0.5, 0.5}, at equilibriu,
eq <- m_r1$single_equilibrium(matrix(0.5, nrow=2))

## Then consider a range of mutants with different K or different C
## values.  The species with the smallest K value will win.
x.mutant <- seq(0, 1, length=101)[-1]
x.K <- rbind(x.mutant, eq$x[2], deparse.level=0)
x.C <- rbind(eq$x[1], x.mutant, deparse.level=0)

## We can compute the fitness (per capita rate of increase).  Fitness
## increases as species have lower K values.
## rstar_1_fitness_K
w.mutant.K <- m_r1$fitness(x.K, eq$x, eq$y, eq$R)
plot(x.mutant, w.mutant.K, type="l",
     xlab="Trait (K)", ylab="Fitness")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

## In contrast, fitness does not depend on C for invasion (given all
## invaders have the same K value):
## rstar_1_fitness_C
w.mutant.C <- m_r1$fitness(x.C, eq$x, eq$y, eq$R)
plot(x.mutant, w.mutant.C, type="l",
     xlab="Trait (K)", ylab="Fitness")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

## Then compute the initial maximum growth rate and carrying capacity
## for the invading species.
r.mutant.K <- m_r1$max_growth_rate(x.K)
K.mutant.K <- m_r1$carrying_capacity(x.K)
r.mutant.C <- m_r1$max_growth_rate(x.C)
K.mutant.C <- m_r1$carrying_capacity(x.C)

## Maximum growth rate:
plot(x.mutant, r.mutant.K, type="l")
lines(x.mutant, r.mutant.C, type="l", col="blue")

## Carrying capacity:
plot(x.mutant, K.mutant.K, type="l")
lines(x.mutant, K.mutant.C, type="l", col="blue")

## And the density of the resident species.
N.resident <- eq$y

compute_alpha <- function(w.i, r.i, K.i, N.r, type=1) {
  if (type == 1) {
    (K.i - 1 / r.i * w.i) / N.r
  } else if (type == 2) {
    (1 - 1 / r.i * w.i) * K.i / N.r
  } else {
    stop("no such type")
  }
}

alpha1.K <- compute_alpha(w.mutant.K, r.mutant.K, K.mutant.K,
                          N.resident, 1)
alpha2.K <- compute_alpha(w.mutant.K, r.mutant.K, K.mutant.K,
                          N.resident, 2)

alpha1.C <- compute_alpha(w.mutant.C, r.mutant.C, K.mutant.C,
                          N.resident, 1)
alpha2.C <- compute_alpha(w.mutant.C, r.mutant.C, K.mutant.C,
                          N.resident, 2)

## So, above 1 competition is more extreme than for 
plot(x.mutant, alpha2.K, type="l", xlab="Invader trait (K)",
     ylab="Competitive effect")
lines(x.mutant, alpha1.K, type="l", lty=2)
abline(h=c(0, 1), v=eq$x[[1]], lty=3)

plot(x.mutant, alpha2.C, type="l", xlab="Invader trait (C)",
     ylab="Competitive effect")
lines(x.mutant, alpha1.C, type="l", lty=2)
abline(h=c(0, 1), v=eq$x[[1]], lty=3)

## We can also compute the full trait-wide competition landscape and
## fitness landscape:
xx <- unname(t(as.matrix(expand.grid(K=x.mutant, C=x.mutant))))
ww <- matrix(m_r1$fitness(xx, eq$x, eq$y, eq$R), length(x.mutant))
rr <- m_r1$max_growth_rate(xx)
kk <- m_r1$carrying_capacity(xx)
a1 <- compute_alpha(ww, rr, kk, N.resident, 1)
a2 <- compute_alpha(ww, rr, kk, N.resident, 2)
a1[!is.finite(a1)] <- a2[!is.finite(a2)] <- NA

persp(x.mutant, x.mutant, ww, xlab="K", ylab="C", zlab="w",
      theta=30, shade=.2, border=NA, col="dodgerblue4")
## OK, this is really surprising.  The different alpha calculations
## don't give the same shape
persp(x.mutant, x.mutant, a1, xlab="K", ylab="C", zlab="alpha (1)",
      theta=30, shade=.2, border=NA, col="dodgerblue4")
persp(x.mutant, x.mutant, a2, xlab="K", ylab="C", zlab="alpha (2)",
      theta=30, shade=.2, border=NA, col="dodgerblue4")

alpha <- function(x, m, sys, type=1) {
  r.i <- m$max_growth_rate(x)
  K.i <- m$carrying_capacity(x)
  w.i <- m$fitness(x, sys$x, sys$y, sys$R)
  N.r <- sys$y
  compute_alpha(w.i, r.i, K.i, N.r, type)
}

## So for relatively low C values, we get a flip in the trait slope of
## alpha.  Oh dear.
x.K1 <- rbind(x.mutant, .1, deparse.level=0)
plot(x.mutant, alpha(x.K1, m_r1, eq, 2), type="l")
lines(x.mutant, alpha(x.K1, m_r1, eq, 1), lty=2)

## Also note that the height of competition here has increased above
## 1; this is OK because the resident trait *does not appear in this
## graph*: we're looking elsewhere in C space.

## We should alway have intraspecific alpha = 1 because we have
##
## (1 - w/r) * K / Nr
##
## And K == Nr so
##
## (1 - w/r)
##
## and at equilibrium w = 0 so we have 1.

## Invasion depends only on the K value and not on the C value (see
## the fitness 3d plot above.

## Based on this, it seems that the solid line form (alpha2) is
## perhaps better suited to things here; it doesn't seem to vary with
## C, except as a linear rescaling (so no change to *shape*, just to
## the intensity of competition).  This might not be so obc

## Next, let's see how this varies with density of the resident. To do
## this we need to let the resources equilibrate at the level we'd
## expect if we subsidise/penalise the resident population.

eq$R
m_r1$single_equilibrium_R(eq)
m_r1$single_equilibrium_R(eq, TRUE)
m_r1$equilibrium_R(eq, method="nleqslv")

## Now, halve the population size:
tmp <- modifyList(eq, list(y=eq$y * 0.5))
tmp$R <- m_r1$single_equilibrium_R(tmp)

tt <- seq(0, 30, length=101)
plot(R ~ t, m_r1$run_fixed_density(tmp, tt), type="l")
abline(h=m_r1$single_equilibrium_R(tmp), lty=3)

## Now, vary the density and see how the estimates of competition
## change:
x.K1 <- rbind(x.mutant, .5, deparse.level=0)
plot(x.mutant, alpha(x.K1, m_r1, tmp, 2), type="l")
lines(x.mutant, alpha(x.K1, m_r1, tmp, 1), lty=2)

## To do this efficiently we might need to have a system for computing
## the equilibrium resource levels as a function of density, but for
## now, let's just do that outselves.
N.r.varying <- eq$y * 2^seq(-6, 2, length=51)
## Then compute the resource level at each of these:
f <- function(N, sys) {
  sys$y <- N
  sys$R <- m_r1$single_equilibrium_R(sys)
  sys
}
sys.varying <- lapply(N.r.varying, f, eq)

m1 <- sapply(sys.varying, function(sys) alpha(x.K1, m_r1, sys, 1))
m2 <- sapply(sys.varying, function(sys) alpha(x.K1, m_r1, sys, 2))

## OK, so that's an interesting shape, roughly similar to what I saw
## from the derivatives:
persp(x.mutant, log(N.r.varying), m2, theta=60,
      xlab="Trait (K)", ylab="Resident density (log scale)",
      zlab="Competition", shade=.2, border=NA, col="dodgerblue4")

## Same data, as an image plot:
image(x.mutant, N.r.varying, m2,
      xlab="Trait (K)", ylab="Resident density (log scale)",
      log="y")

## So we have competition that is nonlinear in the resident density to
## the point of being *non-monotonic*, the same basic result as using
## the derivatives.  In fact this plot is almost the same.  So that's
## cool.  Multiplying this through by the density of the resident then
## gives the total competitive impact on the invader.

## In contrast, the other version shows quite different results:
persp(x.mutant, log(N.r.varying), m1, theta=60,
      xlab="Trait (K)", ylab="Resident density (log scale)",
      zlab="Competition", shade=.2, border=NA, col="dodgerblue4")

## Same data, as an image plot:
image(x.mutant, N.r.varying, m1,
      xlab="Trait (K)", ylab="Resident density (log scale)",
      log="y")

## Is this just because of the Monod equation?
## 
## Let's rip that out.  Instead of growth being p(R) -> r R / (K + R)
## let's have p(R) -> r R / K (this choice being so that the
## sign-relatedness of growth rate and K remains the same).

## We should be able to do that by inheritance, but R's reference
## classes don't seem to like me.  It would actually break a few
## places, so I'm going to leave that be for a little bit.

## # Two resources

## Using the tradeoff version again
mat_r2 <- rstar_matrices(rstar_mat_2_tradeoff,
                         rstar_mat_2_tradeoff)
m_r2 <- rstar(mat_r2, S=c(0.5, 0.5))
eq <- m_r2$single_equilibrium(matrix(0.3, nrow=2))

x.K <- rbind(x.mutant, eq$x[2], deparse.level=0)
x.C <- rbind(eq$x[1], x.mutant, deparse.level=0)

## I need to add in the special points to get these lines nice and
## smooth.
w.mutant.K <- m_r2$fitness(x.K, eq$x, eq$y, eq$R)
plot(x.mutant, w.mutant.K, type="l",
     xlab="Trait (K)", ylab="Fitness")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

w.mutant.C <- m_r2$fitness(x.C, eq$x, eq$y, eq$R)
plot(x.mutant, w.mutant.C, type="l",
     xlab="Trait (K)", ylab="Fitness")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

r.mutant.K <- m_r2$max_growth_rate(x.K)
K.mutant.K <- m_r2$carrying_capacity(x.K)
r.mutant.C <- m_r2$max_growth_rate(x.C)
K.mutant.C <- m_r2$carrying_capacity(x.C)

## Maximum growth rate:
plot(x.mutant, r.mutant.K, type="l")
lines(x.mutant, r.mutant.C, type="l", col="blue")

## Carrying capacity:
plot(x.mutant, K.mutant.K, type="l")
lines(x.mutant, K.mutant.C, type="l", col="blue")

## And the density of the resident species.
N.resident <- eq$y

alpha1.K <- compute_alpha(w.mutant.K, r.mutant.K, K.mutant.K,
                          N.resident, 1)
alpha2.K <- compute_alpha(w.mutant.K, r.mutant.K, K.mutant.K,
                          N.resident, 2)

alpha1.C <- compute_alpha(w.mutant.C, r.mutant.C, K.mutant.C,
                          N.resident, 1)
alpha2.C <- compute_alpha(w.mutant.C, r.mutant.C, K.mutant.C,
                          N.resident, 2)

## This has an erorr in it: the intraspecific competition value is not
## equal to 1, but it should be.
K.mutant

## So the alpha1 (dashed line) is really weird; it gives values that
## are negative (facilitation) for large K values - super weird!
plot(x.mutant, alpha2.K, type="l", ylim=range(alpha1.K, alpha2.K),
     xlab="Invader trait (K)", ylab="Competitive effect")
lines(x.mutant, alpha1.K, type="l", lty=2)
abline(h=c(0, 1), v=eq$x[[1]], lty=3)

## Check that that's actually OK -- yes the difference is what we'd
## predict.
plot(x.mutant, alpha1.K - alpha2.K, type="l")
points(x.mutant,
       w.mutant.K * (K.mutant.K - 1) / (N.resident * r.mutant.K))

## Again, with the 
plot(x.mutant, alpha2.C, type="l",
     ylim=range(alpha1.C, alpha2.C), xlab="Invader trait (C)",
     ylab="Competitive effect")
lines(x.mutant, alpha1.C, type="l", lty=2)
abline(h=c(0, 1), v=eq$x[[1]], lty=3)

## And we can make a little 3d picture of what the competition kernel
## looks like:
ww <- matrix(m_r2$fitness(xx, eq$x, eq$y, eq$R), length(x.mutant))
rr <- m_r2$max_growth_rate(xx)
kk <- m_r2$carrying_capacity(xx)
a1 <- compute_alpha(ww, rr, kk, N.resident, 1)
a2 <- compute_alpha(ww, rr, kk, N.resident, 2)
a1[!is.finite(a1)] <- a2[!is.finite(a2)] <- NA

## Fitness varies with K but not with C
persp(x.mutant, x.mutant, ww, xlab="K", ylab="C", zlab="w",
      theta=30, shade=.2, border=NA, col="dodgerblue4")

## Well, that does not look Gaussian:
persp(x.mutant, x.mutant, a1, xlab="K", ylab="C", zlab="alpha (1)",
      theta=30, shade=.2, border=NA, col="dodgerblue4")
persp(x.mutant, x.mutant, a2, xlab="K", ylab="C", zlab="alpha (2)",
      theta=30, shade=.2, border=NA, col="dodgerblue4")

image(x.mutant, x.mutant, a1, xlab="K", ylab="C")
image(x.mutant, x.mutant, a2, xlab="K", ylab="C")

## It's probably worth overlaying this with the different adaptive
## regions -- species that can invade, coexist and outcompete.

## In this section here, invasion is possible where |K - 0.5| is less
## than the resident value.  That can be generalised to different S
## values, and to different slopes of the K1/K2 tradeoff (at least it
## looks easy graphically).
abline(v=0.5 + c(-1, 1) * abs(eq$x[[1]] - 0.5))
abline(h=eq$x[[2]], lty=2)
slope <- x.mutant / (1 - x.mutant)

image(x.mutant, x.mutant, a2, xlab="K", ylab="C")
