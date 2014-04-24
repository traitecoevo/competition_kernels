library(Revolve)
library(numDeriv)

# First consider what happens when we have species that vary in their
# preference for resources, but have identical uptake levels.
mat <- rstar_matrices(rstar_mat_2_tradeoff,
                        make_rstar_mat_constant(matrix(0.5, 2)))
m <- make_rstar(mat, S=c(1,1))

sys0 <- sys(matrix(0.3), 1)

# Here is the ZNGI plot:
rstar_plot(m, sys0)

# Ultimately, in this case the evolutionary race will be determined by
# which side of the line we start on.  Changing C would allow
# coexistence, but does not actually figure in here!

# So, let's proceed as if this was not the case.
t <- seq(0, 30, length=201)

res <- m$run(sys0, t)
matplot(res$t, cbind(res$R, res$y), type="l", lty=c(2,2,1), xlab="Time",
        ylab="Abundance (solid), Resource (dashed)")

eq <- m$single_equilibrium(sys0$x)
abline(h=unlist(eq[c("R", "y")]), col=c(1:3), lty=3)

# Displace the solution from equilibrium and look at the new level of
# resources:
dy <- eq$y * 0.1
sys.y1 <- modifyList(eq, list(y=eq$y - dy, R=eq$R))
sys.y2 <- modifyList(eq, list(y=eq$y + dy, R=eq$R))

res.y1 <- m$run_fixed_density(sys.y1, t)
res.y2 <- m$run_fixed_density(sys.y2, t)

eq.y1 <- m$equilibrium_R(sys.y1)
eq.y2 <- m$equilibrium_R(sys.y2)

matplot(res.y1$t, cbind(res.y1$R, res.y2$R), type="l",
        col=c(2,2,3,3))
abline(h=eq$R, lty=3, col=1)
abline(h=eq.y1$R, lty=3, col=2)
abline(h=eq.y2$R, lty=3, col=3)
# Determined analytically:
points(rep(max(res.y1$t), 2), m$single_equilibrium_R(sys.y1), col=2)
points(rep(max(res.y1$t), 2), m$single_equilibrium_R(sys.y2), col=3)

# Next, we start working towards the instantaneous growth rate of a
# new type at this equilibrium

# This computes all the same things as in the 1 resource case, given a
# model, species traits, and mutant x values:
f <- function(m, x.res, x.mutant, pars=list()) {
  op <- m$parameters$set(pars)
  on.exit(m$parameters$set(op))
  eq <- m$single_equilibrium(x.res)
  S <- m$parameters$get()[["S"]]
  w.empty  <- m$fitness(x.mutant, x.res, 0, S)
  w.mutant <- m$fitness(x.mutant, x.res, eq$y)
  z <- jacobian(function(y) m$fitness(x.mutant, x.res, y), eq$y)
  list(x.mutant=x.mutant, w.empty=w.empty, w.mutant=w.mutant, z=z,
       x.res=x.res)
}

# And this plots them:
plt <- function(obj) {
  op <- par(mfrow=c(5, 1), mar=c(.5, 4.1, .5, .5), oma=c(3, 0, 0, 0))
  on.exit(par(op))
  plot(w.empty ~ drop(x.mutant), obj, type="l", xaxt="n")
  abline(v=obj$x.res, lty=2)
  abline(h=0, col="grey", lty=3)

  plot(w.mutant ~ drop(x.mutant), obj, type="l", xaxt="n")
  abline(v=obj$x.res, lty=2)
  abline(h=0, col="grey", lty=3)

  plot(w.mutant - w.empty ~ drop(x.mutant), obj, type="l", xaxt="n")
  abline(h=0, col="grey", lty=3)
  abline(v=obj$x.res[1], lty=2)

  plot(z ~ drop(x.mutant), obj, type="l", xaxt="n")
  abline(h=0, col="grey", lty=3)
  abline(v=obj$x.res[1], lty=2)

  plot(z / w.empty ~ drop(x.mutant), obj, type="l")
  abline(v=obj$x.res[1], lty=2)
}

xm <- matrix(seq(0, 1, length=301), nrow=1)

# Changing the supply vector:
plt(f(m, sys0$x, xm))
plt(f(m, sys0$x, xm, pars=list(S=c(.8, 1))))
plt(f(m, sys0$x, xm, pars=list(S=c(.5, 1))))
plt(f(m, sys0$x, xm, pars=list(S=c(.3, 1))))

# Changing the resident K values
plt(f(m, matrix(.2), xm, pars=list(S=c(.5, 1))))
plt(f(m, matrix(.1), xm, pars=list(S=c(.5, 1))))
plt(f(m, matrix(.05), xm, pars=list(S=c(.5, 1))))
plt(f(m, matrix(.025), xm, pars=list(S=c(.5, 1))))
plt(f(m, matrix(.0125), xm, pars=list(S=c(.5, 1))))

# Changing the resident C values
g <- function(C, K, S, xm) {
  mat <- rstar_matrices(rstar_mat_2_tradeoff,
                        make_rstar_mat_constant(matrix(C, 2)))
  sys0$x <- matrix(K)
  m <- make_rstar(mat, S=c(1,1))
  f(m, sys0$x, xm, pars=list(S=S))
}

plt(g(c(0.5, 0.5), 0.3, c(1, 1), xm))
plt(g(c(0.6, 0.4), 0.3, c(1, 1), xm))
plt(g(c(0.7, 0.3), 0.3, c(1, 1), xm))
plt(g(c(0.8, 0.2), 0.3, c(1, 1), xm))
plt(g(c(0.9, 0.1), 0.3, c(1, 1), xm))
