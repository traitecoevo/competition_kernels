library(Revolve)
library(deSolve)
library(numDeriv)

# This should move into the model:
singular_kisdi <- function(m) {
  p <- m$parameters$get()
  alpha <- function(xp) drop(m$competition(0, 0 - xp))
  xp <- p$beta / p$b + alpha(0) / grad(alpha, 0)
  list(x=xp, y=m$capacity(xp) / alpha(0))
}

m_d <- make_dieckmann_1999()
m_k <- make_kisdi_1999()

step_d <- make_step_equilibrium(m_d$fitness, method="runsteady")
step_k <- make_step_equilibrium(m_k$fitness, method="nleqslv")

## The D+D model, stepped to demographic equilibrium
sys_d <- step_d(list(x=c(-1, 1) * 0.5, y=rep(1, 2), t=0))
sys_k <- with(singular_kisdi(m_k),
              step_k(list(x=x + c(-1, 1) * 0.1, y=rep(y / 2, 2), t=0)))

# Trait space to look over.  This range is actually useful for both
# models.
xx <- seq(-2, 2, length=101)
plot(xx, m_d$fitness(xx, sys_d$x, sys_d$y), type="l",
     xlab="Trait", ylab="Fitness")
abline(v=sys_d$x, col=1:2)

# We're interested in changes in fitness as the densities of residents
# changes.
target <- function(y, x, sys, m) {
  m$fitness(x, sys$x, y)
}

z_d <- jacobian(function(y) target(y, xx, sys_d, m_d), sys_d$y)

matplot(xx, z_d, type="l", lty=1, xlab="Trait", ylab="Fitness derivative")
abline(v=sys_d$x, lty=3, col=1:2)

scal <- - m_d$parameters$get()$r / m_d$capacity(xx)
c_real_d <- m_d$competition(sys_d$x, xx)

matplot(xx, c_real_d, type="l", lty=1, xlab="Trait", ylab="True competition")
abline(v=sys_d$x, lty=3, col=1:2)

# This is why we're close:
all.equal(z_d / scal, c_real_d)

xx <- seq(-2, 2, length=101)
plot(xx, m_k$fitness(xx, sys_k$x, sys_k$y), type="l",
     xlab="Trait", ylab="Fitness")
abline(v=sys_k$x, col=1:2)

c_real_k <- m_k$competition(sys_k$x, xx)
matplot(xx, c_real_k, type="l", lty=1, xlab="Trait", ylab="True competition")
abline(v=sys_k$x, lty=3, col=1:2)

z_k <- jacobian(function(y) target(y, xx, sys_k, m_k), sys_k$y)
matplot(xx, z_k, type="l", lty=1, xlab="Trait", ylab="Fitness derivative")
abline(v=sys_k$x, lty=3, col=1:2)

plot(xx, m_k$capacity(xx), type="l",
     xlab="Trait", ylab="Carrying capacity")

m_k <- make_kisdi_1999(capacity="gaussian")
sys_k <- with(singular_kisdi(m_k),
              step_k(list(x=x + c(-1, 1) * 0.1, y=rep(y / 2, 2), t=0)))






i <- 1
x <- 0
g <- function(y) f(y, i, sys, x, m)
z <- grad(g, sys$y[i])

h <- function(y) f(y, i, sys, xx, m)
zz <- drop(jacobian(h, sys$y[i]))

plot(xx, zz, type="l")
zz[which.min(abs(xx - x))] - z

# I think that we can actually do better than this; we can analyse all
# species at the same time!
i <- 1:2
h <- function(y) f(y, i, sys, xx, m)
zz <- jacobian(h, sys$y[i])


## These aren't *quite* where we'd expect them.
pdf("competition_dd.pdf")
matplot(xx, zz, type="l", lty=1)
abline(v=sys$x, lty=3, col=1:2)
dev.off()


## One species
sys <- step(list(x=c(-1, 1) * 0.5, y=c(1, 0), t=0))

sys$y[2] <- 0
h <- function(y) f(y, i, sys, xx, m)
zz <- drop(jacobian(h, sys$y[i]))
plot(xx, zz, type="l")

## There's an offset here.  PRobably due to both indirect competition
## and the carrying capacity.  Subtracting off the carrying capacity
## effect would be useful here.

## Then by simple finite differencing; looks like we agree pretty well
## now.
y0 <- sys$y[i]
dy <- y0 * 0.001
(g(y0 + dy) - g(y0 - dy)) / (2 * dy) - z

## So, now for the kisdi model:
k <- make_kisdi_1999()
singular_kisdi <- function(m) {
  p <- m$parameters$get()
  alpha <- function(xp) drop(m$competition(0, 0 - xp))
  xp <- p$beta / p$b + alpha(0) / grad(alpha, 0)
  list(x=xp, y=m$capacity(xp) / alpha(0))
}
step_k <- make_step_equilibrium(k$fitness, method="nleqslv")
sys_k <- with(singular_kisdi(k),
              step_k(list(x=x + c(-1, 1) * 0.1, y=rep(y / 2, 2), t=0)))

i <- 1:2
xx <- seq(-2, 2, length=101)
h <- function(y) f(y, i, sys_k, xx, k)
zz <- jacobian(h, sys_k$y[i])

pdf("competition_kisdi.pdf")
matplot(xx, zz, type="l", lty=1)
abline(v=sys_k$x, lty=3, col=1:2)
dev.off()


