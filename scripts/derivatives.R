library(Revolve)
library(deSolve)
library(numDeriv)

# Trait space to look over.  This range is actually useful for both
# models.
xx <- seq(-2, 2, length=101)

# We're interested in changes in fitness as the densities of residents
# changes.  This little wrapper function will help for all the models
# that are implemented in Revolve.
model_jacobian_density <- function(x, sys, m) {
  jacobian(function(y) m$fitness(x, sys$x, y), sys$y)
}

# # Dieckmann and Dobeli 1999:

m_d <- make_dieckmann_1999()
sys_d <- m_d$equilibrium(x=c(-1, 1) * 0.5, y=rep(1, 2))

# Fitness landscape for the D+D model:
plot(xx, m_d$fitness(xx, sys_d$x, sys_d$y), type="l",
     xlab="Trait", ylab="Fitness")
abline(v=sys_d$x, col=1:2)

z_d <- model_jacobian_density(xx, sys_d, m_d)

matplot(xx, z_d, type="l", lty=1, xlab="Trait", ylab="Fitness derivative")
abline(v=sys_d$x, lty=3, col=1:2)

scal <- - m_d$parameters$get()$r / m_d$capacity(xx)
c_true_d <- m_d$competition(sys_d$x, xx)

matplot(xx, c_true_d, type="l", lty=1, xlab="Trait", ylab="True competition")
abline(v=sys_d$x, lty=3, col=1:2)

# This is why we're close:
all.equal(z_d / scal, c_true_d)

# # Kisdi 1999

m_k <- make_kisdi_1999()
step_k <- make_step_equilibrium(m_k$fitness, method="nleqslv")
sys_k0 <- m_k$single_equilibrium()
sys_k0$x <- sys_k0$x + c(-1, 1) * 0.1
sys_k0$y <- rep(sys_k0$y/2, 2)
sys_k <- step_k(sys_k0)

xx <- seq(-2, 2, length=101)
plot(xx, m_k$fitness(xx, sys_k$x, sys_k$y), type="l",
     xlab="Trait", ylab="Fitness")
abline(v=sys_k$x, col=1:2)

z_k <- model_jacobian_density(xx, sys_k, m_k)
matplot(xx, z_k, type="l", lty=1, xlab="Trait", ylab="Fitness derivative")
abline(v=sys_k$x, lty=3, col=1:2)

c_true_k <- m_k$competition(sys_k$x, xx)
matplot(xx, c_true_k, type="l", lty=1, xlab="Trait", ylab="True competition")
abline(v=sys_k$x, lty=3, col=1:2)

plot(xx, m_k$intrinsic_growth(xx), type="l",
     xlab="Trait", ylab="Carrying capacity")

# # All together:

par(mfcol=c(3, 2), mar=c(2.1, 2.1, .5, .5), oma=c(2, 2, 2, 2))

plot(xx, m_d$fitness(xx, sys_d$x, sys_d$y), type="l",
     xlab="", ylab="")
abline(v=sys_d$x, col=1:2)
mtext("D+D 1999", 3, line=1)
mtext("Fitness", 2, line=2.5)

matplot(xx, z_d, type="l", lty=1, xlab="", ylab="")
abline(v=sys_d$x, lty=3, col=1:2)
mtext("Fitness derivative", 2, line=2.5)

matplot(xx, c_true_d, type="l", lty=1, xlab="", ylab="")
abline(v=sys_d$x, lty=3, col=1:2)
mtext("True competition", 2, line=2.5)
mtext("Trait", 1, line=2.5)


plot(xx, m_k$fitness(xx, sys_k$x, sys_k$y), type="l",
     xlab="", ylab="")
abline(v=sys_k$x, col=1:2)
mtext("Kisdi 1999", 3, line=1)

matplot(xx, z_k, type="l", lty=1, xlab="", ylab="")
abline(v=sys_k$x, lty=3, col=1:2)

matplot(xx, c_true_k, type="l", lty=1, xlab="", ylab="")
abline(v=sys_k$x, lty=3, col=1:2)
mtext("Trait", 1, line=2.5)
