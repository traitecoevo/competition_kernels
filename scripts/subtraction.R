## It might be best to think about this in the two-species
## equilibrium; then the system should collapse down nicely.
library(numDeriv)
library(Revolve)

k <- make_kisdi_1999()
d <- make_dieckmann_1999()

step_k <- make_step_equilibrium(k$fitness, method="nleqslv")
step_d <- make_step_equilibrium(d$fitness, method="runsteady")

## Identify the singlular points. For the Kisdi model:
singular_kisdi <- function(m) {
  p <- m$parameters$get()
  alpha <- function(xp) drop(m$competition(0, 0 - xp))
  xp <- p$beta / p$b + alpha(0) / grad(alpha, 0)
  list(x=xp, y=m$capacity(xp) / alpha(0))
}

## Set up a two species system that is stable:
sys_k <- with(singular_kisdi(k),
              step_k(list(x=x + c(-1, 1) * 0.1, y=rep(y / 2, 2), t=0)))

## Similarly for the D+D model:
sys_d <- step_d(list(x=c(-1, 1) * 0.5, y=rep(1, 2), t=0))

## Ways of looking at things; subtract away one strategy to see how
## we'd do:
f <- function(i, x, sys, fitness, step) {
  w.with <- fitness(x, sys$x, sys$y)
  w.without <- fitness(x, sys$x[-i], sys$y[-i])
  tmp <- step(list(x=sys$x[-i], y=sys$y[-i], t=0))
  w.without.eq <- fitness(x, tmp$x, tmp$y)
  cbind(with=w.with,
        without=w.without,
        without.eq=w.without.eq)
}

xx_k <- seq(-1, 2, length=101)
yy_k <- f(1, xx_k, sys_k, k$fitness, step_k)
matplot(xx_k, yy_k, type="l")
matplot(xx_k, yy_k - yy_k[,1], type="l", ylim=c(0, 2))

## What does the actual competition function look like there?
lines(xx_k, k$competition(xx_k, sys_k$x[1]), col="blue")
abline(v=sys_k$x, col="grey", lty=3)

xx_d <- seq(-3, 2, length=101)
yy_d <- f(1, xx_d, sys_d, d$fitness, step_d)
matplot(xx_d, yy_d, type="l")
matplot(xx_d, yy_d - yy_d[,1], type="l")

## What does the actual competition function look like there?
lines(xx_d, d$competition(xx_d, sys_d$x[1]), col="blue")
abline(v=sys_d$x, col="grey", lty=3)
