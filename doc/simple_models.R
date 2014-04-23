## # Simple models and competition

##+ echo=FALSE
knit_hooks$set(small_mar=function(before, options, envir) {
  if (before) par(mar=c(4, 4, .1, .1)) # smaller margin on top and right
})
opts_chunk$set(tidy=FALSE, fig.height=5, error=FALSE, small_mar=TRUE,
               fig.path="figure/simple_models_")

## The idea here is to see if we can infer the shape of competition
## from simple models.  I'm trying this on two types of models:
##
## * Models with *explicit* competition, where the competition
##   "kernel" is given as an ingredient of the model and is therefore
##   known.
## * Models with *implicit* competition, where competition emerges
##   from the model, for example through the shared use of a
##   resource.

## There are a number of live issues that I'm unresolved about how to
## deal with:
##
## 1. Are we interested in measuring competition *exerted* (species X
##    depresses how much species Y are capable of growing) or in
##    competition *felt* (species X's growth is depressed by species
##    Y).  In both is that X is the "focal" species and exists at
##    nontrivial density, while Y is something that we want to
##    evaluate at any point along a trait axis.
## 2. Indirect competition and second order effects.  As a result
##    mostly considering single resident cases at the moment.

## I'm skeptical that we can come up with one "true" measure of
## competition.  What I am hoping to be able to do is to describe
## competition *qualitatively*
##
## * is it symmetric or asymmetric?
## * is it linear with respect to population density?  sublinear or
##   superlinear?

## All the models are in the Revolve package:
library(Revolve)

## And we'll need tools in numDeriv
library(numDeriv)

## The basic approach is centred around the idea that we should be
## able to get a handle on competition by looking to see how
## performance of one species varies with the density of another.  So
## if we increase a resident species' density slightly, then if mutant
## fitness decreases a lot then the resident has a strong competitive
## effect on the mutant (so this is competition *exerted* by the
## mutant)

## # Explicit competition

## There are two models here: Dieckmann and Dobeli 1999 and Kisdi
## 1999.  Both these models are related to Lokta Volterra models with
## continuous competition functions that vary as a function of species
## traits rather than being particular to a species.  The species
## competitive effect scales linearly with density.

## ## Dieckmann and Dobeli 1999

### TODO: Probably add in the model definitions in here so that we can
### at least see the general form.

## This makes a copy of the model using the same parameters that the
## paper uses; in particular, these are parameters that allow
## coexistance of multiple species.
m_d <- make_dieckmann_1999()

## There is a single equilibrium at x=0 and population density at the
## carrying capacity.
sys_d0 <- m_d$single_equilibrium()

## Displace this equilibrium by 0.5 to each side and run out to
## demographic eqilibrium:
sys_d <- m_d$equilibrium(sys_split(sys_d0, 0.5))

## Here is the fitness landscape
##+ dd_fitness_landscape
x.mutant <- seq(-2, 2, length=101)
plot(x.mutant, m_d$fitness(x.mutant, sys_d$x, sys_d$y), type="l",
     xlab="Trait", ylab="Fitness")
abline(v=sys_d$x, col=1:2, lty=1:2)

## The resident competitive effects are explicitly given in this
## model; they are:
##+ dd_competition_true
matplot(x.mutant, m_d$competition(sys_d$x, x.mutant), type="l", lty=1,
        xlab="Trait", ylab="Competitive effect")
abline(v=sys_d$x, col=1:2, lty=1:2)

## The idea is that we can compute the derivative of this fitness
## landscape with respect to resident density.

## First, at a single point.  Consider a mutant at position x=0.2; it
## has fitness:
mutant <- 0.2
m_d$fitness(mutant, sys_d$x, sys_d$y)

## As we vary the density of resident individuals and see how fitness
## changes.  Varying just the density of the *first* species we have:
##+ dd_mutant_1
dy <- seq(-10, 10, length=101)
f1 <- function(dy) {
  m_d$fitness(mutant, sys_d$x, sys_d$y + c(dy, 0))
}
plot(dy, sapply(dy, f1), type="l",
     xlab="Change in density of resident", ylab="Fitness")
abline(h=f1(0), lty=3)
abline(v=0, col="grey", lty=1)

## Add the second resident to this:
##+ dd_mutant_2
f2 <- function(dy) {
  m_d$fitness(mutant, sys_d$x, sys_d$y + c(0, dy))
}
plot(dy, sapply(dy, f1), type="l",
     xlab="Change in density of resident", ylab="Fitness")
lines(dy, sapply(dy, f2), type="l", col=2, lty=2)
abline(h=f1(0), lty=3)
abline(v=0, col="grey", lty=1)
legend("topright",
       c("Species 1 (further trait)", "Species 2 (closer trait)"),
       lty=1:2, col=1:2, bty="n")

## The slope of the red line is steeper indicating that the red
## species has a larger competitive effect on the mutant (a unit
## increase of density of the red species decreases fitness faster
## than a unit increase of density in the black species).  The mutant
## is closer to the red species so that makes sense.

## The slopes can be directly computed by using the jacobian function
## in numDeriv:
z <- jacobian(function(y) m_d$fitness(mutant, sys_d$x, y), sys_d$y)
z

## This is the slope of the lines in the previous plot:
##+ dd_mutant_2_slope
plot(dy, sapply(dy, f1), type="l",
     xlab="Change in density of resident", ylab="Fitness")
lines(dy, sapply(dy, f2), type="l", col=2, lty=2)
abline(h=f1(0), lty=3)
abline(v=0, col="grey", lty=1)
legend("topright",
       c("Species 1 (further trait)", "Species 2 (closer trait)"),
       lty=1:2, col=1:2, bty="n")
abcline(0, f1(0), z[[1]], lwd=10, col=make_transparent(1, .2))
abcline(0, f2(0), z[[2]], lwd=10, col=make_transparent(2, .2))

## This approach can be generalised over the whole mutant trait space:
model_jacobian_density <- function(x, sys, m) {
  jacobian(function(y) m$fitness(x, sys$x, y), sys$y)
}

## This computes the derivative over a range of mutants (`x.mutant`)
## with respect to the density of both resident types.
z_d <- model_jacobian_density(x.mutant, sys_d, m_d)

## Plotted, showing the resident types (vertical dotted lines) and the
## mutant evaluated in the previous calculations (at `r mutant`) as
## solid points).
##+ dd_derivative
matplot(x.mutant, z_d, type="l", lty=1,
        xlab="Trait", ylab="Fitness derivative")
abline(v=sys_d$x, lty=3, col=1:2)
points(rep(mutant, 2), z, pch=19, col=1:2)

## This seems to capture most of the important bits about competition:
## it's roughly symmetric around the resident, with the strongest
## impact felt at the same trait value as the resident.  Note that in
##
## However, the competition functions are offset away from the
## resident locations slightly.  If we rescale the derivatives by
## multiplying through by the carrying capacity *in an empty
## environment* that offset goes away:
##+ dd_derivative_scaled
capacity.mutant <- m_d$capacity(x.mutant)
matplot(x.mutant, z_d * capacity.mutant, type="l", lty=1,
        xlab="Trait", ylab="Fitness derivative")
abline(v=sys_d$x, lty=3, col=1:2)
points(rep(mutant, 2), z * m_d$capacity(mutant), pch=19, col=1:2)

## Now the competitive effects are exactly centred around the resident
## densities.

## Note that in this model this is *different* to the fitness in an
## empty environment, which we'll use below.  The fitness in an empty
## environment is 1 for all traits:
m_d$fitness(x.mutant, numeric(0), numeric(0))

## # Appendix

## Package version information from `sessionInfo()`:
##+ echo=FALSE
sessionInfo()
