## # Simple models and competition

### TODO:
### * Subtraction approach for the simple models.

##+ echo=FALSE
knitr::knit_hooks$set(small_mar=function(before, options, envir) {
  if (before) par(mar=c(4, 4, .1, .1)) # smaller margin on top and right
})
knitr::opts_chunk$set(tidy=FALSE, fig.height=5, error=FALSE,
                      small_mar=TRUE,
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
abline(h=0, lty=3)
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

## ## Kisdi 1999

### TODO: Probably add in the model definitions in here so that we can
### at least see the general form.
m_k <- make_kisdi_1999()
sys_k0 <- m_k$single_equilibrium()
sys_k <- m_k$equilibrium(sys_split(sys_k0, 0.1))

## Here is the fitness landscape
##+ k_fitness_landscape
x.mutant <- seq(-0.5, 1.5, length=301)
plot(x.mutant, m_k$fitness(x.mutant, sys_k$x, sys_k$y), type="l",
     xlab="Trait", ylab="Fitness")
abline(v=sys_k$x, col=1:2, lty=1:2)

## It's not totally obvious here, but there is some fluctuations
## around the local traits:
##+ k_fitness_landscape_detail
plot(x.mutant, m_k$fitness(x.mutant, sys_k$x, sys_k$y), type="l",
     xlab="Trait", ylab="Fitness",
     xlim=range(sys_k$x) + c(-.05, .05), ylim=c(-0.005, 0.005))
abline(v=sys_k$x, col=1:2, lty=1:2)
abline(h=0, lty=3)

## Here is the true competitive effect.

## This is how much competition is felt by mutants at different
## mutant trait values
##+ k_competition_felt
matplot(x.mutant, t(m_k$competition(x.mutant, sys_k$x)), type="l", lty=1,
        xlab="Trait", ylab="Competitive effect")
abline(v=sys_k$x, col=1:2, lty=1:2)

## In the D+D model, we didn't really have to specify much about the
## species exerting and recieving competition, but here we do.  This
## is the amount that a mutant whose trait is at a position on the x
## axis has its growth reduced by the black and the red species; so
## this is competition *felt*.

## In contrast, this is how much competition is *exerted* by the
## mutants on the black and the red species:
##+ k_competition_true
matplot(x.mutant, m_k$competition(sys_k$x, x.mutant), type="l", lty=1,
        xlab="Trait", ylab="Competitive effect")
abline(v=sys_k$x, col=1:2, lty=1:2)

## We can jump straight to the derivatives:
##+ k_derivative
z_k <- model_jacobian_density(x.mutant, sys_k, m_k)
matplot(x.mutant, z_k, type="l", lty=1,
        xlab="Trait", ylab="Fitness derivative")
abline(v=sys_k$x, lty=3, col=1:2)
abline(h=0, lty=3)

## Here, a competitive effect of zero means that the species is
## unaffected by changing density of the resident.  So a mutant
## species around 1.5 is basically unaffected by the black species,
## while a species at -0.5 has its growth strongly suppressed by both
## species.  This is consistent with the competitive hierarchy.

## In this model, the intrinsic growth rate is the fitness in an empty
## environment:
##+ k_intrinsic_growth
plot(x.mutant, m_k$intrinsic_growth(x.mutant), type="l",
     xlab="Trait", ylab="Intrinsic growth rate")
abline(h=0, lty=3)

##+ k_fitness_empty
plot(x.mutant, m_k$fitness(x.mutant, numeric(0), numeric(0)), type="l",
     xlab="Trait", ylab="Growth in empty environment")
abline(h=0, lty=3)

## We can scale the derivatives by this value, but nothing meaningful
## emerges.
##+ k_derivative_scaled_1
intrinsic.mutant <- m_k$intrinsic_growth(x.mutant)
matplot(x.mutant, z_k * intrinsic.mutant, type="l", lty=1,
        xlab="Trait", ylab="Scaled fitness derivative")
abline(v=sys_k$x, lty=3, col=1:2)

##+ k_derivative_scaled_2
matplot(x.mutant, z_k / intrinsic.mutant, type="l", lty=1,
        xlab="Trait", ylab="Scaled fitness derivative")
abline(v=sys_k$x, lty=3, col=1:2)

## # Implicit competition

## This is where things get more interesting.

## ## Rstar model

## This version of the model is based loosely off Huisman and Weissing
## 2001.

## ### Single resource

## I'm going to consider the single resource case first because it's a
## bunch easier to think about than the two (or more!) resource case.
## It's also confusing enough to be worth thinking about.
##
## In this model, there can only be a single dominant species, so we
## can't study what happens during coexistance.  However, there's
## still plenty of competition going on.

## In the single resource case species are characterised by two
## traits: K (the constant of half-saturation of the Monod curve --
## higher K indicates higher requirements for the resource) and C (the
## rate of consumption of the resource).  All else being equal, the
## species with the lowest K will win because they will drive K below
## the minimum level required for positive net growth by the other
## species.
m_r1 <- make_rstar(rstar_matrices(rstar_mat_1, rstar_mat_1), S=1)
sys_r1 <- list(x=matrix(0.5, nrow=2), y=1)

eq <- m_r1$single_equilibrium(sys_r1$x)

## Look at the fitness landscape: how does the instantaneous growth
## rate look with respect to mutant C and K.
x.mutant <- seq(0, 1, length=101)
x.K <- rbind(x.mutant, eq$x[2], deparse.level=0) # K varying

## Here is the fitness landscape with respect to mutant K (changing C
## doesn't change the fitness landscape because it only affects
## species at nontrivial densities).
##+ r1_fitness_landscape
plot(x.mutant, m_r1$fitness(x.K, eq$x, eq$y, eq$R), type="l",
     xlab="Trait (K)", ylab="Fitness")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

## As with the Dieckmann and Dobeli model, vary the density of the
## resident a little and see what the response of the fitness
## landscape is.

x.mutant1 <- rbind(0.3, eq$x[[2]])
y.resident <- eq$y + seq(-1, 1, length=101)
w.mutant1 <- sapply(y.resident, function(y)
                    m_r1$fitness(x.mutant1, eq$x, y))
w.mutant0 <- m_r1$fitness(x.mutant1, eq$x, eq$y) # at equil. density

##+ r1_mutant_1
plot(w.mutant1 ~ y.resident, type="l",
     xlab="Resident density", ylab="Fitness")
points(eq$y, w.mutant0, pch=19)

## Then compute the slope
##+ r1_mutant_1_slope
plot(w.mutant1 ~ y.resident, type="l",
     xlab="Resident density", ylab="Fitness")
points(eq$y, w.mutant0, pch=19)
z.mutant1 <- jacobian(function(y) m_r1$fitness(x.mutant1, eq$x, y), eq$y)
abcline(eq$y, w.mutant0, z.mutant1,
        lwd=10, col=make_transparent(1, .2))

## First, note that in contrast with the D+D model, fitness is
## nonlinearly related to resident density (i.e., the line has a
## nonzero second derivative).  As resident density increases, mutant
## fitness decreases, but this happens *slower than linearly*.

## Computing this slope over the entire range of mutant trait values:
z <- jacobian(function(y) m_r1$fitness(x.K, eq$x, y), eq$y)

##+ r1_derivative
plot(x.mutant, z, type="l", xlab="Mutant trait (K) value",
     ylab="Fitness derivative")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

## As above, the more negative this is, the more a species growth rate
## is depressed by a unit increase in redsident density.  The resident
## is indicated by the vertical dashed line.  In contrast with the
## explicit competition models, the line is neither centred on the
## resident nor monotonic.

## The greatest sensitivity to the density of the resident species is
## around 0.17.  Counter-intuitively though, this species can
## competitively displace the resident (invasion is always possible
## for smaller values of K).

## However, species with a smaller K value also have a higher
## "intrinsic" growth rate -- a growth rate in a totally empty
## environment.
##+ r1_fitness_empty
w.empty <- m_r1$fitness(x.K, eq$x, 0, m_r1$parameters$get()[["S"]])
plot(x.mutant, w.empty, type="l", xlab="Mutant trait (K) value",
     ylab="Fitness in an empty environment")

## If we scale the fitness derivative by this maximum potential
## fitness we get the fractional reduction in fitness by a unit
## increase in density of the resident species:
##+ r1_derivative_scaled
plot(x.mutant, z / w.empty, type="l",
     xlab="Mutant trait (K) value", ylab="Fitness derivative, scaled")
abline(v=eq$x[1], lty=2)
abline(h=0, col="grey", lty=3)

## This function is monotonic and shows the greatest reduction in
## fitness for species with greater K values -- the same species that
## would be competitively replaced.  The function is flattest around
## the resident species (vertical dashed line).  And the impact on
## species that have very low K values is very small.  This seems to
## agree with the intuitive behaviour of the model.

## # Unresolved things

## * I've scaled by different things all over the show, until things
##   look "about right" -- would be nice to know what was going on
##   here.  I suspect that working out what units things are in will
##   enlightening.



## # Appendix

## Package version information from `sessionInfo()`:
##+ echo=FALSE
sessionInfo()
