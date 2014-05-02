source("common.R")

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
## demographic equilibrium:
sys_d <- m_d$equilibrium(sys_split(sys_d0, 0.5))

## Here is the fitness landscape
##+ dd_fitness_landscape
x_mutant <- seq(-2, 2, length=101)
plot(x_mutant, m_d$fitness(x_mutant, sys_d$x, sys_d$y), type="l",
     xlab="Trait", ylab="Fitness")
abline(v=sys_d$x, col=1:2, lty=1:2)

## The resident competitive effects are explicitly given in this
## model; they are:
##+ dd_competition_true
matplot(x_mutant, m_d$competition(sys_d$x, x_mutant), type="l", lty=1,
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

## This approach can be generalised over the whole mutant trait space,
## which is what the `model_jacobian_density` function, defined in
## `common.R` does:
model_jacobian_density

## This computes the derivative over a range of mutants (`x_mutant`)
## with respect to the density of both resident types.
z_d <- model_jacobian_density(x_mutant, sys_d, m_d)

## Plotted, showing the resident types (vertical dotted lines) and the
## mutant evaluated in the previous calculations (at `r mutant`) as
## solid points).
##+ dd_derivative
matplot(x_mutant, z_d, type="l", lty=1,
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
capacity_mutant <- m_d$capacity(x_mutant)
matplot(x_mutant, z_d * capacity_mutant, type="l", lty=1,
        xlab="Trait", ylab="Fitness derivative")
abline(v=sys_d$x, lty=3, col=1:2)
points(rep(mutant, 2), z * m_d$capacity(mutant), pch=19, col=1:2)

## Now the competitive effects are exactly centred around the resident
## densities.

## Carrying capacity in an empty environment is a gaussian centred on
## zero:

##+ dd_carrying_capacity
capacity_mutant <- m_d$capacity(x_mutant)
matplot(x_mutant, capacity_mutant, type="l", lty=1,
        xlab="Trait", ylab="Carrying capacity (empty environment)")
abline(v=sys_d$x, lty=3, col=1:2)

## The rate of decline away from zero, relative to the widths of
## competition, allows the system to maintain two species.  However,
## they will be constantly driven apart by disruptive selection.

## Note that in this model this is *different* to the **fitness** in
## an empty environment, which we'll use below.  The fitness in an
## empty environment is 1 for all traits:
m_d$fitness(x_mutant, numeric(0), numeric(0))

