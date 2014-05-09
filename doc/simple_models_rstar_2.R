source("common.R")

## # Rstar model, two resource

## This version of the model is based loosely off Huisman and Weissing
## 2001.  With two resources, things are significantly more
## complicated than the [single resource case](simple_models_rstar_1).
## Some issues (mostly unresolved at the moment)
##
## * A species requires four parameters ($K_1$, $K_2$, $C_1, $C_2$),
##   but should these all be free?  We can imagine a trade-off where
##   $K_1 = 1 - K_2$ (as done in some R star papers).
## * Generalisation of the resource use relationships (e.g., following
##   Schreiber and Tobiason
##   [2003](http://doi.org/10.1007/s00285-003-0195-9), and the
##   treatment in Tilman's book.
##
## This case is also more difficult to think about because it's harder
## to boil down the entire intuition of the model into a few sentences
## like we can with the 1 dimensional model (partly because there is
## not one possible model).

## First consider what happens when we have species that vary in their
## preference for resources, but have identical uptake levels (uptake
## levels won't affect mutant fitness so these can be ignored for
## now).  These, and the resident C values will be set to 0.5 for both
## resources for now.
mat_r2 <- rstar_matrices(rstar_mat_2_tradeoff,
                         make_rstar_mat_constant(matrix(0.5, 2)))
m_r2 <- make_rstar(mat_r2, S=c(0.5, 0.5))

sys_r2 <- sys(matrix(0.3), 1)
eq <- m_r2$single_equilibrium(sys_r2$x)

## Here is the ZNGI plot, with the resource consumption vector
## shown as a dashed diagonal line.
##+ r2_zngi_1sp
rstar_plot(m_r2, sys_r2)

## Ultimately, in this case, long term coexistence depends on
## consumption vectors. With every species having the same consumption
## vector, as here, there will be only one possible resident species
## and the outcome will depend on which side of the diagonal line we
## start from (the supply point, S).

## This formulation of the model boils down to a single parameter: how
## much of resource 1 is needed for growth relative to resource 2.  So
## we can compute fitness over this range:
x_mutant <- matrix(seq(0, 1, length=301), 1)
w_mutant <- m_r2$fitness(x_mutant, eq$x, eq$y)
##+ r2_fitness
plot(x_mutant, w_mutant, type="l", xlab="Trait (K1)", ylab="Fitness")
abline(h=0, lty=3)
abline(v=eq$x, lty=2)

## Fitness has a peak at 0.5 -- the point at which both resources are
## required equally.  A species with a trait value closer to 0.5 than
## the resident can invade -- fitness is positive there.

## The fitness landscape in the empty environment looks similar,
## though it is positive for all trait values.
## (TODO: retrieving S should not be needed here -- bug in Revolve)
get_S <- function(m) {
  m$parameters$get()[["S"]]
}
w_empty <- m_r2$fitness(x_mutant, eq$x, 0, get_S(m_r2))
##+ r2_fitness_empty
plot(x_mutant, w_empty, type="l", xlab="Trait (K1)", ylab="Fitness")
abline(v=eq$x, lty=2)

## Then compute the derivative of fitness with respect to resident
## density.
z_r2 <- model_jacobian_density(x_mutant, eq, m_r2)

## As usual, this shows the scaling issues -- trait values that have
## relatively high fitnesses show the greatest derivatives.
##+ r2_derivative
plot(x_mutant, z_r2, type="l", xlab="Mutant trait (K1) value",
     ylab="Fitness derivative")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

## Rescale the derivative by fitness in the empty environment (OK, I
## think that this is the same as looking at the derivative of
## w/w_empty, which will be well behaved everywhere that w_empty is).
##+ r2_derivative_scaled
plot(x_mutant, z_r2 / w_empty, type="l",
     xlab="Mutant trait (K1) value", ylab="Fitness derivative, scaled")
abline(v=eq$x[1], lty=2)
abline(h=0, col="grey", lty=3)

## This tells the opposite story -- the resident exerts the strongest
## proportional competitive effect on invading strategies with K1
## values *further* from 0.5 than itself.
