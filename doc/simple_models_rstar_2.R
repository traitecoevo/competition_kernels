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
w_empty <- m_r2$fitness(x_mutant, eq$x, 0, m_r2$parameters$get()[["S"]])
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

## ## Varying supply rates

## Things start getting a bit more complicated once start allowing the
## resource supply rates, S, to change.

## Following Tilman, I'll move the resources along a gradient holding
## S1 + S2 = 1.  What matters for invasion is which side of the C
## vector the starting point lands on, so this should be sufficient.
## Note that this does not include the previous example at this point...

## The resident species will survive wherever the S line falls above
## the ZNGIs.
Rstar <- m_r2$Rstar(sys_r2$x)

c.slope <- 1
c.intercept <- Rstar[2] - c.slope * Rstar[1]
S.crit <- (1 - c.intercept) / (1 + c.slope)

S1 <- sort(c(seq(0, 1, length=51), Rstar[1], 1 - Rstar[2], S.crit))
S <- cbind(S1=S1, S2=1 - S1)

##+ r2_rstar_S
rstar_plot(m_r2, sys_r2)
abline(1, -1, lty=3, col="blue")
segments(Rstar[1], 1 - Rstar[1], 1 - Rstar[2], Rstar[2],
         col="blue")
points(S.crit , 1 - S.crit, col="blue")

m_r2_S <- apply(S, 1, function(S) make_rstar(mat_r2, S=S))
eq_r2_S <- lapply(m_r2_S, function(m) m$single_equilibrium(eq$x))
y_r2_S <- sapply(eq_r2_S, "[[", "y")

## Resident density as a function of the supply rate of the first
## resource (when it directly trades off with the second).  The
## optimum resident density occurs when the supply point is at S.crit
##+ r2_fitness_resident_S
plot(S1, y_r2_S, type="l",
     xlab="Supply rate of resource 1", ylab="Resident density")
abline(v=c(Rstar[1], 1 - Rstar[2]), lty=3)
abline(v=0.5, lty=2)
abline(v=S.crit, lty=2, col="blue")

## So, outside of the critical S values determined by the resident
## Rstar, measures of competition don't make any sense because there
## is no resident species to have a competitive effect on invaders.
## This might cause issues later...
w_r2_S <- mapply(function(m, eq)
                 m$fitness(x_mutant, eq$x, eq$y, eq$R),
                 m_r2_S, eq_r2_S)

red.blue <- colorRampPalette(RColorBrewer::brewer.pal(4, "RdBu"))
cols <- red.blue(length(S1))

## Here is a plot of mutant fitness as a function of K1 value for a
## range of different resource supply ratios.  Red lines have
## relatively high rates of S1, enabling the carrying capacity of the
## resident to be higher.  Blue lines have S2 > S1, and relatively
## lower densities of the resident.
##+ r2_fitness_S
matplot(drop(x_mutant), w_r2_S, type="l", col=cols, lty=1,
        xlab="Mutant trait (K1) value",
        ylab="Fitness")
abline(h=0, lty=3)
abline(v=eq$x[1], lty=2)
lines(x_mutant, w_mutant, lty=2)

## Empty fitness
get_S <- function(m) {
  m$parameters$get()[["S"]]
}
w_r2_empty_S <- mapply(function(m, eq)
                       m$fitness(x_mutant, eq$x, 0, get_S(m)),
                       m_r2_S, eq_r2_S)

##+ r2_fitness_empty_S
matplot(drop(x_mutant), w_r2_empty_S, type="l", col=cols, lty=1,
        xlab="Mutant trait (K1) value",
        ylab="Fitness in empty environment")
abline(h=0, lty=3)
abline(v=eq$x[1], lty=2)
lines(x_mutant, w_empty, lty=2)

ok <- apply(t(S) > drop(Rstar), 2, all)

## These are going to be easier to draw if we code the lines up by
## resident surviving.  Basically where the resident can exist it
## *will* draw the resource down to Rstar, and the fitness function
## will then pass through that point.
##+ r2_fitness_compare_S
op <- par(mfrow=c(2, 1))
matplot(drop(x_mutant), w_r2_S, type="l",
        col=ifelse(ok, cols, make_transparent(cols, .25)),
        lty=ifelse(ok, 1, 2),
        xlab="Mutant trait (K1) value",
        ylab="Fitness")
abline(h=0, lty=3)
abline(v=eq$x[1], lty=2)
lines(x_mutant, w_mutant, lty=2)

matplot(drop(x_mutant), w_r2_empty_S, type="l",
        col=ifelse(ok, cols, make_transparent(cols, .25)),
        lty=ifelse(ok, 1, 2),
        xlab="Mutant trait (K1) value",
        ylab="Fitness in empty environment")
abline(h=0, lty=3)
abline(v=eq$x[1], lty=2)
lines(x_mutant, w_empty, lty=2)
par(op)

## Compute the fitness derivatives with respect to resident density.
## This can only be done easily where the equilibrium resident density
## is at least zero.  We can use the same hack as above for zero
## equilibrium though.

##+ r2_derivative_S_calculate, cache=TRUE
z_r2_S <- mapply(function(m, eq)
                 model_jacobian_density(x_mutant, eq, m),
                 m_r2_S[ok], eq_r2_S[ok])

##+ r2_derivative_S
matplot(drop(x_mutant), z_r2_S, col=cols[ok], type="l", lty=1,
        xlab="Mutant trait (K1) value", ylab="Fitness derivative")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)
lines(x_mutant, model_jacobian_density(x_mutant, eq, m_r2),
      lty=2)

## Then we can rescale this by the fitness in an empty environment.

z_r2_scaled_S <- z_r2_S / w_r2_empty_S[,ok]
z_r2_scaled_twice_S <- t(t(z_r2_scaled_S) * y_r2_S[ok])

## There are some really weird discontinuities here that I need to
## chase up...
##+ r2_derivative_scaled_S
matplot(drop(x_mutant), z_r2_scaled_S, type="l", lty=1,
        col=cols[ok],
        xlab="Mutant trait (K1) value",
        ylab="Fitness derivative, scaled",
        ylim=c(-5,0))
abline(v=eq$x[1], lty=2)
abline(h=0, col="grey", lty=3)
lines(x_mutant, z_r2 / w_empty, lty=2)

## And this really makes no sense at all.  Back to the drawing board
## here, I think.
##+ r2_derivative_scaled_twice_S
matplot(drop(x_mutant), z_r2_scaled_twice_S, type="l", lty=1,
        col=cols[ok],
        xlab="Mutant trait (K1) value",
        ylab="Fitness derivative, scaled twice",
        ylim=c(-5,0))
abline(v=eq$x[1], lty=2)
abline(h=0, col="grey", lty=3)
lines(x_mutant, z_r2 / w_empty * eq$y, lty=2)
