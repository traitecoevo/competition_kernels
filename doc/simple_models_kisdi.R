source("common.R")

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
## species exerting and receiving competition, but here we do.  This
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
