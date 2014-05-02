source("common.R")

## # Rstar model, single resource

## This version of the model is based loosely off Huisman and Weissing
## 2001, which for the 1 resource case appears basically the same as
## the model on p. 46 (equation 5) of Tilman 1982.


## I'm going to consider the single resource case first because it's a
## bunch easier to think about than the two (or more!) resource case.
## It's also confusing enough to be worth thinking about.
##
## In this model, there can only be a single dominant species, so we
## can't study what happens during coexistence.  However, there's
## still plenty of competition going on.

## In the single resource case species are characterised by two
## traits: K (the constant of half-saturation of the Monod curve --
## higher K indicates higher requirements for the resource) and C (the
## rate of consumption of the resource).  All else being equal, the
## species with the lowest K will win because they will drive K below
## the minimum level required for positive net growth by the other
## species.

## +[MW: So ... are Huisman and Weissing modeling the resource here,
## rather than the population? My vague memory of Tilman 82 treatment
## was that the Monod curve described sensitivity of cell
## multiplication rate to resource concentration (50% at K), and there
## was a constant cell death rate that in effect defined R*. Plausible
## trade-off in this system was for strategies with low K to have high
## background death rate.]
##
## RF: Both the resource and the population are being modelled, and I
## think the treatment here is the same as your memory from Tilman
## (this is basically the same as the model on p. 46 (equation 5) of
## Tilman 1982).  The specific growth rate of each species depends on
## it's K, and then the consumption rate of the resource depends on
## how much is needed for growth, so K appears in both the equations
## for change in N and R.  See p. 2683 of the paper.  I do need to sit
## down and transcribe the equations into Revolve at some point
## though.
mat_r1 <- rstar_matrices(rstar_mat_1, rstar_mat_1)
m_r1 <- make_rstar(mat_r1, S=1)
sys_r1 <- list(x=matrix(0.5, nrow=2), y=1)

eq <- m_r1$single_equilibrium(sys_r1$x)

## Look at the fitness landscape: how does the instantaneous growth
## rate look with respect to mutant C and K.
x.mutant <- seq(0, 1, length=101)
x.K <- rbind(x.mutant, eq$x[2], deparse.level=0) # K varying
x.C <- rbind(eq$x[1], x.mutant, deparse.level=0) # C varying

## Here is the fitness landscape with respect to mutant K
##+ r1_fitness_landscape_K
plot(x.mutant, m_r1$fitness(x.K, eq$x, eq$y, eq$R), type="l",
     xlab="Trait (K)", ylab="Fitness")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

## Changing C doesn't change the fitness landscape because it only
## affects species at nontrivial densities.
##+ r1_fitness_landscape_C
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

## This is the fitness of a *single invader/mutant strategy*, as a
## function of the density of the resident.  A resident density
## implies a resource density (increasing residents decreases
## resources) and the fitness of the invading strategy decreases.

## Then compute the slope at the equilibrium density (indicated by the
## solid circle).
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
## is depressed by a unit increase in resident density.  The resident
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
w.empty <- m_r1$fitness(x.K, eq$x, 0)
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

## The other thing that changes the exact behaviour here is the
## resource supply rate, S, but that will only scale things up and
## down.

m_r1_Slo <- make_rstar(mat_r1, S=0.5)
m_r1_Shi <- make_rstar(mat_r1, S=2)

## The equilibrium density of the species varies, but equilibrium of
## the resource stays the same.
eq_Slo <- m_r1_Slo$single_equilibrium(sys_r1$x)
eq_Shi <- m_r1_Shi$single_equilibrium(sys_r1$x)

##+ r1_fitness_S
plot(x.mutant, m_r1$fitness(x.K, eq$x, eq$y, eq$R), type="l",
     xlab="Trait (K)", ylab="Fitness")
lines(x.mutant, m_r1_Slo$fitness(x.K, eq_Slo$x, eq_Slo$y, eq_Slo$R), col="red")
lines(x.mutant, m_r1_Shi$fitness(x.K, eq_Shi$x, eq_Shi$y, eq_Shi$R), col="blue")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

## The derivative with respect to density does change though, due to
## the nonlinear effect of density on fitness.
z_Slo <- jacobian(function(y) m_r1_Slo$fitness(x.K, eq$x, y), eq_Slo$y)
z_Shi <- jacobian(function(y) m_r1_Shi$fitness(x.K, eq$x, y), eq_Shi$y)

##+ r1_derivative_S
plot(x.mutant, z, type="l", ylim=range(z, z_Slo, z_Shi),
     xlab="Mutant trait (K) value",
     ylab="Fitness derivative")
lines(x.mutant, z_Slo, col="red")
lines(x.mutant, z_Shi, col="blue")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)
legend("topright", c("Low S", "Medium S", "High S"),
       col=c("red", "black", "blue"), lty=1, bty="n")

## As resources become less available, fitness becomes more absolutely
## sensitive to resident density.  This is consistent with the
## sublinear decrease of fitness with mutant density.

## Fitness in an empty environment also varies.
w.empty_Slo <- m_r1_Slo$fitness(x.K, eq$x, 0)
w.empty_Shi <- m_r1_Shi$fitness(x.K, eq$x, 0)

## As resources become more available, fitness increases.
##+ r1_empty_S
plot(x.mutant, w.empty, type="l", xlab="Mutant trait (K) value",
     ylab="Fitness in an empty environment")
lines(x.mutant, w.empty_Slo, col="red")
lines(x.mutant, w.empty_Shi, col="blue")
legend("topright", c("Low S", "Medium S", "High S"),
       col=c("red", "black", "blue"), lty=1, bty="n")

## Rescaling the derivative by fitness in an empty environment is a
## bit more surprising.  At high resource density (blue) this metric
## is also not monotonic (like the derivative calculation itself).  So
## perhaps this is not great after all (or perhaps that is actually
## OK).
##+ r1_derivative_scaled_S
plot(x.mutant, z / w.empty, type="l", ylim=range(z_Slo / w.empty_Slo),
     xlab="Mutant trait (K) value", ylab="Fitness derivative, scaled")
lines(x.mutant, z_Slo / w.empty_Slo, col="red")
lines(x.mutant, z_Shi / w.empty_Shi, col="blue")
abline(v=eq$x[1], lty=2)
abline(h=0, col="grey", lty=3)
legend("bottomleft", c("Low S", "Medium S", "High S"),
       col=c("red", "black", "blue"), lty=1, bty="n")

## ## Varying resident density

## Another assumption made by most explicit competition models is that
## "competitive effect" scales linearly with density of a species.
## From the nonlinear fitness / resident density plot above, this is
## clearly not the case.  We'd still expect that fitness should at
## least be monotonic decreasing with density, though.
mat_r1 <- rstar_matrices(rstar_mat_1, rstar_mat_1)
m_r1 <- make_rstar(mat_r1, S=1)
sys_r1 <- list(x=matrix(0.5, nrow=2), y=1)
eq <- m_r1$single_equilibrium(sys_r1$x)

## Build a vector of resident densities as a logarithmic series around
## the equilibrium density.
y.resident <- 2^seq(-6, 2) * eq$y

## Then we want to know what the fitness is at these levels.  I'm
## using `x.K` and `x.C` from above.

## Fitness in an empty environment, and at equilibrium resident
## density:
w.empty    <- m_r1$fitness(x.K, eq$x, 0)
w.resident <- m_r1$fitness(x.K, eq$x, eq$y)

## Fitness
w.varying <- sapply(y.resident, function(y)
                    m_r1$fitness(x.K, eq$x, y))

## A vector of colours for plotting
cols <- RColorBrewer::brewer.pal(length(y.resident), "Blues")

## Here, as the density of residents increases (darker blues),
## invasion fitness decreases, and the decline in fitness becomes
## increasingly concave up with respect to the mutant K.  The upper
## dashed line is fitness in an empty environment, and the lower
## blue/black dashed line is at equilibrium resident density.
##+ r1_fitness_varying
plot(w.empty ~ x.mutant, type="l", lty=2, ylim=range(w.varying),
     xlab="Mutant K", ylab="Fitness")
abline(h=0, lty=3, col="grey")
abline(v=eq$x[[1]], lty=2)
matlines(x.mutant, w.varying, lty=1, col=cols)
lines(x.mutant, w.resident, lty=2)

## Scaling fitness against that in an empty environment (so this is
## the proportional decrease in fitness).  As resident density
## increases, the biggest fitness decreases are felt by mutants with
## the biggest K values.
##+ r1_fitness_varying_scaled
w.varying.scaled <- w.varying / w.empty
plot(w.empty/w.empty ~ x.mutant, type="l", lty=2,
     ylim=range(w.varying.scaled), xlab="Mutant K", ylab="Fitness")
abline(h=0, lty=3, col="grey")
abline(v=eq$x[[1]], lty=2)
matlines(x.mutant, w.varying.scaled, lty=1, col=cols)
lines(x.mutant, w.resident / w.empty, lty=2)

## Next, compute the competition estimate (derivatrive) at each of
## these resident denities:
z.varying <- sapply(y.resident, function(y.res)
                    jacobian(function(y) m_r1$fitness(x.K, eq$x, y), y.res))
## And at equilibrium
z.resident <- jacobian(function(y) m_r1$fitness(x.K, eq$x, y), eq$y)
## (The gradient in the empty environment is harder to compute because
## we have to move up from zero only.  We should get close enough with
## the hack below)
eps <- sqrt(.Machine$double.eps)
z.empty <- jacobian(function(y) m_r1$fitness(x.K, eq$x, y), eps,
                    method="simple", method.args=list(eps=eps))

## The dashed line is the derivative in an empty environment, while
## the dashed black/blue line is derivative at resident equilibrium.
##
## As the resident density increases, the mutant type that is most
## strongly affected by the resident (most negative derivative) has a
## lower and lower K value.
##+ r1_derivative_varying
plot(z.empty ~ x.mutant, type="l", lty=2, ylim=range(z.varying),
     xlab="Mutant K", ylab="Fitness derivative")
abline(h=0, lty=3, col="grey")
abline(v=eq$x[[1]], lty=2)
matlines(x.mutant, z.varying, lty=1, col=cols)
lines(x.mutant, z.resident, lty=2)

## Earlier, normalising the derivative by fitness in an empty
## environment seemed to make sense:
z.varying.scaled <- z.varying / w.empty

##+ r1_derivative_varying_scaled
plot(z.empty / w.empty ~ x.mutant, type="l", lty=2,
     ylim=range(z.varying.scaled),
     xlab="Mutant K", ylab="Scaled fitness derivative")
abline(h=0, lty=3, col="grey")
abline(v=eq$x[[1]], lty=2)
matlines(x.mutant, z.varying.scaled, lty=1, col=cols)
lines(x.mutant, z.resident/w.empty, lty=2)

## We can rescale these values by the total resident density, which
## seems to put everything onto a similar scale.
### TODO: Go through and work out units here so that I can work out
### wnat this actually means.
##+ r1_derivative_varying_scaled_resident
z.varying.scaled.resident <- t(t(z.varying) * y.resident)
plot(z.empty * 0 ~ x.mutant, type="l", lty=2,
     ylim=range(z.varying.scaled.resident),
     xlab="Mutant K",
     ylab="Fitness derivative, scaled against resident density")
abline(h=0, lty=3, col="grey")
abline(v=eq$x[[1]], lty=2)
matlines(x.mutant, z.varying.scaled.resident, lty=1, col=cols)
lines(x.mutant, z.resident * eq$y, lty=2)

## Note that this scaling only matters in the context of comparing
## different resident densities, because we've not really worked out
## what the y axis means.

## And we can scale both ways:
z.varying.scaled.twice <- t(t(z.varying.scaled) * y.resident)

##+ r1_derivative_varying_scaled_twice
plot(z.empty / w.empty * 0 ~ x.mutant, type="l", lty=2,
     ylim=range(z.varying.scaled.twice),
     xlab="Mutant K",
     ylab="Fitness derivative, scaled twice")
abline(h=0, lty=3, col="grey")
abline(v=eq$x[[1]], lty=2)
matlines(x.mutant, z.varying.scaled.twice, lty=1, col=cols)
lines(x.mutant, z.resident / w.empty * eq$y, lty=2)

## If these plots capture something reasonable about competition, then
## it means that competition in this model is not monotonic in
## resident density, let alone linear.

## Let's take this a little further.
y.resident.hires <- eq$y * 2^seq(-6, 2, length=51)
z.varying.hires <- sapply(y.resident.hires, function(y)
                          model_jacobian_density(x.K, sys(eq$x, y), m_r1))
z.varying.scaled.hires <- z.varying.hires / w.empty
z.varying.scaled.twice.hires <- t(t(z.varying.scaled.hires) *
                                  y.resident.hires)

## Heatmap of the derivative (dark red is most negative, light yellow
## is zero).
##+ r1_derivative_image
image(x.mutant, y.resident.hires, z.varying.hires, log="y",
      xlab="Mutant trait (K)", ylab="Resident density")
abline(v=eq$x[1], h=eq$y, lty=2)

## Heatmap of the derivative, scaled against fitness in the empty
## landscape.
##+ r1_derivative_scaled_image
image(x.mutant, y.resident.hires, z.varying.scaled.hires, log="y",
      xlab="Mutant trait (K)", ylab="Resident density")
abline(v=eq$x[1], h=eq$y, lty=2)

## Heatmap of the derivative, scaled against fitness in the empty
## landscape, and against resident density.
##+ r1_derivative_scaled_twice_image
image(x.mutant, y.resident.hires, z.varying.scaled.twice.hires, log="y",
      xlab="Mutant trait (K)", ylab="Resident density")
abline(v=eq$x[1], h=eq$y, lty=2)

##+ r1_derivative_3d, fig.height=7
persp(x.mutant, log(y.resident.hires), z.varying.hires,
      theta=30, phi=30, shade=0.5, col="green3", border=NA,
      xlab="Mutant trait (K)", ylab="Resident density (log scale)",
      zlab="Fitness derivative")

##+ r1_derivative_scaled_3d, fig.height=7
persp(x.mutant, log(y.resident.hires), z.varying.scaled.hires,
      theta=30, phi=30, shade=0.5, col="green3", border=NA,
      xlab="Mutant trait (K)", ylab="Resident density (log scale)",
      zlab="Fitness derivative, scaled")

## (note the change in perspective here)
##+ r1_derivative_scaled_twice_3d, fig.height=7
persp(x.mutant, log(y.resident.hires), z.varying.scaled.twice.hires,
      theta=55, phi=40, shade=0.5, col="green3", border=NA,
      xlab="Mutant trait (K)", ylab="Resident density (log scale)",
      zlab="Fitness derivative, scaled twice")

knitr::knit_hooks$set(document=function(x) browser())
knitr::knit("test.Rmd")
writeLines(readLines("test.md"))
