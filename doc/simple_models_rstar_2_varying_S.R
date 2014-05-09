source("common.R")

## # Rstar model, two resource, with varying supply rates

## Once I get a simple story here, I'll move the relevant details into
## [the main file](simple_models_rstar_2).  But this is sufficiently
## complicated to warrant it's own file.  Similar explorations are
## going to be needed for the resident consumption vector
## exploration, which will be in
## [this file](simple_models_rstar_2_varying_S).

## This will be easier to deal with if we consider a couple of
## different vectors first.

## Picking up from the end of [the main file](simple_models_rstar_2):
get_S <- function(m) {
  m$parameters$get()[["S"]]
}
mat_r2 <- rstar_matrices(rstar_mat_2_tradeoff,
                         make_rstar_mat_constant(matrix(0.5, 2)))
m_r2 <- make_rstar(mat_r2, S=c(0.5, 0.5))
sys_r2 <- sys(matrix(0.3), 1)
eq <- m_r2$single_equilibrium(sys_r2$x)
x_mutant <- matrix(seq(0, 1, length=301), 1)
w_mutant <- m_r2$fitness(x_mutant, eq$x, eq$y)
w_empty <- m_r2$fitness(x_mutant, eq$x, 0, get_S(m_r2))
z_r2 <- model_jacobian_density(x_mutant, eq, m_r2)

## ## Varying supply rates

## Things start getting a bit more complicated once start allowing the
## resource supply rates, S, to change.

## Following Tilman, I'll move the resources along a gradient holding
## S1 + S2 = 1.  What matters for invasion is which side of the C
## vector the starting point lands on, so this should be sufficient.

## The resident species will survive wherever the S line falls above
## the ZNGIs.
Rstar <- m_r2$Rstar(sys_r2$x)

c.slope <- 1 # exp(drop(diff(log(mat_r2$C(eq$x))))), I think.
c.intercept <- Rstar[2] - c.slope * Rstar[1]
S.max <- (1 - c.intercept) / (1 + c.slope)
S.crit <- c(Rstar[1], 1 - Rstar[2])

S1 <- sort(c(seq(0, 1, length=51), S.crit, S.max))
S <- cbind(S1=S1, S2=1 - S1)

## I want a vector of colours along this supply vector, with red being
## very low supply of R1 and blue being very high supply of R1.  
red.blue <- colorRampPalette(RColorBrewer::brewer.pal(4, "RdBu"))
cols <- red.blue(length(S1))

## It's going to be useful to know if the resident can survive, so
## code up line types, colours and point styles for this:
ok <- apply(t(S) > drop(Rstar), 2, all)
cols.crit <- ifelse(ok, cols, make_transparent(cols, .25))
lty.crit  <- ifelse(ok, 1, 2)
pch.crit  <- ifelse(ok, 19, 1)

##+ r2_rstar_S
rstar_plot(m_r2, sys_r2)
abline(1, -1, lty=3)
points(S1, 1 - S1, col=cols, pch=pch.crit, lwd=1.5)

## Next, consider the equilibrium density of residents as the resource
## supply vector changes:
m_r2_S <- apply(S, 1, function(S) make_rstar(mat_r2, S=S))
eq_r2_S <- lapply(m_r2_S, function(m) m$single_equilibrium(eq$x))
y_r2_S <- sapply(eq_r2_S, "[[", "y")

## Resident density as a function of the supply rate of the first
## resource (when it directly trades off with the second).  The
## optimum resident density occurs when the supply point is at S.max
##+ r2_fitness_resident_S
plot(S1, y_r2_S, type="l",
     xlab="Supply rate of resource 1", ylab="Resident density")
points(S1, y_r2_S, pch=pch.crit, col=cols, lwd=1.5)
abline(v=S.crit, lty=3)
abline(v=0.5, lty=2)
abline(v=S.max, lty=2, col="blue")

## There are a couple of cases worth considering:
##   * a: S1 below the the smallest critical value
##   * b: ...above the critical value but above the optimum
##   * c: ...at optimum
##   * d: ...at S1 = S2 = 0.5, the same as the previous case
##   * e: ...above 0.5, but below the larger critical value
##   * f: ...above the critical value

S1 <- c(a=mean(c(S.crit[1], 0)),
        b=mean(c(S.crit[1], min(S.max, 0.5))),
        c=min(S.max, 0.5),
        d=max(S.max, 0.5),
        e=mean(c(S.crit[2], max(S.max, 0.5))),
        f=mean(c(S.crit[2], 1)))
S <- cbind(D1=S1, S2 = 1 - S1)

cols <- red.blue(length(S1))
cols[S1 == 0.5] <- 1

rstar_plot(m_r2, sys_r2)
abline(1, -1, lty=3)
abline(0, 1, lty=3, col="grey")
points(S1, 1 - S1, col=cols, pch=19, lwd=1.5)
text(S1, 1 - S1, names(S1), adj=c(.5, 2))

## Lots of colour tweaking here...
ok <- apply(t(S) > drop(Rstar), 2, all)
cols.crit <- ifelse(ok, cols, make_transparent(cols, .25))
lty.crit  <- ifelse(ok, 1, 2)
pch.crit  <- ifelse(ok, 19, 1)
i_max <- S1 == S.max
lwd.crit <- ifelse(i_max, 10, 1)
cols.crit[i_max] <- make_transparent(cols[i_max], 0.3)

## Rebuild the models, equilibriums and get the resident densities:
m_r2_S <- apply(S, 1, function(S) make_rstar(mat_r2, S=S))
eq_r2_S <- lapply(m_r2_S, function(m) m$single_equilibrium(eq$x))
y_r2_S <- sapply(eq_r2_S, "[[", "y")

## So, outside of the critical S values determined by the resident
## Rstar, measures of competition don't make any sense because there
## is no resident species to have a competitive effect on invaders.
## This might cause issues later...
w_r2_S <- mapply(function(m, eq)
                 m$fitness(x_mutant, eq$x, eq$y, eq$R),
                 m_r2_S, eq_r2_S)

## Here is a plot of mutant fitness as a function of K1 value for a
## range of different resource supply ratios.  Red lines have
## relatively high rates of S1, enabling the carrying capacity of the
## resident to be higher.  Blue lines have S2 > S1, and relatively
## lower densities of the resident.  Basically Where the resident can
## exist it *will* draw the resource down to Rstar, and the fitness
## function will then pass through that point.
##+ r2_fitness_S
matplot(drop(x_mutant), w_r2_S, type="l", col=cols.crit, lty=lty.crit,
        lwd=lwd.crit,
        xlab="Mutant trait (K1) value",
        ylab="Invader fitness")
abline(h=0, lty=3)
abline(v=eq$x[1], lty=2)
abline(v=0.5, lty=2, col="grey")
lines(x_mutant, w_mutant, lty=2)

## Fitness in an empty landscape
w_r2_empty_S <- mapply(function(m, eq)
                       m$fitness(x_mutant, eq$x, 0, get_S(m)),
                       m_r2_S, eq_r2_S)

## As the S1 rate decreases (redder lines), the optimal K1 value
## decreases.  The dotted lines do are negative at the point that our
## resident is at -- the resident cannot survive here.
##+ r2_fitness_empty_S
matplot(drop(x_mutant), w_r2_empty_S, type="l",
        col=cols.crit, lty=lty.crit, lwd=lwd.crit,
        xlab="Mutant trait (K1) value",
        ylab="Fitness in empty environment")
abline(h=0, lty=3)
abline(v=eq$x[1], lty=2)

## Here's the *depression* in fitness -- the total decrease in fitness
## that is felt in the presence of the resident.  There are only four
## lines here because the two lines that correspond to the resident
## being absent don't belong.
##+ fitness_depression
matplot(drop(x_mutant), w_r2_empty_S - w_r2_S, type="l",
        col=cols.crit, lty=lty.crit, lwd=lwd.crit,
        xlab="Mutant trait (K1) value",
        ylab="Depression in fitness")
abline(h=0, lty=3)
abline(v=eq$x[1], lty=2)
abline(v=S.crit, lty=3)

## This is the *relative depression* in fitness -- the fractional
## decrease in fitness due to the presence of the resident.  It takes
## into account the fact that relatively fit species have further to
## fall.
matplot(drop(x_mutant), w_r2_S / w_r2_empty_S , type="l",
        col=cols.crit, lty=lty.crit, lwd=lwd.crit,
        xlab="Mutant trait (K1) value",
        ylab="Relative depression in fitness", ylim=c(-4, 1))
abline(h=0, lty=3)
abline(v=eq$x[1], lty=2)
abline(v=S.crit, lty=3)

## Where a value is above zero, the mutant can invade in the presence
## of the resident - that's why all lines have a y value of zero at
## the resident K1 value.  If a line was above 1, that would represent
## facilitation -- the mutant does better in the presence of the
## resident.  That is not really possible in this model though, but
## there are some pathalogical cases here where negative
## empty-environment fitness can look like very strong facilitation.

## At an even supply rate (S1 = S2), the fitness depression is zero
## around K1 = 0.5.  As the supply rate moves towards the optimum for
## the resident fitness depression becomes (roughly) symetric around
## the resident trait.

## At S1 supply rates less than the resident optimum, or greater than
## 0.5, fitness depression becomes increasingly asymmetric until the
## resident will outcompete all species with a K1 bigger than it (for
## very low S1) or outcompete all species with a K1 smaller than it
## (for very high S1).

## Compute the fitness derivatives with respect to resident density.
## This can only be done easily where the equilibrium resident density
## is at least zero.  We can use the same hack as above for zero
## equilibrium though.
z_r2_S <- mapply(function(m, eq)
                 model_jacobian_density(x_mutant, eq, m),
                 m_r2_S[ok], eq_r2_S[ok])

##+ r2_derivative_S1
matplot(drop(x_mutant), z_r2_S, col=cols.crit[ok], type="l", lty=1,
        lwd=lwd.crit[ok],
        xlab="Mutant trait (K1) value", ylab="Fitness derivative")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

## Below here is unedited from before, but needs some work...


## Before considering all the points at once, let's just consider a
## single point.
S_cmp <- c(0.6, 0.4)
m_cmp <- make_rstar(mat_r2, S=S_cmp)
eq_cmp <- m_cmp$single_equilibrium(eq$x)
z_cmp <- model_jacobian_density(x_mutant, eq_cmp, m_cmp)

##+ fitness_derivative_S_1
matplot(drop(x_mutant), cbind(z_r2, z_cmp), type="l",
        xlab="Mutant trait (K1) value",
        ylab="Fitness derivative")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

## This is really weird - what's the discontinuity in the
## uneven-supply case coming from?  The derivative is always negative
## (so the resident is always decreasing growth) but it does so in a
## way that species with a K1 value that is smaller than some value
## are *much* more affected by the resident than species with a K1
## value just below the threshhold.

## This discontinuity happens at small deviations from the even supply
## vector (here a deviation of 0.01, rather than the 0.1 used above).
##+ fitness_derivative_S_2
S_cmp2 <- c(0.5, 0.5) - c(-1, 1) * 0.01
m_cmp2 <- make_rstar(mat_r2, S=S_cmp2)
eq_cmp2 <- m_cmp2$single_equilibrium(eq$x)
z_cmp2 <- model_jacobian_density(x_mutant, eq_cmp2, m_cmp2)
matplot(drop(x_mutant), cbind(z_r2, z_cmp, z_cmp2), type="l",
        xlab="Mutant trait (K1) value",
        ylab="Fitness derivative")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

## So, what is going on there?
##+ cache=TRUE
yy_resident <- eq_cmp2$y + seq(-0.1, 0.1, length=201)
xx_mutant <- matrix(seq(0.46, 0.56, length=221), 1)
ww_mutant <- sapply(yy_resident, function(y)
                    m_cmp2$fitness(xx_mutant, eq_cmp2$x, y))
zz_mutant <- model_jacobian_density(xx_mutant, eq_cmp2, m_cmp2)
image(drop(xx_mutant), yy_resident, ww_mutant, xlab="Mutant trait",
      ylab="Resident density")
abline(h=eq_cmp2$y)

## Hmm, I need to think more about that.  But it looks like there's a
## step change from one regime to another.

## I think what is going on is based on this view of the same data:
matplot(drop(xx_mutant), ww_mutant, type="l",
        col=red.blue(length(yy_resident)), lty=1,
        xlab="Mutant trait (K1) value", ylab="Mutant fitness")

## Proceeding across the entire range of supply values anyway:

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

## I think I need to think about this more graphically -- who are the
## species that are getting most affected, and who are being less
## affected.

## Then we can rescale this by the fitness in an empty environment.
## We might have to exclude these trait values (which vary as a
## function of S).
z_r2_scaled_S <- z_r2_S / w_r2_empty_S[,ok]
z_r2_scaled_twice_S <- t(t(z_r2_scaled_S) * y_r2_S[ok])

## There are some really weird discontinuities here that I need to
## chase up.  I think that they're caused by calculations where
## fitness switches sign.
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

mask <- function(x, col="#ffffffcc", lty=3) {
  usr <- par("usr")
  rect(x, usr[3], usr[1:2], usr[4], col=col, border=NA)
  abline(v=x, lty=lty)
}
