## Fundamental result that we need: Species can coexist *at
## equilibrium* only if each species becomes limited by the resource
## for which it has, compared to its competitors, the highest
## requirement.

## Several species competing for one resource (p. 2684)
##
##   The species that has the lowest resource requirement for the
##   limiting resouce (i.e., the species with the lowest R^*_{1i}) will
##   displace all other species.
##
## So let's look at the behaviour of the system with a single resource
## and move towards understanding how the competition function might
## look like.
source("rstar_new.R")
source("rstar_support.R")

p <- rstar_parameters(rstar_mat_1, rstar_mat_1, S=1.0)
x <- matrix(0.5, nrow=2)
N0 <- 1

t <- seq(0, 30, length.out=201)

obj_tr <- rstar_run(t, x, N0, p)
obj_eq <- rstar_equilibrium(x, N0, p)

matplot(obj_tr$t, cbind(obj_tr$R, obj_tr$N), type="l", lty=1:2,
        xlab="Time", ylab="Abundance (red), Resource (black)")
abline(h=unlist(obj_eq[c("R", "N")]), col=1:2, lty=3)

## Test the resource equiliubrium:
dN <- obj_eq$N * 0.1

N1 <- obj_eq$N - dN
N2 <- obj_eq$N + dN

## I've made a mistake here somewhere: the equilibrium values should
## be 0.1861407 and 0.150721; surrounding the original value.

rstar_dRdt(x, obj_eq$N, obj_eq$R, p)
## This shows both increasing, which should not happen...
rstar_dRdt(x, N1, obj_eq$R, p) #  0.02083333
rstar_dRdt(x, N2, obj_eq$R, p) # -0.02083333

obj_tr_0 <- rstar_run_R(t, x, obj_eq$N, p)
obj_tr_1 <- rstar_run_R(t, x, N1, p)
obj_tr_2 <- rstar_run_R(t, x, N2, p)

obj_eq_1 <- rstar_equilibrium_R(x, N1, p)
obj_eq_2 <- rstar_equilibrium_R(x, N2, p)

matplot(t, cbind(obj_tr_0$R, obj_tr_1$R, obj_tr_2$R), type="l", col=1:3)
abline(h=c(obj_eq$R, obj_eq_1$R, obj_eq_2$R), col=1:3, lty=3)

## Next, we start working towards the instantaneous growth rate of a
## new type at this equilibrium

## Look at the fitness landscape: how does the instantaneous growth
## rate look with respect to K:
x_mutant <- seq(0, 1, length.out=101)
x_K <- rbind(x_mutant, obj_eq$x[2], deparse.level=0) # K varying
x_C <- rbind(obj_eq$x[2], x_mutant, deparse.level=0) # C varying

## We need R for the resident community: obj_eq$R
w_K <- rstar_fitness_given_R(x_K, obj_eq$R, p)
w_C <- rstar_fitness_given_R(x_C, obj_eq$R, p)

plot(x_mutant, w_K, type="l", xlab="Trait (K)", ylab="Fitness")
abline(h=0, col="grey", lty=3)
abline(v=obj_eq$x[1], lty=2)

## Fitness does not vary with respect to c when rare:
plot(x_mutant, w_C, type="l", xlab="Trait (C)", ylab="Fitness")
abline(h=0, col="grey", lty=3)

## Fitness in an empty environment:
plot(x_mutant, rstar_fitness_given_R(x_K, p$S, p), type="l",
     xlab="Trait (K)", ylab="Fitness in an empty environment")
lines(x_mutant, rstar_max_growth_rate(x_K, p), col="red")
abline(h=0, col="grey", lty=3)
abline(v=obj_eq$x[1], lty=2)
