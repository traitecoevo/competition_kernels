## Issues and improvements:
##   1: Why is the optimal strategy all the way over at 0.9
##   2: Fix the break point in simulations
##   3: Different measures in main paper or on SM?
##   4: Wire up C as well as K so that coexistence is possible, so
##   that figures are more directly comparable with tree.

## 6 panels:
## column 1: R* zngi and trajectories with residents
## column 2: CC for resident 1
## column 3: CC for resident 2
## row1: symmetric C
## row2: asymmetric C

source("rstar_new.R")
source("rstar_support.R")
source("rstar_plot.R")
source("support.R")

rstar_competition <- function(x_resident, p, n=301L) {
  x_invade <- rbind(seq(0, 1, length.out=n))

  eq <- rstar_single_equilibrium(x_resident, p)
  K_invade <- rstar_carrying_capacity(x_invade, p)
  r_invade <- rstar_max_growth_rate(x_invade, p)
  w_invade <- rstar_fitness_given_R(x_invade, eq$R, p)
  
  list(p=p,
       x_resident=x_resident,
       x_invade=x_invade,
       N_resident=eq$N,
       R_eq=eq$R,
       K_invade=K_invade,
       w_invade=w_invade,
       alpha=competition_LV(w_invade, r_invade, K_invade, eq$N))
}

C1 <- c(.2, .2) * 10
C2 <- c(.3, .7)
x_resident1_1 <- matrix(0.6, nrow=1)
x_resident1_2 <- matrix(0.5, nrow=1)
x_resident2_1 <- matrix(0.7, nrow=1)
x_resident2_2 <- matrix(0.9, nrow=1)

x_invade <- rbind(seq(0, 1, length.out=301))

p1 <- rstar_parameters(rstar_mat_2_tradeoff, matrix(C1, nrow=2), S=0.25)
p2 <- rstar_parameters(rstar_mat_2_tradeoff, matrix(C2, nrow=2), S=0.5)

## obji_j is model i, resident j
obj1_1 <- rstar_competition(x_resident1_1, p1)
obj1_2 <- rstar_competition(x_resident1_2, p1)
obj2_1 <- rstar_competition(x_resident2_1, p2)
obj2_2 <- rstar_competition(x_resident2_2, p2)

xx <- seq(0, 1, length.out=6)

## TODO: probably better in absolute size?
black_bar <- function(x, y, col="black") {
  usr <- par("usr")
  rect(x[1], usr[3], x[2], usr[3] + y, col=col, border=col)
}

ylim_alpha <- c(.45, 1.75)
ybar <- 0.05

par(mfrow=c(2, 3), mar=rep(.5, 4), oma=c(2, 2, 0, 0))

plot(NA, xlim=0:1, ylim=0:1, xlab="", ylab="", xaxt="n", yaxt="n")
for (i in xx) {
  rstar_lines(rbind(i), p1, col="grey")
}
rstar_lines(obj1_1$x_resident, obj1_1$p, col="red")
rstar_lines(obj1_2$x_resident, obj1_2$p, col="blue")
abline(0, 1, lty=2)

plot(obj1_1$x_invade, obj1_1$alpha, type="l", ylim=ylim_alpha)
abline(h=1.0, v=obj1_1$x_resident, lty=2, col="darkgrey")
points(obj1_1$x_resident, 1.0, pch=19)
black_bar(range(obj1_1$x_invade[obj1_1$w_invade > 0]), ybar)

plot(obj1_2$x_invade, obj1_2$alpha, type="l", ylim=ylim_alpha)
abline(h=1.0, v=obj1_2$x_resident, lty=2, col="darkgrey")
points(obj1_2$x_resident, 1.0, pch=19)

plot(NA, xlim=0:1, ylim=0:1, xlab="", ylab="", xaxt="n", yaxt="n")
for (i in xx) {
  rstar_lines(rbind(i), p2, col="grey")
}
rstar_lines(obj2_1$x_resident, obj2_1$p, col="red")
rstar_lines(obj2_2$x_resident, obj2_2$p, col="blue")
R <- drop(rstar_Rstar(x_resident2_2, p2))
C <- p2$C(x_resident2_2)
abcline(R[1], R[2], C[2] / C[1], lty=2)

plot(obj2_1$x_invade, obj2_1$alpha, type="l", ylim=ylim_alpha)
abline(h=1.0, v=obj2_1$x_resident, lty=2, col="darkgrey")
points(obj2_1$x_resident, 1.0, pch=19)
black_bar(range(obj2_1$x_invade[obj2_1$w_invade > 0]), ybar)

plot(obj2_2$x_invade, obj2_2$alpha, type="l", ylim=ylim_alpha)
abline(h=1.0, v=obj2_2$x_resident, lty=2, col="darkgrey")
points(obj2_2$x_resident, 1.0, pch=19)
