library(Revolve)
library(plyr)
library(numDeriv)
library(RColorBrewer)

mat <- huisman_matrices(huisman_mat_1, huisman_mat_1)

m <- make_huisman_2001(mat, S=1)
x <- matrix(0.5, nrow=2)

eq <- m$single_equilibrium(x)

# First, change in the density of the *resident* species.  This one is
# comparatively easy; everthing to do this already sorted out.

# Here is a vector of resident densities; the third to last of these
# is the actual equilibrium density.
yy <- 2^seq(-6, 2) * eq$y
cols <- RColorBrewer::brewer.pal(length(yy), "Blues")

# Then we want to know what the fitness is at these levels.  That's
# actually quite easy too.  As with derivatives-huisman-1.R, look at
# varying the K parameter and hold the C parameter.
xx <- seq(0, 1, length=201)
x.K <- rbind(xx, x[2], deparse.level=0) # K varying

w.empty <- m$fitness(x.K, x, 0)
w.resident <- m$fitness(x.K, x, eq$y)

ans.w <- sapply(yy, function(y) m$fitness(x.K, x, y))

plot(w.empty ~ xx, type="l", lty=2, ylim=range(ans.w),
     xlab="Mutant K", ylab="Fitness")
abline(h=0, lty=3, col="grey")
abline(v=x[[1]], lty=2)
matlines(xx, ans.w, lty=1, col=cols)
lines(xx, w.resident, lty=2)

# Scaled against the fitness in an empty environment
plot(w.empty/w.empty ~ xx, type="l", lty=2, ylim=range(ans.w/w.empty),
     xlab="Mutant K", ylab="Fitness")
abline(h=0, lty=3, col="grey")
abline(v=x[[1]], lty=2)
matlines(xx, ans.w/w.empty, lty=1, col=cols)
lines(xx, w.resident / w.empty, lty=2)

# Now, compute the competition estimate (derivative) at each of these
# densities:
ans.z <- sapply(yy, function(y.res)
                jacobian(function(y) m$fitness(x.K, x, y), y.res))
z.resident <- jacobian(function(y) m$fitness(x.K, x, y), eq$y)

plot(rep(0, length(xx)) ~ xx, type="l", lty=2, ylim=range(ans.z),
     xlab="Mutant K", ylab="Fitness derivative")
abline(h=0, lty=3, col="grey")
abline(v=x[[1]], lty=2)
matlines(xx, ans.z, lty=1, col=cols)
lines(xx, z.resident, lty=2)

# Rescale against resident density
tmp <- t(t(ans.z) * yy)

plot(rep(0, length(xx)) ~ xx, type="l", lty=2, ylim=range(tmp),
     xlab="Mutant K", ylab="Fitness derivative")
abline(h=0, lty=3, col="grey")
abline(v=x[[1]], lty=2)
matlines(xx, tmp, lty=1, col=cols)
lines(xx, z.resident * eq$y, lty=2)

# Scaled by the growth in an empty environment
plot(rep(0, length(xx)) ~ xx, type="l", lty=2, ylim=range(ans.z/w.empty),
     xlab="Mutant K", ylab="Scaled fitness derivative")
abline(h=0, lty=3, col="grey")
abline(v=x[[1]], lty=2)
matlines(xx, ans.z/w.empty, lty=1, col=cols)
lines(xx, z.resident/w.empty, lty=2)

tmp2 <- t(t(ans.z / w.empty) * yy)
plot(rep(0, length(xx)) ~ xx, type="l", lty=2, ylim=range(tmp2),
     xlab="Mutant K", ylab="Fitness derivative")
abline(h=0, lty=3, col="grey")
abline(v=x[[1]], lty=2)
matlines(xx, tmp2, lty=1, col=cols)
lines(xx, z.resident / w.empty * eq$y, lty=2)

# So if this is true, then competition is not actually monotonic in
# resident density, let alone linear.

# Let's take this a little further.
yyy <- tree::seq_log(eq$y * 2^(-6), eq$y * 2^2, length=51)

res.z <- sapply(yyy, function(y.res)
                jacobian(function(y) m$fitness(x.K, x, y), y.res))
image(xx, yyy, res.z, log="y")
image(xx, yyy, res.z / w.empty, log="y")

persp(xx, yyy, res.z, theta=30, phi=30, shade=TRUE,
      col="steelblue", border=NA)

# Then, see what happens as we change the density of the mutant type,
# holding the resident at equilibrium.  This one is a bit more messed
# up.


