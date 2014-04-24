library(Revolve)
library(numDeriv)
library(RColorBrewer)

mat <- rstar_matrices(rstar_mat_1, rstar_mat_1)

m <- make_rstar(mat, S=1)
sys0 <- sys(matrix(0.5, nrow=2), 1)

eq <- m$single_equilibrium(sys0$x)

# First, change in the density of the *resident* species.  This one is
# comparatively easy; everthing to do this already sorted out.

# Here is a vector of resident densities; the third to last of these
# is the actual equilibrium density.
yy <- 2^seq(-6, 2) * eq$y
cols <- RColorBrewer::brewer.pal(length(yy), "Blues")

# Then we want to know what the fitness is at these levels.  That's
# actually quite easy too.  As with derivatives-rstar-1.R, look at
# varying the K parameter and hold the C parameter.
x.mutant <- seq(0, 1, length=201)
x.K <- rbind(x.mutant, eq$x[2], deparse.level=0) # K varying

w.empty <- m$fitness(x.K, eq$x, 0)
w.resident <- m$fitness(x.K, eq$x, eq$y)

ans.w <- sapply(yy, function(y) m$fitness(x.K, eq$x, y))

plot(w.empty ~ x.mutant, type="l", lty=2, ylim=range(ans.w),
     xlab="Mutant K", ylab="Fitness")
abline(h=0, lty=3, col="grey")
abline(v=eq$x[[1]], lty=2)
matlines(x.mutant, ans.w, lty=1, col=cols)
lines(x.mutant, w.resident, lty=2)

# Scaled against the fitness in an empty environment
plot(w.empty/w.empty ~ x.mutant, type="l", lty=2, ylim=range(ans.w/w.empty),
     xlab="Mutant K", ylab="Fitness")
abline(h=0, lty=3, col="grey")
abline(v=eq$x[[1]], lty=2)
matlines(x.mutant, ans.w/w.empty, lty=1, col=cols)
lines(x.mutant, w.resident / w.empty, lty=2)

# Now, compute the competition estimate (derivative) at each of these
# densities:
ans.z <- sapply(yy, function(y.res)
                jacobian(function(y) m$fitness(x.K, eq$x, y), y.res))
z.resident <- jacobian(function(y) m$fitness(x.K, eq$x, y), eq$y)

plot(rep(0, length(x.mutant)) ~ x.mutant, type="l", lty=2, ylim=range(ans.z),
     xlab="Mutant K", ylab="Fitness derivative")
abline(h=0, lty=3, col="grey")
abline(v=eq$x[[1]], lty=2)
matlines(x.mutant, ans.z, lty=1, col=cols)
lines(x.mutant, z.resident, lty=2)

# Rescale against resident density
tmp <- t(t(ans.z) * yy)

plot(rep(0, length(x.mutant)) ~ x.mutant, type="l", lty=2, ylim=range(tmp),
     xlab="Mutant K", ylab="Fitness derivative")
abline(h=0, lty=3, col="grey")
abline(v=eq$x[[1]], lty=2)
matlines(x.mutant, tmp, lty=1, col=cols)
lines(x.mutant, z.resident * eq$y, lty=2)

# Scaled by the growth in an empty environment
plot(rep(0, length(x.mutant)) ~ x.mutant, type="l", lty=2, ylim=range(ans.z/w.empty),
     xlab="Mutant K", ylab="Scaled fitness derivative")
abline(h=0, lty=3, col="grey")
abline(v=eq$x[[1]], lty=2)
matlines(x.mutant, ans.z/w.empty, lty=1, col=cols)
lines(x.mutant, z.resident/w.empty, lty=2)

tmp2 <- t(t(ans.z / w.empty) * yy)
plot(rep(0, length(x.mutant)) ~ x.mutant, type="l", lty=2, ylim=range(tmp2),
     xlab="Mutant K", ylab="Fitness derivative")
abline(h=0, lty=3, col="grey")
abline(v=eq$x[[1]], lty=2)
matlines(x.mutant, tmp2, lty=1, col=cols)
lines(x.mutant, z.resident / w.empty * eq$y, lty=2)

# So if this is true, then competition is not actually monotonic in
# resident density, let alone linear.

# Let's take this a little further.
seq_log <- function(from, to, length) {
  exp(seq(log(from), log(to), length.out=length))
}
yyy <- seq_log(eq$y * 2^(-6), eq$y * 2^2, length=51)

res.z <- sapply(yyy, function(y.res)
                jacobian(function(y) m$fitness(x.K, eq$x, y), y.res))
image(x.mutant, yyy, res.z, log="y")
image(x.mutant, yyy, res.z / w.empty, log="y")

persp(x.mutant, yyy, res.z, theta=30, phi=30, shade=TRUE,
      col="steelblue", border=NA)
persp(x.mutant, yyy, res.z / w.empty, theta=30, phi=30, shade=TRUE,
      col="steelblue", border=NA)

# Then, see what happens as we change the density of the mutant type,
# holding the resident at equilibrium.  This one is a bit more messed
# up.
