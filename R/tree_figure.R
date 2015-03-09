fig_tree_lma <- function(d1, d2) {
  ybar <- 0.05
  par(mfrow=c(1, 2), mar=c(4.1, 1, .5, .5), oma=c(0, 2, 0, 0))
  plot(d1$x_mutant, d1$alpha, log="x", type="l", ylim=c(0, 1.5), las=1,
       xlab="Leaf mass per area (lma; kg / m^2)",
       ylab="Competition")
  points(d1$x_resident, 1.0, pch=19)
  abline(v=d1$x_resident, h=1.0, lty=2)
  black_bar(range(d1$x_resident, d1$x_mutant[d1$w > 0]), ybar)

  ## TODO: The bifurcation will get missed on this plot because of
  ## using range()
  plot(d2$x_mutant, d2$alpha, log="x", type="l", ylim=c(0, 1.5), las=1,
       xlab="Leaf mass per area (lma; kg / m^2)",
       ylab="", yaxt="n")
  axis(2, labels=FALSE)
  points(d2$x_resident, 1.0, pch=19)
  abline(v=d2$x_resident, h=1.0, lty=2)
  black_bar(range(d2$x_resident, d2$x_mutant[d2$w > 0]), ybar)
}

fig_tree_hmat <- function(d1, d2) { 
  ybar <- 0.05 
  ylim <- range(d1$alpha, d2$alpha)
  par(mfrow=c(1, 2), mar=c(4.1, 1, .5, .5), oma=c(0, 2, 0, 0))
  plot(d1$x_mutant, d1$alpha, log="x", type="l", ylim=ylim, las=1,
       xlab="Height at maturation (m)",
       ylab="Competition")
  points(d1$x_resident, 1.0, pch=19)
  abline(v=d1$x_resident, h=1.0, lty=2)
  black_bar(range(d1$x_resident, d1$x_mutant[d1$w > 0]), ybar)

  plot(d2$x_mutant, d2$alpha, log="x", type="l", ylim=ylim, las=1,
       xlab="Height at maturation (m)",
       ylab="", yaxt="n")
  axis(2, labels=FALSE)
  points(d2$x_resident, 1.0, pch=19)
  abline(v=d2$x_resident, h=1.0, lty=2)
  black_bar(range(d2$x_resident, d2$x_mutant[d2$w > 0]), ybar)
}
