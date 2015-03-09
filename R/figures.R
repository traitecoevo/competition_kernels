fig_kernel <- function() {
  par(mar=c(2, 2, .5, .5))
  xlim <- c(-5, 5)

  x <- seq(xlim[1] * .9, xlim[2] * .9, length.out=5)
  xx <- seq(xlim[1], xlim[2], length.out=400)

  plot(NA, xlim=xlim, ylim=c(0, 1.1), xaxt="n", yaxt="n",
       xlab="Resource", ylab="Utilisation", mgp=c(.5, 0, 0))
  for (xi in x) {
    lines(xx, dnorm(xx, xi) / dnorm(0))
  }
  text(x, 1.05)
}

fig_components <- function() {
  mat <- rbind(c(0, 1, 1, 0),
               c(0, 2, 2, 0),
               c(3, 3, 4, 4),
               c(5, 5, 6, 6))
  xlim <- c(-2, 2)
  xx <- seq(xlim[[1]], xlim[[2]], length.out=200)
  p <- dd99_parameters()
  cex <- .66

  x1 <- -1.0
  x2 <- 0.0
 
  ## First, r (it's a straight line)
  r <- dd99_max_growth_rate(xx, p)
  K <- dd99_carrying_capacity(xx, p)
  N1 <- dd99_carrying_capacity(x1, p)
  N2 <- dd99_carrying_capacity(x2, p)
  w1 <- dd99_fitness(xx, x1, N1, p)
  w2 <- dd99_fitness(xx, x2, N2, p)
  a1 <- compute_alpha(w1, r, K, N1)
  a2 <- compute_alpha(w2, r, K, N2)

  layout(mat)
  par(mar=c(2, 2, .5, .5), oma=c(2, 2, 0, 0))
  
  plot(xx, r, las=1, ylim=c(0, 1.2), type="l")
  mtext("Max growth rate (r)", 2, xpd=NA, line=3, cex=cex)

  plot(xx, K, las=1, type="l", ylim=c(0, max(K) * 1.1))
  mtext("Carrying capacity (K)", 2, xpd=NA, line=3, cex=cex)

  ## TODO: Add 5% buffer.
  ylim_w <- c(min(w1, w2), 1.2)
  plot(xx, w1, las=1, type="l", ylim=ylim_w)
  mtext("Fitness (w)", 2, xpd=NA, line=3, cex=cex)
  points(x1, 0.0, pch=19)
  abline(h=0, v=x1, lty=2)
  
  plot(xx, w2, las=1, type="l", ylim=ylim_w)
  points(x2, 0.0, pch=19)
  abline(h=0, v=x2, lty=2)

  ylim_a <- c(0, 1)
  plot(xx, a1, las=1, type="l", ylim=ylim_a)
  mtext("Competition (a)", 2, xpd=NA, line=3, cex=cex)
  mtext("Ecological character (trait)", 1, xpd=NA, line=2, cex=cex)    
  points(x1, 1.0, pch=19)
  abline(h=1.0, v=x1, lty=2)
  
  plot(xx, a2, las=1, type="l", ylim=ylim_a)
  mtext("Ecological character (trait)", 1, xpd=NA, line=2, cex=cex)
  points(x2, 1.0, pch=19)
  abline(h=1.0, v=x2, lty=2)
}
