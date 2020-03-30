fig_kernel <- function() {
  par(mar=c(2, 2, .5, .5))
  xlim <- c(-5, 5)

  x <- seq(xlim[1] * .9, xlim[2] * .9, length.out=5)
  xx <- seq(xlim[1], xlim[2], length.out=400)

  plot(NA, xlim=xlim, ylim=c(0, 1.1), xaxt="n", yaxt="n",
       xlab="Resource continuum", ylab="Utilisation", mgp=c(.5, 0, 0))
  for (xi in x) {
    lines(xx, dnorm(xx, xi) / dnorm(0))
  }
  text(x, 1.05)
}

fig_components <- function() {
  p <- dd99_parameters()
  xx <- seq(-2, 2, length.out=200)
  x1 <- -1.0
  x2 <- 0.0

  r <- dd99_max_growth_rate(xx, p)
  K <- dd99_carrying_capacity(xx, p)
  N1 <- dd99_carrying_capacity(x1, p)
  N2 <- dd99_carrying_capacity(x2, p)
  w1 <- dd99_fitness(xx, x1, N1, p)
  w2 <- dd99_fitness(xx, x2, N2, p)
  a1 <- compute_alpha(w1, r, K, N1)
  a2 <- compute_alpha(w2, r, K, N2)

  plot_components(xx, r, K, x1, x2, N1, N2, w1, w2, a1, a2)
}

plot_components <- function(x, r, K, x1, x2, N1, N2, w1, w2, a1, a2,
                            xlab="Ecological character (trait)",
                            ylim_r=NULL, ylim_K=NULL, ylim_w=NULL,
                            ylim_a=NULL, log="") {
  cex <- .8

  if (is.null(ylim_r)) {
    ylim_r <- expand_range(c(0, max(r)), c(0, 0.1))
  }
  if (is.null(ylim_K)) {
    ylim_K <- expand_range(c(0, max(K)), c(0, 0.1))
  }
  if (is.null(ylim_w)) {
    ylim_w <- expand_range(range(w1, w2), 0.1)
  }
  if (is.null(ylim_a)) {
    ylim_a <- expand_range(range(a1, a2, na.rm=TRUE), 0.1)
  }

  ylab_r <- expression("Maximum growth rate (" * italic(r) * ")")
  ylab_K <- expression("Carrying capacity (" * italic(K) * ")")
  ylab_w <- expression("Fitness (" * italic(w) * ")")
  ylab_a <- expression("Competition (" * alpha * ")")

  par(mfrow=c(2, 3), mar=c(2.5, 4.0, .5, .5), oma=c(2, 1, 2, 0),
      mgp=c(2.5, 1, 0))

  plot(x, r, las=1, ylim=ylim_r, type="l", log=log, ylab="")
  mtext(ylab_r, 2, 3, cex=0.8*cex)
  label_panel(1)

  plot(x, w1, las=1, type="l", ylim=ylim_w, log=log, ylab="")
  mtext(ylab_w, 2, 3, cex=0.8*cex)
  abline(h=0, v=x1, lty=2, col="grey")
  points(x1, 0.0, pch=19)
  label_panel(3)

  mtext("Away from attractor", 3, 1, xpd=NA, cex=cex)


  plot(x, w2, las=1, type="l", ylim=ylim_w, log=log, ylab="")#, yaxt="n")
#  axis(2, labels=FALSE)
  abline(h=0, v=x2, lty=2, col="grey")
  points(x2, 0.0, pch=19)
  label_panel(5)

  mtext("At attractor", 3, 1, xpd=NA, cex=.8)


  plot(x, K, las=1, type="l", ylim=ylim_K, log=log, ylab="")
  mtext(ylab_K, 2, 3, cex=0.8*cex)
  label_panel(2)

  plot(x, a1, las=1, type="l", ylim=ylim_a, log=log, ylab="")
  mtext(ylab_a, 2, 3, cex=0.8*cex)
  abline(h=1.0, v=x1, lty=2, col="grey")
  points(x1, 1.0, pch=19)
  mtext(xlab, 1, xpd=NA, line=3, cex=cex)
  label_panel(4)
  add_black_bar(x, w1)

  plot(x, a2, las=1, type="l", ylim=ylim_a, log=log, ylab="")
  axis(2, labels=FALSE)
  abline(h=1.0, v=x2, lty=2, col="grey")
  points(x2, 1.0, pch=19)
#  mtext(xlab, 1, xpd=NA, line=2.3, cex=cex)
  label_panel(6)
  add_black_bar(x, w2)
}

add_black_bar <- function(x, w) {
  usr <- par("usr")
  ## Identify changes in sign in w, counting 0 as non-invasible:
  sw <- sign(w)
  sw[sw == 0] <- -1
  i <- which(diff(sw) != 0)
  xx <- (x[i] + x[i + 1]) / 2
  if (sw[1] > 0) {
    xx <- c(usr[1], xx)
  }
  if (sw[length(sw)] > 0) {
    xx <- c(xx, usr[2])
  }
  xx <- matrix(xx, 2)
  for (i in seq_len(ncol(xx))) {
    black_bar(xx[, i])
  }
}
