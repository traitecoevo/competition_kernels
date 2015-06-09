dat_rstar <- function() {
  C1 <- c(.2, .2) * 10
  C2 <- c(.3, .7)
  x1 <- c(0.6, 0.5)
  x2 <- c(0.7, 0.9)
  p1 <- rstar_parameters(rstar_mat_2_tradeoff, matrix(C1, nrow=2), S=0.25)
  p2 <- rstar_parameters(rstar_mat_2_tradeoff, matrix(C2, nrow=2), S=0.5)
  list(dat_rstar1(p1, x1),
       dat_rstar1(p2, x2))
}

dat_rstar1 <- function(p, x_resident) {
  dat <- lapply(x_resident,
                function(xi) rstar_competition(matrix(xi, 1, 1), p))
  list(p=p,
       x=drop(dat[[1]]$x_invade),
       r=dat[[1]]$r_invade,
       K=dat[[1]]$K_invade,
       x1=x_resident[[1]],
       x2=x_resident[[2]],
       N1=dat[[1]]$N_resident,
       N2=dat[[2]]$N_resident,
       w1=dat[[1]]$w_invade,
       w2=dat[[2]]$w_invade,
       a1=dat[[1]]$alpha,
       a2=dat[[2]]$alpha)
}

fig_rstar <- function() {
  dat <- dat_rstar()
  dat1 <- dat[[1]]
  dat2 <- dat[[2]]

  xx <- seq(0, 1, length.out=6)
  ylim_alpha <- c(.45, 1.75)

  par(mfrow=c(2, 2), mar=rep(1, 4), oma=c(3, 3, 1, 1))

  plot(dat1$x, dat1$a1, type="l", ylim=ylim_alpha, las=1, xaxt="n")
  axis(1, labels=FALSE)
  abline(h=1.0, v=dat1$x1, lty=2, col="darkgrey")
  points(dat1$x1, 1.0, pch=19)
  add_black_bar(dat1$x, dat1$w1)
  label_panel(1)
  mtext(expression("Competition (" * alpha * ")"), 2, 1.6, outer=TRUE)
  mtext("Trait value", 1, 1.8, outer=TRUE)
  mtext("Away from attractor", 3, 0.5, xpd=NA)

  plot(dat1$x, dat1$a2, type="l", ylim=ylim_alpha, xaxt="n", yaxt="n")
  axis(1, labels=FALSE)
  axis(2, labels=FALSE)
  abline(h=1.0, v=dat1$x2, lty=2, col="darkgrey")
  points(dat1$x2, 1.0, pch=19)
  label_panel(2)
  mtext("At attractor", 3, 0.5, xpd=NA)
  mtext("Equal resource supply", 4, 0.5, xpd=NA)

  plot(dat2$x, dat2$a1, type="l", ylim=ylim_alpha, las=1)
  abline(h=1.0, v=dat2$x1, lty=2, col="darkgrey")
  points(dat2$x1, 1.0, pch=19)
  add_black_bar(dat2$x, dat2$w1)
  label_panel(3)

  plot(dat2$x, dat2$a2, type="l", ylim=ylim_alpha, yaxt="n")
  axis(2, labels=FALSE)
  abline(h=1.0, v=dat2$x2, lty=2, col="darkgrey")
  points(dat2$x2, 1.0, pch=19)
  label_panel(4)
  mtext("Unequal resource supply", 4, 0.5, xpd=NA)
}

fig_rstar_components <- function(i) {
  d <- dat_rstar()[[i]]
  plot_components(d$x, d$r, d$K, d$x1, d$x2, d$N1, d$N2,
                  d$w1, d$w2, d$a1, d$a2)
}

rstar_add_abrams <- function() {
  dat <- dat_rstar()
  lapply(dat, rstar_add_abrams1)
}

rstar_add_abrams1 <- function(d) {
  p <- d$p
  x_invade <- d$x

  f_dgi_dNj <- function(N_res, x_res) {
    rstar_fitness_given_N(rbind(x_invade), rbind(x_res), N_res, p)
  }
  dgi_dNj1 <- drop(numDeriv::jacobian(f_dgi_dNj, d$N1, x_res=d$x1))
  dgi_dNj2 <- drop(numDeriv::jacobian(f_dgi_dNj, d$N2, x_res=d$x2))

  f_dgi_dNi <- function(N_invade, N_res, x_res) {
    f <- function(x) {
      rstar_fitness_given_N(cbind(x),
                            cbind(x_res, x),
                            c(N_res, N_invade), p)
    }
    sapply(x_invade, f)
  }
  dN <- 1e-4
  dgi_dNi1 <- (f_dgi_dNi(dN, d$N1, d$x1) - d$w1) / dN
  dgi_dNi2 <- (f_dgi_dNi(dN, d$N2, d$x2) - d$w2) / dN

  d$X1 <- dgi_dNj1 / dgi_dNi1
  d$X2 <- dgi_dNj2 / dgi_dNi2

  d
}

fig_rstar_abrams <- function(d) {
  d1 <- d[[1]]
  d2 <- d[[2]]

  ylim <- c(.45, 1.75)

  par(mfrow=c(2, 2), mar=rep(.5, 4), oma=c(2, 2, 0, 0))

  plot(d1$x, d1$X1, type="l", ylim=ylim, las=1, xaxt="n")
  axis(1, labels=FALSE)
  lines(d1$x, d1$a1, lty=2)
  abline(v=d1$x1, h=1, lty=3)
  points(d1$x1, 1, pch=19)

  plot(d1$x, d1$X2, type="l", ylim=ylim, las=1, xaxt="n", yaxt="n")
  axis(1, labels=FALSE)
  axis(2, labels=FALSE)
  lines(d1$x, d1$a2, lty=2)
  abline(v=d1$x2, h=1, lty=3)
  points(d1$x2, 1, pch=19)

  plot(d2$x, d2$X1, type="l", ylim=ylim, las=1)
  lines(d2$x, d2$a1, lty=2)
  abline(v=d2$x1, h=1, lty=3)
  points(d2$x1, 1, pch=19)

  plot(d2$x, d2$X2, type="l", ylim=ylim, las=1, yaxt="n")
  axis(1, labels=FALSE)
  lines(d2$x, d2$a2, lty=2)
  abline(v=d2$x2, h=1, lty=3)
  points(d2$x2, 1, pch=19)
}

## R* density dependence
## I might keep this one as an additional supplementary figure?
fig_rstar_density_lines <- function() {
  d <- dat_rstar()
  d1 <- d[[1]]
  d2 <- d[[2]]

  ## Come up with a vector of densities:
  scal <- c(0.25, 0.5, 1, 2)

  d1_scal1 <- lapply(scal, function(s)
    rstar_competition(rbind(d1$x1), d1$p, d1$N1 * s))
  d1_scal2 <- lapply(scal, function(s)
    rstar_competition(rbind(d1$x2), d1$p, d1$N2 * s))
  d2_scal1 <- lapply(scal, function(s)
    rstar_competition(rbind(d2$x1), d2$p, d2$N1 * s))
  d2_scal2 <- lapply(scal, function(s)
    rstar_competition(rbind(d2$x2), d2$p, d2$N2 * s))

  a11 <- sapply(d1_scal1, "[[", "alpha")
  a12 <- sapply(d1_scal2, "[[", "alpha")
  a21 <- sapply(d2_scal1, "[[", "alpha")
  a22 <- sapply(d2_scal2, "[[", "alpha")

  x <- d1$x
  ylim <- c(.45, 1.75)
  col <- grey(seq(.8, 0, length.out=4))

  par(mfrow=c(2, 2), mar=rep(.5, 4), oma=c(2, 2, 0, 0))

  ## Hmm: competition should always be 1 for self, no?  Suggests I've
  ## got something a bit wrong?
  ##
  ## This suggests that per-capita competition *decreases* with
  ## increasing density (decelerating competition).
  matplot(x, a11, type="l", lty=1, col=col,
          ylim=ylim, las=1, xaxt="n")
  axis(1, labels=FALSE)
  abline(v=d1$x1, h=1, lty=3)
  points(d1$x1, 1, pch=19)

  matplot(x, a12, type="l", lty=1, col=col,
          ylim=ylim, las=1, xaxt="n", yaxt="n")
  axis(1, labels=FALSE)
  axis(2, labels=FALSE)
  abline(v=d1$x2, h=1, lty=3)
  points(d1$x2, 1, pch=19)

  matplot(x, a21, type="l", lty=1, col=col,
          ylim=ylim, las=1)
  abline(v=d2$x1, h=1, lty=3)
  points(d2$x1, 1, pch=19)

  matplot(x, a22, type="l", lty=1, col=col,
          ylim=ylim, las=1, yaxt="n")
  axis(2, labels=FALSE)
  abline(v=d2$x2, h=1, lty=3)
  points(d2$x2, 1, pch=19)
}

fig_rstar_density <- function(type) {
  d <- dat_rstar()[[type]]

  ## Come up with a vector of densities:
  scal <- seq_log(0.25, 2.0, 50)

  d_scal1 <- lapply(scal, function(s)
                    rstar_competition(rbind(d$x1), d$p, d$N1 * s))
  d_scal2 <- lapply(scal, function(s)
                    rstar_competition(rbind(d$x2), d$p, d$N2 * s))

  ## TODO: massive overlap here.
  a1 <- sapply(d_scal1, "[[", "alpha")
  a2 <- sapply(d_scal2, "[[", "alpha")

  i <- which.min(abs(scal - 1))
  n <- length(scal)
  b1 <- a1 / matrix(rep(a1[, i], n), ncol=n)
  b2 <- a2 / matrix(rep(a2[, i], n), ncol=n)

  cols <- c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7",
            "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC")
  cols2 <- colorRampPalette(rev(cols))(51)

  r <- range(a1, a2, na.rm=TRUE)
  zmax <- max(abs(r - 1))
  zlim <- c(-zmax, zmax) + 1

  r <- range(b1[i, ], b2[i, ], na.rm=TRUE)
  bmax <- max(abs(r - 1))
  blim <- c(-bmax, bmax) + 1

  x <- d$x

  par(mfrow=c(2, 2), mar=rep(1, 4), oma=c(3, 3, 0, 0))
  image(x, scal, a1, log="y", col=cols2, zlim=zlim, xaxt="n", las=1)
  axis(1, labels=FALSE)
  contour(x, scal, a1, add=TRUE, levels=1, drawlabels=FALSE)
  abline(h=1, v=d$x1, lty=3)
  points(d$x1, 1, pch=19)
  label_panel(1)

  image(x, scal, a2, log="y", col=cols2, zlim=zlim, xaxt="n", yaxt="n")
  axis(1, labels=FALSE)
  axis(2, labels=FALSE)
  contour(x, scal, a2, add=TRUE, levels=1, drawlabels=FALSE)
  abline(h=1, v=d$x2, lty=3)
  points(d$x2, 1, pch=19)
  label_panel(2)

  image(x, scal, b1, log="y", col=cols2, zlim=zlim, las=1)
  abline(h=1, v=d$x1, lty=3)
  points(d$x1, 1, pch=19)
  label_panel(3)

  image(x, scal, b2, log="y", col=cols2, zlim=zlim, yaxt="n")
  axis(2, labels=FALSE)
  abline(h=1, v=d$x2, lty=3)
  points(d$x2, 1, pch=19)
  label_panel(4)

  mtext("Relative population size", 2, 1.8, outer=TRUE)
  mtext("Trait value", 1, 1.8, outer=TRUE)
}
