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
  ybar <- 0.05

  par(mfrow=c(2, 3), mar=rep(.5, 4), oma=c(2, 2, 0, 0))

  plot(NA, xlim=0:1, ylim=0:1, xlab="", ylab="", xaxt="n", yaxt="n")
  for (i in xx) {
    rstar_lines(rbind(i), dat1$p, col="grey")
  }
  rstar_lines(rbind(dat1$x1), dat1$p, col="red")
  rstar_lines(rbind(dat1$x2), dat1$p, col="blue")
  abline(0, 1, lty=2)

  plot(dat1$x, dat1$a1, type="l", ylim=ylim_alpha)
  abline(h=1.0, v=dat1$x1, lty=2, col="darkgrey")
  points(dat1$x1, 1.0, pch=19)
  black_bar(range(dat1$x[dat1$w1 > 0]), ybar)

  plot(dat1$x, dat1$a2, type="l", ylim=ylim_alpha)
  abline(h=1.0, v=dat1$x2, lty=2, col="darkgrey")
  points(dat1$x2, 1.0, pch=19)

  plot(NA, xlim=0:1, ylim=0:1, xlab="", ylab="", xaxt="n", yaxt="n")
  for (i in xx) {
    rstar_lines(rbind(i), dat2$p, col="grey")
  }
  rstar_lines(rbind(dat2$x1), dat2$p, col="red")
  rstar_lines(rbind(dat2$x2), dat2$p, col="blue")
  R <- drop(rstar_Rstar(rbind(dat2$x2), dat2$p))
  C <- dat2$p$C(rbind(dat2$x2))
  abcline(R[1], R[2], C[2] / C[1], lty=2)

  plot(dat2$x, dat2$a1, type="l", ylim=ylim_alpha)
  abline(h=1.0, v=dat2$x1, lty=2, col="darkgrey")
  points(dat2$x1, 1.0, pch=19)
  black_bar(range(dat2$x[dat2$w1 > 0]), ybar)

  plot(dat2$x, dat2$a2, type="l", ylim=ylim_alpha)
  abline(h=1.0, v=dat2$x2, lty=2, col="darkgrey")
  points(dat2$x2, 1.0, pch=19)
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
