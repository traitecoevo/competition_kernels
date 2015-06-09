plant_trait_lab <- function(trait) {
  switch(trait,
         lma="Leaf mass per area (lma; kg / m^2)",
         hmat="Height at maturation (m)",
         stop("Unknown trait ", trait))
}

fig_plant <- function(d1, d2) {
  if (d1$trait == "lma") {
    ylim1 <- ylim2 <- expand_range(range(range(d1$alpha, d2$alpha)), 0.05)
    xlim <- exp(expand_range(log(range(d1$x_mutant)), c(0, 0)))
    xlab <- "Leaf mass per area (lma; kg / m^2)"
    xlab <- expression("Leaf mass per area (lma; kg / " ~  m ^ 2 * ")")
    xlab_line <- 2.4
  } else {
    ylim1 <- expand_range(range(range(d1$alpha)), 0.05)
    ylim2 <- expand_range(range(range(d2$alpha)), 0.05)
    xlim <- exp(expand_range(log(range(d1$x_mutant)), c(0, 0)))
    xlab <- "Height at maturation (m)"
    xlab_line <- 2.1
  }

  ylab <- expression("Competition (" * alpha * ")")
  par(mfrow=c(2, 1), mar=c(1, 3.7, .5, .5), oma=c(2.3, 0, 0, 0),
      mgp=c(2.3, 1, 0))
  plot(d1$x_mutant, d1$alpha, log="x", type="l",
       xlim=xlim, ylim=ylim1, las=1, xaxt="n", xlab="", ylab=ylab)
  axis(1, labels=FALSE)
  abline(v=d1$x_resident, h=1.0, lty=2, col="grey")
  points(d1$x_resident, 1.0, pch=19)
  add_black_bar(d1$x_mutant, d1$w)
  label_panel(1)

  plot(d2$x_mutant, d2$alpha, log="x", type="l",
       xlim=xlim, ylim=ylim2, las=1, ylab=ylab)
  abline(v=d2$x_resident, h=1.0, lty=2, col="grey")
  points(d2$x_resident, 1.0, pch=19)
  add_black_bar(d2$x_mutant, d2$w)
  mtext(xlab, 1, xlab_line, xpd=NA)
  label_panel(2)
}

fig_plant_components <- function(d1, d2) {
  xlab <- plant_trait_lab(d1$trait)
  x <- d1$x_mutant
  if (d1$trait == "hmat") {
    ylim_a <- c(-1, 40)
    ylim_w <- c(-5, 3)
  } else {
    ylim_a <- c(0, 2.5)
    ylim_w <- c(-8, 8)
  }

  plot_components(x, d1$r(x), d1$K(x),
                  d1$x_resident, d2$x_resident,
                  d1$N_resident, d2$N_resident,
                  d1$w, d2$w, d1$alpha, d2$alpha,
                  log="x", xlab=xlab, ylim_a=ylim_a, ylim_w=ylim_w)
}

fig_plant_density <- function(d1, d2) {
  a1 <- sapply(d1$dat, "[[", "alpha")
  a2 <- sapply(d2$dat, "[[", "alpha")

  i <- which.min(abs(d1$scal - 1))
  n <- length(d1$scal)
  b1 <- a1 / matrix(rep(a1[, i], n), ncol=n)
  b2 <- a2 / matrix(rep(a2[, i], n), ncol=n)

  ## Need to drop the edges as they're extreme (inviable most of the
  ## time):
  i <- -c(1, nrow(a1))
  r <- range(a1[i, ], a2[i, ], na.rm=TRUE)
  zmax <- max(abs(r - 1))
  zlim <- c(-zmax, zmax) + 1

  r <- range(b1[i, ], b2[i, ], na.rm=TRUE)
  bmax <- max(abs(r - 1))
  blim <- c(-bmax, bmax) + 1

  cols <- c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7",
            "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC")
  cols2 <- colorRampPalette(rev(cols))(51)

  par(mfrow=c(2, 2), mar=rep(1, 4), oma=c(3, 3, 0, 0))
  image(d1$x[i], d1$scal, a1[i, ], log="xy", col=cols2, zlim=zlim,
        las=1, xaxt="n")
  axis(1, labels=FALSE)
  contour(d1$x[i], d1$scal, a1[i, ], add=TRUE, levels=1, drawlabels=FALSE)
  abline(h=1, v=d1$dat[[1]]$x_resident, lty=3)
  points(d1$dat[[1]]$x_resident, 1, pch=19)
  label_panel(1)

  image(d2$x[i], d2$scal, a2[i, ], log="xy", col=cols2, zlim=zlim,
        las=1, xaxt="n", yaxt="n")
  axis(1, labels=FALSE)
  axis(2, labels=FALSE)
  contour(d2$x[i], d2$scal, a2[i, ], add=TRUE, levels=1, drawlabels=FALSE)
  abline(h=1, v=d2$dat[[1]]$x_resident, lty=3)
  points(d2$dat[[1]]$x_resident, 1, pch=19)
  label_panel(2)

  image(d1$x[i], d1$scal, b1[i, ], log="xy", col=cols2, zlim=blim, las=1)
  abline(h=1, v=d1$dat[[1]]$x_resident, lty=3)
  points(d1$dat[[1]]$x_resident, 1, pch=19)
  label_panel(3)

  image(d2$x[i], d2$scal, b2[i, ], log="xy", col=cols2, zlim=blim,
        las=1, yaxt="n")
  axis(2, labels=FALSE)
  abline(h=1, v=d2$dat[[1]]$x_resident, lty=3)
  points(d2$dat[[1]]$x_resident, 1, pch=19)
  label_panel(4)

  mtext("Relative population size", 2, 1.8, outer=TRUE)
  mtext(plant_trait_lab(d1$dat[[1]]$trait), 1, 1.8, outer=TRUE)
}
