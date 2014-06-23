library(Revolve)

compute_alpha <- function(w.i, r.i, K.i, N.r) {
  (1 - 1 / r.i * w.i) * K.i / N.r
}

## Rstar, 1 resource case.
prep_rstar_1 <- function() {
  mat_r1 <- rstar_matrices(rstar_mat_1,
                           make_rstar_mat_constant(matrix(.5)))
  m_r1 <- rstar(mat_r1, S=1)
  x.r <- 0.5
  eq_r1 <- m_r1$single_equilibrium(matrix(x.r))
  x.mutant <- seq(0, 1, length=101)[-1]
  x.K <- rbind(x.mutant, deparse.level=0)

  r.mutant <- m_r1$max_growth_rate(x.K)
  K.mutant <- m_r1$carrying_capacity(x.K)
  N.resident <- eq_r1$y

  w.mutant <- m_r1$fitness(x.K, eq_r1$x, eq_r1$y, eq_r1$R)
  alpha <- compute_alpha(w.mutant, r.mutant, K.mutant, N.resident)

  list(x=x.mutant,
       r=r.mutant,
       K=K.mutant,
       w=w.mutant,
       alpha=alpha,
       x.resident=x.r,
       N.resident=N.resident,
       n.resource=1)
}

prep_rstar_2 <- function(x.r=.3, S=c(.7, .5), C=c(.3, .7)) {
  mat_r2 <- rstar_matrices(rstar_mat_2_tradeoff,
                           make_rstar_mat_constant(matrix(C, nrow=2)))
  m_r2 <- rstar(mat_r2, S=S)
  eq_r2 <- m_r2$single_equilibrium(matrix(x.r))
  x.mutant <- seq(0, 1, length=301)
  x.K <- rbind(x.mutant, deparse.level=0)

  r.mutant <- m_r2$max_growth_rate(x.K)
  K.mutant <- m_r2$carrying_capacity(x.K) 
  N.resident <- eq_r2$y

  w.mutant <- m_r2$fitness(x.K, eq_r2$x, eq_r2$y, eq_r2$R)
  alpha <- compute_alpha(w.mutant, r.mutant, K.mutant, N.resident)

  list(x=x.mutant,
       r=r.mutant,
       K=K.mutant,
       w=w.mutant,       
       alpha=alpha,
       x.resident=x.r,
       N.resident=N.resident,
       n.resource=2)
}

fig_rstar_alpha <- function(res, plt) {
  leaf <- vector_read("pics/leaf.svg")
  ## So draw this as "resource required for growth"
  mar <- c(.18, .1, .18, .02)

  h <- 1 - sum(mar[c(1, 3)])
  w <- 1 - sum(mar[c(2, 4)])
  x <- mar[2] + w/2
  y <- mar[1] + h/2

  ylim <- expand_range(range(res$alpha, 0, 1), down=FALSE)
  
  if (interactive()) {
    grid.newpage()
    popViewport(0)
  }

  pushViewport(viewport(w=w, h=h, x=x, y=y, yscale=ylim))
  cols <- colours_use()
  grid.rect(gp=gpar(col=NA, fill=cols$cream))

  at <- unit(-.5, "lines")
  if (res$n.resource == 1) {
    grid_xaxis_simple(c(0, 1), at)
  } else {
    grid.text("Specialist (A)", 0, at, just=c("left", "top"))
    grid.text("Generalist", .5, at, just=c("center", "top"))
    grid.text("Specialist (B)", 1, at, just=c("right", "top"))
  }
  grid_xlab_simple("Resource requirement",
                   unit(-1.2 - (res$n.resource == 2) * .5, "lines"))

  if (plt == 1) {
    y <- 0.4
    len <- if (res$n.resource == 1) 0.25 else 0.2
    dx <- 0.05
    dy <- 0.075
    arr <- arrow(type="closed", angle=15)
    arr.gp <- gpar(col="black", fill="black", lwd=2)

    if (res$n.resource == 1) {
      grid.lines(res$x.r - (dx + c(0, len)), rep(y, 2),
                 arrow=arr, gp=arr.gp)
      grid.lines(res$x.r + (dx + c(0, len)), rep(y, 2),
                 arrow=arr, gp=arr.gp)
      grid.text("Can invade", res$x.r - dx, y + dy,
                just=c("right", "bottom"))
      grid.text("Cannot invade", res$x.r + dx, y + dy,
                just=c("left", "bottom"))
    } else {
      grid.lines(res$x.r - (dx + c(0, len)), rep(y, 2),
                 arrow=arr, gp=arr.gp)
      grid.lines(res$x.r + (dx + c(0, len*2)), rep(y, 2),
                 arrow=arr, gp=arr.gp)
      grid.text("Cannot\ninvade", res$x.r - dx, y + dy,
                just=c("right", "bottom"))
      grid.text("Can invade", res$x.r + dx, y + dy,
                just=c("left", "bottom"))
    }
  }

  if (plt > 1) {
    grid_yaxis_simple(unit(-.5, "lines"), c(0, 1))
    grid_ylab_simple("Competition", unit(-1.2, "lines"))
    ## Work out a nicer way of doing this.
    grid.rect(y=0, height=unit(1, "native"), just="bottom",
              gp=gpar(col=NA, fill=make_transparent("white", .5)))
  }

  ## Add the resident:
  grid.lines(rep(res$x.r, 2), c(0, 1),
             gp=gpar(lty=2, col=cols$orange, lwd=2, lineend="butt"))
  
  if (plt > 2) {
    grid.lines(res$x, unit(res$alpha, "native"),
               gp=gpar(lwd=2, col=cols$blue, lineend="butt"))
  }
  if (plt > 3) {
    grid.lines(res$x,
               unit(gaussian(res$x, mean=res$x.r, sd=1/8), "native"),
               gp=gpar(col=cols$yellow_dk, lwd=1.5))
  }
  popViewport()

  pushViewport(viewport(w=w, h=mar[3], x=x, y=1, just="top"))
  grid.picture(colour_picture(leaf, cols$orange), x=res$x.r,
               y=unit(0, "npc"), height=unit(.5, "npc"), just="bottom")
  popViewport()
}
