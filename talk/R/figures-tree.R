prep_tree_lma <- function() {
  tree:::run.cached(tree_competition("lma",  c(0.03, 3.8), 0.2),
                    "res_tree_lma.rds")
}

prep_tree_hmat <- function() {
  tree:::run.cached(tree_competition("hmat", c(10, 50), 16),
                    "res_tree_hmat.rds")
}

fig_tree_alpha <- function(res, plt) {
  trait <- res$trait
  lab <- c(hmat="Height at maturity",
           lma="Leaf mass per area")[[trait]]
  
  leaf <- vector_read("pics/leaf.svg")
  tree <- vector_read("pics/tree.svg")

  mar <- c(.15, .1, .2, .02)

  h <- 1 - sum(mar[c(1, 3)])
  w <- 1 - sum(mar[c(2, 4)])
  x <- mar[2] + w/2
  y <- mar[1] + h/2

  xx <- log(res$x)

  xlim <- range(xx)
  ylim <- expand_range(range(res$alpha, 0, 1), down=FALSE)
  
  if (interactive()) {
    grid.newpage()
    popViewport(0)
  }

  pushViewport(viewport(w=w, h=h, x=x, y=y, xscale=xlim, yscale=ylim))
  cols <- colours()
  grid.rect(gp=gpar(col=NA, fill=cols$bg))

  # grid_xaxis_simple(c(0, 1), unit(-.5, "lines"))
  grid_xlab_simple(lab, unit(-1.2, "lines"))

  grid_yaxis_simple(unit(-.5, "lines"), c(0, 1))
  grid_ylab_simple("Competition", unit(-1.2, "lines"))
  ## Work out a nicer way of doing this.
  grid.rect(y=0, height=unit(1, "native"), just="bottom",
            gp=gpar(col=NA, fill=make_transparent("white", .5)))
  grid.lines(unit(rep(log(res$x.resident), 2), "native"), c(0, 1),
             gp=gpar(lty=2, col=cols$hl, lwd=2))

  if (plt > 1) {
    grid.lines(unit(xx, "native"), unit(res$alpha, "native"),
               gp=gpar(lwd=2, col=cols$fg, lineend="butt"))
  }
  popViewport()

  pushViewport(viewport(w=w, h=mar[3], x=x, y=1, just="top", xscale=xlim))
  grid.picture(colour_picture(leaf, cols$hl),
               x=unit(log(res$x.resident), "native"),
               y=unit(0, "npc"), height=unit(.5, "npc"), just="bottom")
  popViewport()
}
