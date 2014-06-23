prep_tree_lma <- function() {
  tree:::run.cached(tree_competition("lma",  c(0.03, 3.8), 0.2),
                    "output/res_tree_lma.rds")
}

prep_tree_hmat <- function() {
  tree:::run.cached(tree_competition("hmat", c(10, 50), 18),
                    "output/res_tree_hmat.rds")
}

fig_tree <- function(res, plt, what) {
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
  ylim <- expand_range(range(res[[what]], 0, 1), down=FALSE)
  
  if (interactive()) {
    grid.newpage()
    popViewport(0)
  }

  pushViewport(viewport(w=w, h=h, x=x, y=y, xscale=xlim, yscale=ylim))
  cols <- colours_use()
  grid.rect(gp=gpar(col=NA, fill=cols$cream))

  # grid_xaxis_simple(c(0, 1), unit(-.5, "lines"))
  grid_xlab_simple(lab, unit(-1.2, "lines"))

  grid_ylab_simple(if (what == "alpha") "Competition" else "Fitness",
                   unit(-1.2, "lines"))
  if (what == "alpha") {
    if (res$trait == "hmat") {
      grid.text(format(c(0, 1)), unit(-.5, "lines"),
                unit(c(0, .8), "lines"), just="center")
    } else {
      grid_yaxis_simple(unit(-.5, "lines"), c(0, 1))
    }
    grid.rect(y=0, height=unit(1, "native"), just="bottom",
              gp=gpar(col=NA, fill=make_transparent("white", .5)))
  } else {
    grid.rect(y=0, height=unit(-ylim[1], "native"), just="bottom",
              gp=gpar(col=NA, fill=make_transparent("white", .5)))
  }
  
  grid.lines(unit(rep(log(res$x.resident), 2), "native"), c(0, 1),
             gp=gpar(lty=2, col=cols$orange, lwd=2, lineend="butt"))

  if (plt > 1) {
    grid.lines(unit(xx, "native"), unit(res[[what]], "native"),
               gp=gpar(lwd=2, col=cols$blue, lineend="butt"))
  }
  popViewport()

  pushViewport(viewport(w=w, h=mar[3], x=x, y=1, just="top", xscale=xlim))
  grid.picture(colour_picture(if (trait == "lma") leaf else tree,
                              cols$orange),
               x=unit(log(res$x.resident), "native"),
               y=unit(0.03, "npc"), height=unit(.6, "npc"), just="bottom")
  popViewport()
}

fig_tree_hmat <- function() {
  tree <- vector_read("pics/tree.svg")
  cols <- colours_use()
  pushViewport(viewport(width=unit(1, "snpc"), height=unit(1, "snpc"),
                        xscale=c(-1, 1), yscale=c(-1, 1)))
  on.exit(popViewport())
  grid.roundrect(gp=gpar(col=NA, fill=cols$cream))
  grid.picture(colour_picture(tree, cols$orange),
               .42, .12, height=.8, just="bottom")
  grid.picture(colour_picture(tree, cols$blue),
               .72, .08, height=.3, just="bottom")
}
