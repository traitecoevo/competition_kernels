gaussian <- function(x, mean=0, sd=1, scale=1) {
  dnorm(x, mean, sd) / dnorm(0, sd=sd) * scale
}

## Finches from here
## http://www.science.ca/images/schluter-infographic.jpg

fig_gaussian_competition <- function(plt) {
  leaf <- vector_read("pics/leaf.svg")

  mar <- c(.15, .1, .2, .02)

  h <- 1 - sum(mar[c(1, 3)])
  w <- 1 - sum(mar[c(2, 4)])
  x <- mar[2] + w/2
  y <- mar[1] + h/2
  r <- c(-7, 7)

  p <- .1
  xl <- seq(p, 1-p, length=5)
  xl.native <- r[1] + diff(r) * xl
  hl <- seq_sqrt(.3, .9, length.out=5)

  cols <- colours()

  if (interactive()) {
    grid.newpage()
    popViewport(0)
  }

  pushViewport(viewport(w=w, h=h, x=x, y=y))
  grid.rect(gp=gpar(col=NA, fill=cols$bg))
  if (plt == 1) {
    grid_xlab_simple("Resource", unit(-.75, "lines"))
    grid_ylab_simple("Preference", unit(-1, "lines"))
  } else {
    grid_xlab_simple("Species trait", unit(-1, "lines"))
    grid_ylab_simple("Competition", unit(-1, "lines"))
  }
  pushViewport(viewport(xscale=r, clip="on"))
  xx <- seq(r[1], r[2], length=201)
  if (plt > 2) {
    for (i in seq_len(5)[-3]) {
      grid.lines(xx, gaussian(xx, xl.native[i], scale=.9),
                 gp=gpar(col=cols$fg, lwd=2, lineend="butt"),
                 default.unit="native")
    }
  }
  grid.lines(xx, gaussian(xx, xl.native[3], scale=.9),
             gp=gpar(col=cols$hl, lwd=4, lineend="butt"),
             default.unit="native")
  popViewport(2)

  ## Place to put the leaf pictures:
  pushViewport(viewport(w=w, h=mar[3], x=x, y=1, xscale=r, just="top"))
  if (plt > 2) {
    for (i in seq_along(xl)[-3]) {
      grid.picture(colour_picture(leaf, cols$fg), x=unit(xl[i], "npc"),
                   y=unit(0, "npc"), just="bottom", height=unit(hl[i], "npc"))
    }
  }

  grid.picture(colour_picture(leaf, cols$hl), x=unit(xl.native[3], "native"),
               y=unit(0, "npc"), height=unit(.5, "npc"), just="bottom")
  popViewport()
}

## Need to have a figure so that I can handwave about bridging models
## and data.
fig_bridge <- function(plt) {
  cols <- list(brown="#574134",
               brown_dk="#473629",
               orange="#d54b1a",
               orange_dk="#bf4317",
               yellow="#e3a72f",
               yellow_dk="#c6891d",
               blue="#058789",
               blue_dk="#04797a",
               cream="#f0ecc9",
               cream_dk="#c7c094")

  trees <- vector_read("pics/trees.svg")

  h <- 0.3
  w <- 0.3

  mar <- .05
  x <- c(w/2 + mar, 1 - w/2 - mar)
  y <- c(1 - h/2 - mar, h/2 + mar)
  
  r <- c(-7, 7)
  p <- .1
  xl <- seq(p, 1-p, length=5)
  xl.native <- r[1] + diff(r) * xl
  
  if (interactive()) {
    grid.newpage()
    popViewport(0)
  }

  pushViewport(viewport(x[1], y[1], w, h))
  grid.roundrect(gp=gpar(col=NA, fill=cols$cream))
  pushViewport(viewport(h=.8, w=.8, xscale=r))
  xx <- seq(r[1], r[2], length=201)
  for (i in c(seq_len(5)[-3], 3)) {
    col <- if (i == 3) cols$orange else cols$blue
    lwd <- if (i == 3) 2.5 else 1.5
    grid.lines(xx, gaussian(xx, xl.native[i]),
               gp=gpar(col=col, lwd=lwd, lineend="butt"),
               default.unit="native")
  }
  popViewport(2)

  pushViewport(viewport(last(x), last(y), w, h))
  grid.roundrect(gp=gpar(col=NA, fill=cols$cream))
  grid.picture(colour_picture(trees, cols$blue),
               width=.8)
  popViewport()

  x.arr <- c(x[1] + w/2 + mar/2, last(x) - w/2 - mar/2)
  y.arr <- c(y[1] - h/2 - mar/2, last(y) + h/2 + mar/2)
  grid.lines(x.arr, y.arr, arrow=arrow(type="closed", angle=15),
             gp=gpar(col=cols$brown, fill=cols$brown, lwd=2))
  
  grid.text("?", y=.6, gp=gpar(cex=4))
}
