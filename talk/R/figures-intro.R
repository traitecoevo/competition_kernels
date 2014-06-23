gaussian <- function(x, mean=0, sd=1, scale=1) {
  dnorm(x, mean, sd) / dnorm(0, sd=sd) * scale
}

## Finches from here
## http://www.science.ca/images/schluter-infographic.jpg

fig_gaussian_competition <- function(plt) {
  leaf <- vector_read("pics/leaf.svg")
  finch <- png::readPNG("pics/finch.png")

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

  cols <- colours_use()

  if (interactive()) {
    grid.newpage()
    popViewport(0)
  }

  pushViewport(viewport(w=w, h=h, x=x, y=y))
  grid.rect(gp=gpar(col=NA, fill=cols$cream))
  if (plt %in% c(1, 4)) {
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
                 gp=gpar(col=cols$blue, lwd=2, lineend="butt"),
                 default.unit="native")
    }
  }
  grid.lines(xx, gaussian(xx, xl.native[3], scale=.9),
             gp=gpar(col=cols$orange, lwd=4, lineend="butt"),
             default.unit="native")
  popViewport(2)

  ## Place to put the leaf pictures:
  pushViewport(viewport(w=w, h=mar[3], x=x, y=1, xscale=r, just="top"))
  if (plt > 2) {
    for (i in seq_along(xl)[-3]) {
      if (plt >= 5) {
        grid.picture(colour_picture(leaf, cols$blue),
                     x=unit(xl[i], "npc"),
                     y=unit(0, "npc"), just="bottom",
                     height=unit(hl[i], "npc"))
      } else {
        grid.raster(mask_raster(finch, cols$blue),
                    x=unit(xl[i], "npc"),
                    y=unit(0, "npc"), just="bottom",
                    height=unit(hl[i], "npc"))
      }
    }
  }

  if (plt >= 5) {
    grid.picture(colour_picture(leaf, cols$orange),
                 x=unit(xl.native[3], "native"), y=unit(0, "npc"),
                 height=unit(.5, "npc"), just="bottom")
  } else {
    grid.raster(mask_raster(finch, cols$orange),
                x=unit(xl.native[3], "native"), y=unit(0, "npc"),
                height=unit(.5, "npc"), just="bottom")
  }
  popViewport()
}

## Need to have a figure so that I can handwave about bridging models
## and data.
fig_bridge <- function(plt) {
  cols <- colours_use()

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
  grid.lines(x.arr, y.arr, arrow=arrow(type="closed", angle=15, ends="both"),
             gp=gpar(col="black", fill="black", lwd=2))
}
