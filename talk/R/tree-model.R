fig_tree_model <- function(draft=FALSE, embed=FALSE) {
  cleanup <- function(x) {
    x[is.na(x)] <- 0
    x
  }
  simple.frame <- function(...) {
    grid.segments(c(0, 0), c(0, 0),
                  c(0, 1), c(1, 0), default.units="npc", ...)
  }  

  tree <- vector_read("pics/tree.svg")

  set.seed(1)
  xy <- pack.circles(1000, r0=0.03, w=1.2, max.size=.15)

  res <- readRDS("output/tree.output.rds")
  plant <- make.reference()
  assign("p.eta", 5, environment(plant$leaf.pdf))
  env <- res$light.env[[92]]
  h.max <- max(env[,1])
  height <- seq(0, h.max, length=201)
  light.env <- splinefun(env[,1], env[,2])
  assim <- cleanup(plant$leaf.pdf(height, h.max))
  cols <- colours_use()

  r <- 2/3

  pushViewport(viewport(width=unit(1, "snpc"), height=unit(1, "snpc"),
                        xscale=c(-1, 1), yscale=c(-1, 1)))
  x0 <- unit(0.03, "npc")
  y0 <- unit(0.5, "npc")
  w0 <- unit(0.33, "npc")
  pushViewport(viewport(x0, y0, w0, w0, just="left"))
  grid.roundrect(gp=gpar(col=NA, fill=cols$cream))
  pushViewport(viewport(xscale=c(0, 1.2)))
  cols_trees <- sample(c(cols$blue_dk, cols$blue),
                       ncol(xy), replace=TRUE)
  if (draft)
    grid.circle(unit(xy["x",], "native"), xy["y",], xy["r",])
  else
  for (i in seq_len(ncol(xy)))
    grid.picture(colour_picture(tree, cols_trees[i]),
                 unit(xy["x",i], "native"), xy["y",i],
                 just=c("centre", "centre"), height=xy["r",i] * 2)
  popViewport(2)

  dt <- 2*pi / 6
  grid.points(r * cos(dt), r*sin(dt))

  x1 <- x2 <- unit(r * cos(dt), "native") + unit(0.1, "npc")
  y1 <- unit(r * sin(dt), "native")
  y2 <- unit(r * sin(-dt), "native")

  w1 <- w2 <- 0.33

  pushViewport(viewport(x1, y1, w1, w1))
  grid.roundrect(gp=gpar(col=NA, fill=cols$cream))
  mar <- unit(1.2, "lines")
  pushViewport(viewport(mar, mar,
                        unit(1, "npc") - mar * 1.5,
                        unit(1, "npc") - mar * 1.5,
                        just=c("left", "bottom"),
                        xscale=c(0, 1),
                        yscale=range(env[,1])))
  simple.frame()
  grid.lines(env[,2], env[,1], default.units="native",
             gp=gpar(lwd=3, lineend="butt", col=cols$blue))
  if (!embed) {
    grid.text("Height", x=unit(-0.66, "lines"), rot=90)
    grid.text("Light", y=unit(-0.66, "lines"))
  }
  popViewport(2)

  pushViewport(viewport(x2, y2, w2, w2))
  grid.roundrect(gp=gpar(col=NA, fill=cols$cream))

  pushViewport(viewport(width=0.9))
  pushViewport(viewport(0, width=2/3, just="left"))
  grid.picture(colour_picture(tree, cols$blue), width=unit(1, "npc"))
  popViewport()

  pushViewport(viewport(1, width=1/3, just="right", height=0.76,
                        xscale=range(assim), yscale=range(height)))
  grid.polygon(assim, height,
               gp=gpar(fill=cols$blue, col=NA),
               default.units="native")
  ## grid.lines(assim, height,
  ##            gp=gpar(lwd=if (embed) 1.5 else 3, lineend="butt",
  ##              col=cols$blue),
  ##            default.units="native")
  simple.frame()
  popViewport(3)

  dt <- pi * 0.03
  theta0 <- c(pi - pi / 6,   pi * 0.12, 2*pi - pi * 0.4) - dt
  theta1 <- c(pi * 0.4,    - pi * 0.12, pi + pi / 6) + dt

  theta <- lapply(1:3, function(i)
                  seq(theta0[i], theta1[i], length=101))
  da <- unit(.1, "npc")
  arr <- arrow(type="closed", angle=20,
               length=convertWidth(da * 0.5, "cm"))
  lapply(theta, function(t)
         grid.lines(unit(r * cos(t), "native"),
                    unit(r * sin(t), "native"),
                    arrow=arr, gp=gpar(col="black", lwd=2,
                                 fill="black")))
  if (!embed) {
    theta2 <- (theta0 + theta1) / 2
    grid.text("Growth,\ndeaths",
              unit(1.1*r*cos(theta2[3]), "native"),
              unit(1.1*r*sin(theta2[3]), "native"),
              gp=gpar(col="black"), just=c("right", "top"))
    grid.text("Shading",
              unit(1.1*r*cos(theta2[1]), "native"),
              unit(1.1*r*sin(theta2[1]), "native"),
              gp=gpar(col="black"), just=c("right", "bottom"))
    grid.text("Photosynthesis", unit(r * 0.9, "native"),
              just="right", gp=gpar(col="black"))
  }
  popViewport(1)
}
