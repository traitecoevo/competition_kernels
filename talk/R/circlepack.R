## Simple code for packing circles.
##
## Based on the algorithm here:
##   http://nodebox.net/code/index.php/shared_2008-08-07-12-55-33
##
## Because of how I intend to use this, there is also support for
## avoiding arbitrary rectangles, using the algorithm here:
##   http://stackoverflow.com/a/1879223/1798863
## However, this is optional.
pack.circles <- function(iterations, w=1, h=1, 
                         rects=NULL, nodes=NULL,
                         r0=NULL, g0=NULL, min.sep=NULL, max.size=NULL) {
  wh <- min(w, h)
  if (is.null(r0))
    r0 <- .01 * wh
  if (is.null(g0))
    g0 <- c(.002, .005) * wh
  if (is.null(min.sep))
    min.sep <- min(g0)
  if (is.null(max.size))
    max.size <- .1 * wh
  
  if (length(g0) != 2)
    stop("Expected g0 of length 2")
  random.circle <- make.random.circle(w, h, r0, g0)

  for (i in seq_len(iterations)) {
    nodes <- add.node(nodes, random.circle())
    nodes <- test.node.collision(nodes, min.sep)
    nodes <- test.bounding.box(nodes, w, h, min.sep)
    if (!is.null(rects))
      nodes <- test.rects(nodes, rects, min.sep)
    nodes <- check.new.node(nodes)
    nodes <- grow.nodes(nodes, max.size)
  }

  nodes
}

## Generating function for random circles
make.random.circle <- function(w, h, r0, g0) {
  circ <- function(x, y, r, g)
    c(x=x, y=y, r=r, alive=TRUE, growth=g)
  function()
    circ(runif(1, r0, w-r0), runif(1, r0, h-r0), r0,
         runif(1, g0[1], g0[2]))
}

add.node <- function(nodes, node) {
  cbind(nodes, node)
}

## See which nodes ahve bumped into one another.  Only actually check
## nodes that are still growing.
test.node.collision <- function(nodes, min.sep) {
  is.alive <- as.logical(nodes["alive",])
  alive <- nodes[,is.alive,drop=FALSE]
  f <- function(a, b) square(a - b)
  g <- function(a, b) square(a + b + min.sep)
  dx <- outer(alive["x",], nodes["x",], f)
  dy <- outer(alive["y",], nodes["y",], f)
  rr <- outer(alive["r",], nodes["r",], g)
  kill <- which(is.alive)[rowSums(dx + dy < rr) > 1]
  nodes["alive", kill] <- FALSE
  nodes
}

## Test that all nodes fall within the bounding box
test.bounding.box <- function(nodes, w, h, min.sep) {
  is.alive <- as.logical(nodes["alive",])
  x <- nodes["x",is.alive]
  y <- nodes["y",is.alive]
  r <- nodes["r",is.alive] + min.sep
  ok <- x - r >= 0 & x + r <= w & y - r >= 0 & y + r <= h
  kill <- which(is.alive)[!ok]
  nodes["alive", kill] <- FALSE
  nodes
}

## Test that nodes do not intersect with a set of rectangles.  The
## rectangles (rects) are a matrix with 4 rows (x, y, w, h) and as
## many columns as there are rectangles.
test.rects <- function(nodes, rects, min.sep) {
  test.rect <- function(rect, nodes) {
    circle.x <- nodes["x",]
    circle.y <- nodes["y",]
    circle.r <- nodes["r",] + min.sep
    closest.x <- clamp(circle.x, rect[1], rect[3])
    closest.y <- clamp(circle.y, rect[2], rect[4])
    square(circle.x - closest.x) + square(circle.y - closest.y) <
      square(circle.r)
  }

  is.alive <- as.logical(nodes["alive",])
  if (!any(is.alive))
    return(nodes)

  alive <- nodes[,is.alive,drop=FALSE]
  m <- apply(rects, 2, test.rect, alive)
  if (!is.matrix(m))
    m <- matrix(m, 1)
  if (nrow(m) != sum(is.alive))
    stop("This is getting annoying")
  kill <- which(is.alive)[rowSums(m) > 0]
  nodes["alive", kill] <- FALSE
  nodes
}

## Check that the most recently added node satisfies all criteria,
## otherwise turf it out.
check.new.node <- function(nodes) {
  n <- ncol(nodes)
  if (!nodes["alive",n])
    nodes <- nodes[,-n,drop=FALSE]
  nodes
}

grow.nodes <- function(nodes, max.size) {
  i <- as.logical(nodes["alive",])
  nodes["r", i] <- nodes["r", i] + nodes["growth", i]
  nodes["alive", which(i)[nodes["r", i] > max.size]] <- FALSE
  nodes
}

## Utility function for test.rect; restricts x to range
## [mid - w/2, mid + w/2]
clamp <- function(x, mid, w) {
  w2 <- w / 2
  min <- mid - w2
  max <- mid + w2
  x[x < min] <- min
  x[x > max] <- max
  x
}

## Utility function
square <- function(x)
  x * x
