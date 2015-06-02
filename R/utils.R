## TODO: probably better in absolute size?
black_bar <- function(x, y, col="black") {
  usr <- par("usr")
  rect(x[1], usr[3], x[2], usr[3] + y, col=col, border=col)
}

expand_range <- function(r, p) {
  if (length(p) == 1L) {
    p <- rep_len(p, 2L)
  }
  d <- diff(r)
  r + d * c(-1, 1) * p
}
