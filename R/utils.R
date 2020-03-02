black_bar <- function(x, y, col="black", yp=0.025) {
  usr <- par("usr")
  y <- (usr[4] - usr[3]) * yp
  rect(x[1], usr[3], x[2], usr[3] + y, col=col, border=col)
}

expand_range <- function(r, p) {
  if (length(p) == 1L) {
    p <- rep_len(p, 2L)
  }
  d <- diff(r)
  r + d * c(-1, 1) * p
}

label <- function(text, px=0.03, py=NULL, ..., adj=c(0, 1)) {
  if (is.null(py)) {
    fin <- par("fin")
    r <- fin[[1]] / fin[[2]]
    if (r > 1) { # x is longer.
      py <- 1 - px
      px <- (1 - py) / r
    } else {
      py <- 1 - px * r
    }
  }
  usr <- par("usr")
  x <- usr[1] + px*(usr[2] - usr[1])
  y <- usr[3] + py*(usr[4] - usr[3])

  ## NOTE: base 10 log:
  if (par("xlog")) {
    x <- 10^x
  }
  if (par("ylog")) {
    y <- 10^y
  }

  text(x, y, text, adj=adj, ...)
}

label_panel <- function(i, ...) {
  label(sprintf("%s)", letters[i]), ...)
}

render_pdf <- function(filename) {
  owd <- setwd(dirname(filename))
  on.exit(setwd(owd))
  filename <- basename(filename)
  dest <- paste0(sub("\\.md", "", filename), ".pdf")
  call_system(Sys_which("pandoc"), c(filename, "-o", dest))
}
