rstar_plot <- function(x, parameters,
                       xlim=c(0, 1), ylim=c(0, 1),
                       col=seq_len(ncol(x))) {
  if (parameters$k != 2) {
    stop("Can only plot 2 resource case at the moment")
  }
  n <- ncol(x)
  col <- rep(col, length.out=n)

  rs <- rstar_Rstar(x, parameters)
  plot(NA, xlim=xlim, ylim=ylim, xlab="R1", ylab="R2")
  abline(v=rs[1,], lty=3, col=col)
  abline(h=rs[2,], lty=3, col=col)
  segments(rs[1,], rs[2,], rs[1,], par("usr")[4], col=col)
  segments(rs[1,], rs[2,], par("usr")[2], rs[2,], col=col)

  R.eq <- apply(rs, 1, max)
  # This only makes sense for some cases?
  points(R.eq[1], R.eq[2], pch=19)

  C <- parameters$C(x)
  for (i in seq_len(n)) {
    abcline(R.eq[1], R.eq[2], C[2,i] / C[1,i], col=col[i], lty=2)
  }
}

rstar_trajectory <- function(x, N, parameters,
                             col=seq_len(ncol(x)), eps=1e-6,
                             t_max=100, t_len=101, S=NULL,
                             col_died="grey") {
  if (ncol(x) > 2) {
    stop("Only working for up to two species so far")
  }
  if (!is.null(S)) {
    if (length(S) != length(parameters$S)) {
      stop("Invalid length S")
    }
    parameters$S <- S
  }

  t <- seq(0, t_max, length.out=t_len)
  eq <- rstar_equilibrium(x, N, parameters)
  tr <- rstar_run(t, x, N, parameters)
  survived <- eq$N > eps

  if (sum(survived) == 1) {
    col <- col[which(survived)]
  } else if (sum(survived) == 2) {
    col <- mix(col[1], col[2])
  } else {
    col <- col_died
  }

  lines(tr$R[,1], tr$R[,2], lty=2, col=col)
  points(c(parameters$S[1], eq$R[1]),
         c(parameters$S[2], eq$R[2]), pch=19, col=col, cex=.5)
}
