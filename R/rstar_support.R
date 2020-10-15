## Some basic error checking...
assert_scalar <- function(x, name=deparse(substitute(x))) {
  if (length(x) != 1) {
    stop(sprintf("%s must be a scalar", name), call.=FALSE)
  }
}
colMins <- function(x) {
  apply(x, 2, min)
}

equilibrium <- function(dydt, y, pars, method="nleqslv",
                        init_time=200, max_time=1e5) {
  method <- match.arg(method, c("runsteady", "nleqslv", "simple"))

  if (init_time > 0) {
    y <- unname(.lsoda(y, c(0, init_time), dydt, pars)[2, -1, drop=TRUE])
  }

  if (method == "runsteady") {
    ans <- rootSolve::runsteady(y, derivs_deSolve(dydt), pars,
                                times=c(0, Inf), hmax=1)$y
  } else if (method == "nleqslv") {
    ans <- nleqslv::nleqslv(y, function(y) dydt(0, y, pars), global="none")$x
  } else if (method == "simple") {
    tt <- c(init_time, max_time)
    ans <- unname(.lsoda(y, tt, dydt, pars)[2, -1, drop=FALSE])
  }
  ans
}

.lsoda <- function(y, times, func, ...) {
  deSolve::lsoda(y, times, derivs_deSolve(func), ...)
}

derivs_deSolve <- function(func) {
  function(...) list(func(...))
}

quadratic_roots <- function(a, b, c) {
  (-b + c(-1, 1) * sqrt(b*b - 4*a*c))/(2 * a)
}

make_transparent <- function(col, opacity=.5) {
  alpha <- opacity
  if ( length(alpha) > 1 && any(is.na(alpha)) ) {
    n <- max(length(col), length(alpha))
    alpha <- rep(alpha, length.out=n)
    col <- rep(col, length.out=n)
    ok <- !is.na(alpha)
    ret <- rep(NA, length(col))
    ret[ok] <- make_transparent(col[ok], alpha[ok])
    ret
  } else {
    tmp <- col2rgb(col)/255
    rgb(tmp[1,], tmp[2,], tmp[3,], alpha=alpha)
  }
}

mix <- function(col1, col2, p=0.5) {
  m <- col2rgb(col1, alpha=TRUE)
  m2 <- col2rgb(rep(col2, length.out=length(col1)), alpha=TRUE)
  m3 <- (m * p + m2 * (1-p))/255
  alpha <- if (all(m3[4,] == 1)) NULL else m3[4,]
  rgb(m3[1,], m3[2,], m3[3,], alpha)
}

abcline <- function(x, y, m, ...) {
  abline(y - x * m, m, ...)
}
