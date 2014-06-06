## Functions that we need:

default_parameters <- function() {
  p <- new(Parameters)
  p$set_parameters(list(patch_area=1.0))   # See issue #13
  p$set_control_parameters(fast.control()) # A bit faster
  p
}

default_schedule <- function(p) {
  t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)
  times0 <- cohort.introduction.times(t.max)
  schedule.from.times(times0)
}

default_strategy <- function() {
  new(Strategy)
}

get_or_default <- function(obj, f, ...) {
  if (is.null(obj)) {
    f(...)
  } else {
    obj$copy()
  }
}

max_growth_rate <- function(lma, p=NULL, s=NULL, schedule=NULL,
                            seed_rain=1) {
  p <- get_or_default(p, default_parameters)
  s <- get_or_default(s, default_strategy)
  schedule <- get_or_default(schedule, default_schedule, p)

  if (length(schedule$ode_times) == 0) {
    lma.mid <- exp(mean(log(range(lma))))
    schedule$ode_times <- add_ode_times(lma.mid, p, schedule, s)
  }
  
  landscape.empty(lma, p, schedule)  
}

carrying_capacity <- function(lma, ..., parallel=FALSE) {
  if (parallel) {
    unlist(parallel::mclapply(lma, carrying_capacity1, ...,
                              mc.preschedule=FALSE))
  } else {
    unlist(lapply(lma, carrying_capacity1, ...))
  }
}

carrying_capacity1 <- function(lma, p=NULL, s=NULL, schedule=NULL,
                               seed_rain=1,
                               eps=1e-3, build.args=list()) {
  if (length(lma) != 1) {
    stop("Expected single lma value")
  }
  p <- get_or_default(p, default_parameters)
  s <- get_or_default(s, default_strategy)
  schedule <- get_or_default(schedule, default_schedule, p)

  p <- add_strategy_lma(lma, p, s, seed_rain)

  ## We might want to add the times here?
  res <- equilibrium.seed.rain(p, schedule, 10,
                               build.args=build.args, eps=eps)
  unname(mean(res$seed_rain))
}

make_invasion_fitness <- function(lma, seed_rain,
                                  p=NULL, s=NULL, schedule=NULL,
                                  nsteps=20, eps=1e-3, verbose=TRUE) {
  p <- get_or_default(p, default_parameters)
  s <- get_or_default(s, default_strategy)
  schedule <- get_or_default(schedule, default_schedule, p)
  p <- add_strategy_lma(lma, p, s, seed_rain)
  schedule <- build.schedule(p, schedule, nsteps, eps, verbose=TRUE)
  schedule$ode_times <- attr(schedule, "ebt")$ode_times
  function(lma.mutant) {
    log(landscape(lma.mutant, p, schedule))
  }
}

add_ode_times <- function(lma, p, schedule, s, seed_rain=1) {
  p <- add_strategy_lma(lma, p, s, seed_rain)
  ebt <- new(EBT, p)
  ebt$cohort_schedule <- schedule$copy()
  ebt$run()
  ebt$ode_times
}

add_strategy_lma <- function(lma, p, s, seed_rain=NULL) {
  if (p$size > 0) {
    stop("Must have empty environment")
  }
  p <- p$copy()
  s <- s$copy()
  s$set_parameters(list(lma=lma))
  p$add_strategy(s)
  if (!is.null(seed_rain)) {
    p$seed_rain <- seed_rain
  }
  p
}

## Create a simple spline of x and y in log-log space, but that we can
## use in non-log space.
loglog_splinefun <- function(x, y, log="xy") {
  if (!identical(log, "")) { # wow
    log <- match.arg(log, c("x", "y", "xy"))
  }
  log.x <- grepl("x", log, fixed=TRUE)
  log.y <- grepl("y", log, fixed=TRUE)
  tr.in  <- if (log.x) base::log else identity
  tr.out <- if (log.y) exp else identity
  f <- splinefun(if (log.x) log(x) else x,
                 if (log.y) log(y) else y)
  function(x) {
    tr.out(f(tr.in(x)))
  }
}

## Copied over from direct.R, but only the type that seems more sensible
compute_alpha <- function(w.i, r.i, K.i, N.r) {
  (1 - 1 / r.i * w.i) * K.i / N.r
}
