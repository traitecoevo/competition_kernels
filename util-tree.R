model_jacobian_density <- function(lma, p, schedule) {
  target <- function(N) {
    p <- p$copy()
    p$seed_rain <- N
    landscape(lma, p, schedule)
  }
  jacobian(target, p$seed_rain, method="simple")
}
