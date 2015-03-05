## Functions for computing various bits of the models...

competition_LV <- function(w_i, r_i, K_i, N_r) {
  (1 - 1 / r_i * w_i) * K_i / N_r
}
