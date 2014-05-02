## All the models are in the Revolve package:
library(Revolve)

## And we'll need tools in numDeriv
library(numDeriv)

## This function gets used all over, and is explained in the Dieckmann
## and Dobeli section.
model_jacobian_density <- function(x, sys, m, ...) {
  jacobian(function(y) m$fitness(x, sys$x, y), sys$y, ...)
}
