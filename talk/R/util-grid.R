grid_xaxis_simple <- function(x, y, just="top", ...) {
  grid.text(format(x), unit(x, "native"), y, just=just, ...)
}

grid_yaxis_simple <- function(x, y, just="right", ...) {
  grid.text(format(y), x, unit(y, "native"), just=just, ...)
}

grid_xlab_simple <- function(lab, y, just="top", ...) {
  grid.text(lab, unit(.5, "npc"), y, just=just, ...)
}
grid_ylab_simple <- function(lab, x, just="bottom", ...) {
  grid.text(lab, x, unit(.5, "npc"), just=just, rot=90, ...)
}
