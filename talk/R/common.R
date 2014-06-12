# https://creativemarket.com/Marish/4050-Social-network-infographics-set/screenshots/#screenshot1

# #574134 -- brown (light)
# #473629 -- brown (dark)
# #d54b1a -- orange (light)
# #bf4317 -- orange (dark)
# #e3a72f -- yellow (light)
# #c6891d -- yellow (dark)
# #058789 -- blue (light)
# #04797a -- blue (dark)
# #f0ecc9 -- cream (light)
# #c7c094 -- creak (dark)

colours <- function() {
  ## devtools::install_github("karthik/wesanderson")
  ## wesanderson::wes.palette(4, "Royal1")
  cols.wa <- c("#899DA4", "#C93312", "#FAEFD1", "#DC863B")
  cols <- list(bg=cols.wa[[3]], fg=cols.wa[[1]],
               hl=cols.wa[[2]], con=cols.wa[[4]])
}
