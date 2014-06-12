seq_log <- function (from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}
seq_sqrt <- function(from, to, length.out) {
  (seq(sqrt(from), sqrt(to), length.out = length.out))^2
}

to_dev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if (verbose)
    cat(sprintf("Creating %s\n", filename))
  dev(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
  invisible()
}

to_pdf <- function(expr, filename, ..., cairo=FALSE, pointsize=12) {
  dev <- if (cairo) CairoPDF else pdf
  to_dev(expr, dev, filename, ..., pointsize=pointsize)
}

to_cairo_pdf <- function(expr, filename, ..., pointsize=12,
                         family="Hoefler Text", bg="transparent") {
  to_dev(expr, cairo_pdf, filename, ..., pointsize=pointsize,
         family=family, bg=bg)
}

set_cairo_fonts <- function(name, bold="Bold", italic="Italic",
                            bold.italic="Bold Italic") {
  fmt <- "%s:style=%s"
  CairoFonts(regular    = name,
             bold       = sprintf(fmt, name, bold),
             italic     = sprintf(fmt, name, italic),
             bolditalic = sprintf(fmt, name, bold.italic))
}

## These should be tweaked to match the font pairings in format.tex
set_cairo_fonts_hoefler <- function() {
  set_cairo_fonts("Hoefler Text")
}
set_cairo_fonts_kreon <- function() {
  set_cairo_fonts("Kreon", italic="Light")
}
set_cairo_fonts_cabin <- function() {
  set_cairo_fonts("Cabin")
}
set_cairo_fonts_lato <- function() {
  set_cairo_fonts("Lato")
}

last <- function(x) {
  x[[length(x)]]
}

expand_range <- function(r, p=0.03, down=TRUE, up=TRUE) {
  r + c(if (down) -1 else 0, if (up) 1 else 0) * diff(r) * p
}
