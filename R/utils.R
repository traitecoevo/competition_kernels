## TODO: probably better in absolute size?
black_bar <- function(x, y, col="black") {
  usr <- par("usr")
  rect(x[1], usr[3], x[2], usr[3] + y, col=col, border=col)
}

expand_range <- function(r, p) {
  if (length(p) == 1L) {
    p <- rep_len(p, 2L)
  }
  d <- diff(r)
  r + d * c(-1, 1) * p
}

build_pdf_knit <- function(MD_file){

    require(knitr)
    name_file <- sub(".Rmd", "", MD_file)
    knit(paste0(name_file, '.Rmd'))
    args <- paste0(name_file, ".md", " -o ", name_file,".pdf")
    print(args)
    suppressWarnings(system2("pandoc", args, stdout=TRUE, stderr=TRUE))
  }

