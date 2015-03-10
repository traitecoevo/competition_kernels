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

build_pdf <- function(texfile) {
  owd <- setwd("ms")
  on.exit(setwd(owd))
  texfile <- basename(texfile)

  latex <- Sys.which("pdflatex")
  bibtex <- Sys.which("bibtex")

  if (latex == "" || bibtex == "") {
    stop("latex / bibtex not found in path")
  }

  ## These can be generalised a bunch:
  run_latex <- function() {
    args <- c("-interaction=nonstopmode", texfile)
    res <- suppressWarnings(system2(latex, args, stdout=TRUE, stderr=TRUE))
    code <- attr(res, "status")
    if (!is.null(code) && code != 0L) {
      cat(res, sep="\n")
      stop("pdflatex failed with code ", code, ": message above")
    }
  }
  run_bibtex <- function() {
    args <- sub(".tex", "", texfile)
    res <- suppressWarnings(system2(bibtex, args, stdout=TRUE, stderr=TRUE))
    code <- attr(res, "status")
    if (!is.null(code) && code != 0L) {
      cat(res, sep="\n")
      stop("bibtex failed with code ", code, ": message above")
    }
  }

  run_latex()
  run_bibtex()
  run_latex()
  run_latex()

  aux_files <- paste0(tools::file_path_sans_ext(texfile),
                     c(".log", ".aux", ".bbl", ".blg", ".fls", ".out",
                       ".fdb_latexmk"))
  file.remove(aux_files[file.exists(aux_files)])
  invisible(NULL)
}
