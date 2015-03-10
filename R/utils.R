## TODO: probably better in absolute size?
black_bar <- function(x, y, col="black") {
  usr <- par("usr")
  rect(x[1], usr[3], x[2], usr[3] + y, col=col, border=col)
}

build_pdf <- function(texfile) {
  owd <- setwd("ms")
  on.exit(setwd(owd))
  texfile <- basename(texfile)
  command <- sprintf('latexmk -pdf -pdflatex=\"%s\" %s',
  									"pdflatex -interaction=nonstopmode", texfile)
 	res <- system(command, intern=FALSE)
 	# TODO: how to catch errors? Above does not quite work because res behaves
  # stranegly (no contents, not length)
 	#	res <- system(command, intern=TRUE)
  # code <- attr(res, "status")
  # browser()
  # if (code != 0L) {
  # 	stop("Error running latexmk: ", as.character(attr(res, "errmsg", TRUE)))
  # }
  aux_files <- paste0(tools::file_path_sans_ext(texfile),
                     c(".log", ".aux", ".bbl", ".blg", ".fls", ".out",
                       ".fdb_latexmk"))
  file.remove(aux_files[file.exists(aux_files)])
}
