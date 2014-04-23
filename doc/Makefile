RSCRIPT_PKGS := $(shell Rscript -e 'library(methods);writeLines(Sys.getenv("R_DEFAULT_PACKAGES"))')
RSCRIPT = Rscript --default-packages="${RSCRIPT_PKGS},methods"

all: simple_models.html

%.Rmd: %.R
	${RSCRIPT} -e "library(sowsear); sowsear('$<', 'Rmd')"
%.md: %.Rmd ${DATA_PROCESSED}
	${RSCRIPT} -e "library(knitr); knit('$<')"
%.html: %.md
	${RSCRIPT} -e "library(markdown);\
	  opts <- setdiff(markdownHTMLOptions(TRUE), 'base64_images');\
	markdownToHTML('$<', '$@', options=opts)"

clean:
	rm -f simple_models.{Rmd,html}
	rm -rf figure
