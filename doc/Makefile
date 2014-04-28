RSCRIPT_PKGS := $(shell Rscript -e 'library(methods);writeLines(Sys.getenv("R_DEFAULT_PACKAGES"))')
RSCRIPT = Rscript --default-packages="${RSCRIPT_PKGS},methods"

all: simple_models.html

%.Rmd: %.R
	${RSCRIPT} -e "library(sowsear); sowsear('$<', 'Rmd')"
%.md: %.Rmd
	${RSCRIPT} -e "library(knitr); knit('$<')"
%.pdf: %.md
	pandoc $< -o $@
%.docx: %.md
	pandoc $< -o $@
%.html: %.md
	${RSCRIPT} -e "library(markdown);\
	  opts <- setdiff(markdownHTMLOptions(TRUE), 'base64_images');\
	markdownToHTML('$<', '$@', options=opts)"

clean:
	rm -f simple_models.{Rmd,html}
	rm -rf figure

WIKI_DIR = ../wiki
update_wiki: simple_models.md
	./.update_wiki.sh
publish_wiki: update_wiki
	./.publish_wiki.sh
