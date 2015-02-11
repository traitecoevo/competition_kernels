TARGETS = competition-kernels.pdf

all: ${TARGETS}

# Manual:

# TARGETS_BASE  := $(shell basename -a -s .pdf ${TARGETS})
# TARGETS_BASE2 := $(shell echo ${TARGETS_BASE} | tr "\\n" "," | sed 's/,$$//')
#
# clean:
# 	rm -f ${TARGETS}
# tidy:
# 	rm -f {${TARGETS_BASE2}}.{aux,bbl,blg,log}
# %.pdf: %.tex
# 	pdflatex $<
# 	bibtex $(basename -s.tex $<)
# 	pdflatex $<
# 	pdflatex $<

# Trying latexmk -- see here: http://tex.stackexchange.com/a/40759
%.pdf: %.tex
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make $<

clean:
	latexmk -C

tidy:
	latexmk -c

.PHONY: all tidy clean
