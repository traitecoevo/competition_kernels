all: slides.pdf

figures:
	mkdir figs
	Rscript figures.R

slides.pdf: slides.tex slides_format.tex
	xelatex $<
	xelatex $<

tidy:
	rm -f slides.{aux,log,nav,out,snm,toc}

clean: tidy
	rm -f slides.pdf
	rm -f pics/leaf.{eps,xml} pics/tree.{eps,xml} pics/trees.{eps,xml}
	rm -f figs/*

.PHONY: all tidy clean
