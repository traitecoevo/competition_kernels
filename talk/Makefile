all: slides.pdf

slides.pdf: slides.tex slides_format.tex
	xelatex $<

tidy:
	rm -f slides.{aux,log,nav,out,snm,toc}

clean: tidy
	rm -f slides.pdf

.PHONY: all tidy clean
