## Little dance to do here:
all: ${TARGETS}
TARGETS = lab_meeting.pdf
include $(shell Rscript -e "ghwiki:::rules.mk()")

lab_meeting.pdf: lab_meeting.tex ${TARGETS_MD}
	xelatex $<

history.pdf: history.md the-american-naturalist.csl
	pandoc -s --csl the-american-naturalist.csl --filter pandoc-citeproc --bibliography=refs.bib -o $@ $<

the-american-naturalist.csl:
	curl -k https://raw.githubusercontent.com/citation-style-language/styles/master/the-american-naturalist.csl > $@
