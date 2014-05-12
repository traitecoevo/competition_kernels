## Little dance to do here:
all: ${TARGETS}
TARGETS = lab_meeting.pdf
include $(shell Rscript -e "ghwiki:::rules.mk()")

lab_meeting.pdf: lab_meeting.tex ${TARGETS_MD}
	xelatex $<

