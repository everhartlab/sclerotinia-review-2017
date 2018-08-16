TARGETS  := $(wildcard *.md)
HTML     := $(patsubst %.md,html/%.html,$(TARGETS))
PDF      := $(patsubst %.md,pdf/%.pdf,$(TARGETS))
TEX      := $(patsubst %.md,pdf/%.tex,$(TARGETS))
DOCX     := $(patsubst %.md,docx/%.docx,$(TARGETS))
TEMPLATE := Helvetica-Garamond.docx
ifeq ("$(shell pandoc --version | grep -c 'pandoc 2')","1")
	PANDARGS := --from markdown+smart --standalone -F pandoc-citeproc
else
	PANDARGS := --from markdown --smart --standalone -F pandoc-citeproc
endif

.PHONY: all
all : review.docx review.pdf

pdf html docx :
	mkdir -p $@

pdf/%.tex : %.md pandoc.css references.bib | pdf
	pandoc --to latex \
      $(PANDARGS) \
      --css pandoc.css \
      --output $@ $<

review.tex : pdf/review.tex correct-authors.tex | pdf
	sed "s/^.author/\\\\input{correct-authors.tex}%\\\\author/" $< \
	| sed "s/^.institute/%\\\\institute/" > review.tex

review.docx : review.tex
	pandoc --to docx \
       --reference-doc=$(TEMPLATE) \
			 --output $@ $<

review.pdf : review.tex
	sed "s/.png/.pdf/" $< | sed "s/-small.pdf/.png/" > review-pdf.tex
	pdflatex review-pdf
	pdflatex review-pdf
	mv review-pdf.pdf review.pdf

.PHONY: clean
clean:
	$(RM) -r html \
         pdf \
         docx \
         review-pdf.{tex,out,aux,log} \
         review.{tex,out,aux,log} \
         *.aux

.PHONY: cleanall
cleanall : clean
	$(RM) -r review.pdf review.docx
