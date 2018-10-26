main=genesetdp

all: ${main}.pdf coverletter.pdf response.pdf diff.pdf  # supplement.pdf

${main}.pdf: ${main}.tex # ${main}.bib # img/*

diff.tex: ${main}.tex
	git show plosone_2:./$? > /tmp/$?
	latexdiff /tmp/$? $? > diff.tex


%.pdf: %.tex
	pdflatex -shell-escape $*.tex
	if ( grep -q citation $*.aux ) ; then \
		bibtex $* ; \
		pdflatex -shell-escape $*.tex ; \
	fi
	pdflatex -shell-escape $*.tex

html:
	rm -fr ${main}/
	latex2html -split 0 ${main}.tex

zip:
	zip -r ${main}.zip *.tex *.bib img makefile

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *~ .DS_Store
	rm -fr ${main}/
