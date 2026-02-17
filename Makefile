PAPER = main
TEXFILES = $(PAPER).tex
BIBFILE = references.bib

.PHONY: all clean pdf

all: pdf

pdf: $(PAPER).pdf

$(PAPER).pdf: $(TEXFILES) $(BIBFILE)
	pdflatex $(PAPER)
	bibtex $(PAPER)
	pdflatex $(PAPER)
	pdflatex $(PAPER)

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.synctex.gz *.fdb_latexmk *.fls

distclean: clean
	rm -f $(PAPER).pdf

view: $(PAPER).pdf
	start $(PAPER).pdf

# WSL compilation target
wsl:
	wsl -u agent -- bash -lc "cd /mnt/d/openclaw/biodyn-nmi-paper && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex"