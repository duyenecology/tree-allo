GIT = 047f6e7
GIT2 = 3ac8a60
COVER = ms/cover_letter
MAIN = ms/ms
SI = sup/sup

# all: $(MAIN).pdf   $(MAIN)_diff.pdf
all: $(MAIN).pdf $(MAIN).docx $(SI).pdf $(SI).docx  $(MAIN)_diff.pdf $(SI)_diff.pdf
diff: $(MAIN)_diff.pdf
pdf: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).qmd
	quarto render $< --to pdf

$(MAIN).docx: $(MAIN).qmd ms/my_template.docx
	quarto render $< --to docx

$(SI).pdf: $(SI).qmd
	quarto render $< --to pdf

$(SI).docx: $(SI).qmd
	quarto render $< --to docx

$(COVER).pdf: $(COVER).qmd
	quarto render $< --to pdf

# $(COVER).docx: $(COVER).qmd
# 	quarto render $< --to docx

# because different quarto verions produce different latex files
$(MAIN)_diff.pdf: $(MAIN).tex
	git show $(GIT):$(MAIN).qmd > $(MAIN)_old.qmd
	sed -i 's/link-citations: yes/link-citations: true/' $(MAIN)_old.qmd
	quarto render ms/ms_old.qmd --to pdf && \
	cd ms && \
	latexdiff ms_old.tex ms.tex > ms_diff.tex && \
	xelatex ms_diff.tex
	rm $(MAIN)_old.*

$(SI)_diff.pdf: $(SI).tex
	git show $(GIT2):$(SI).qmd > $(SI)_old.qmd
	sed -i 's/link-citations: yes/link-citations: true/' $(SI)_old.qmd
	quarto render sup/sup_old.qmd --to pdf && \
	cd sup && \
	latexdiff sup_old.tex sup.tex > sup_diff.tex && \
	xelatex sup_diff.tex
	rm $(SI)_old.*

.PHONY: clean
clean:
	rm -f ms/*.tuc \
	ms/*.log \
	rm -rf ms/cache/*

