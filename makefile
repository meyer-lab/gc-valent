fdir = ./Manuscript/Figures
tdir = ./Manuscript/Templates
pan_common = -F pandoc-crossref -F pandoc-citeproc --filter=$(tdir)/figure-filter.py -f markdown ./Manuscript/Text/*.md

.PHONY: clean test all testprofile testcover doc testcpp

all: ckine/ckine.so Manuscript/index.html Manuscript/Manuscript.pdf Manuscript/Manuscript.docx Manuscript/CoverLetter.docx

$(fdir)/figure%.svg: genFigures.py graph_all.svg
	mkdir -p ./Manuscript/Figures
	python3 genFigures.py $*

$(fdir)/figure%pdf: $(fdir)/figure%svg
	rsvg-convert -f pdf $< -o $@

$(fdir)/figure%eps: $(fdir)/figure%svg
	rsvg-convert -f eps $< -o $@

graph_all.svg: ckine/data/graph_all.gv
	dot $< -Tsvg -o $@

Manuscript/Manuscript.pdf: Manuscript/Manuscript.tex $(fdir)/figure1.pdf $(fdir)/figure2.pdf $(fdir)/figure3.pdf $(fdir)/figure4.pdf $(fdir)/figureS1.pdf $(fdir)/figureS2.pdf $(fdir)/figureS3.pdf
	(cd ./Manuscript && latexmk -xelatex -f -quiet)
	rm -f ./Manuscript/Manuscript.b* ./Manuscript/Manuscript.aux ./Manuscript/Manuscript.fls

Manuscript/index.html: Manuscript/Text/*.md
	pandoc -s $(pan_common) -t html5 --mathjax -c ./Templates/kultiad.css --template=$(tdir)/html.template -o $@

Manuscript/Manuscript.docx: Manuscript/Text/*.md
	mkdir -p ./Manuscript/Figures
	cp -R $(fdir) ./
	pandoc -s $(pan_common) -o $@
	rm -r ./Figures

Manuscript/Manuscript.tex: Manuscript/Text/*.md Manuscript/index.html
	pandoc -s $(pan_common) --template=$(tdir)/default.latex --pdf-engine=xelatex -o $@

Manuscript/CoverLetter.docx: Manuscript/CoverLetter.md
	pandoc -f markdown $< -o $@

Manuscript/CoverLetter.pdf: Manuscript/CoverLetter.md
	pandoc --pdf-engine=xelatex --template=/Users/asm/.pandoc/letter-templ.tex $< -o $@

clean:
	rm -f ./Manuscript/Manuscript.* ./Manuscript/index.html Manuscript/CoverLetter.docx Manuscript/CoverLetter.pdf ckine/libckine.debug.so
	rm -f $(fdir)/Figure* profile.p* stats.dat .coverage nosetests.xml coverage.xml ckine.out ckine/cppcheck testResults.xml
	rm -rf html ckine/*.dSYM doxy.log graph_all.svg

test:
	nosetests3 -s --with-timer --timer-top-n 20

testcover:
	nosetests3 --with-xunit --with-xcoverage --cover-inclusive --cover-package=ckine -s --with-timer --timer-top-n 5

stats.dat:
	nosetests3 -s --with-cprofile

testprofile: stats.dat
	pyprof2calltree -i stats.dat -k

doc:
	doxygen Doxyfile
