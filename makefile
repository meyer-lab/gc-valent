SHELL := /bin/bash
fdir = ./Manuscript/Figures
tdir = ./common/templates
pan_common = -F pandoc-crossref -F pandoc-citeproc --filter=$(tdir)/figure-filter.py -f markdown ./Manuscript/Text/*.md

flist = B1 B2 B3 B4 B5

.PHONY: clean test all testprofile testcover testcpp autopep spell leaks profilecpp

all: Manuscript/Manuscript.pdf Manuscript/Manuscript.docx Manuscript/CoverLetter.docx pylint.log

venv: venv/bin/activate

venv/bin/activate: requirements.txt
	test -d venv || virtualenv --system-site-packages venv
	. venv/bin/activate && pip install -Ur requirements.txt
	. venv/bin/activate && python3 -c "import julia; julia.install()"
	. venv/bin/activate && julia -e "using Pkg; Pkg.update(); Pkg.add(PackageSpec(url=\"https://github.com/meyer-lab/gcSolver.jl\"))"
	touch venv/bin/activate

$(fdir)/figure%.svg: venv genFigures.py ckine/ckine.so ckine/figures/figure%.py
	mkdir -p ./Manuscript/Figures
	. venv/bin/activate && THEANO_FLAGS='mode=FAST_COMPILE' ./genFigures.py $*

$(fdir)/figure%pdf: $(fdir)/figure%svg
	rsvg-convert --keep-image-data -f pdf $< -o $@

$(fdir)/figure%eps: $(fdir)/figure%svg
	rsvg-convert --keep-image-data -f eps $< -o $@

Manuscript/Manuscript.pdf: Manuscript/Text/*.md $(patsubst %, $(fdir)/figure%.pdf, $(flist))
	pandoc -s $(pan_common) --template=$(tdir)/default.latex --pdf-engine=xelatex -o $@

Manuscript/Manuscript.docx: Manuscript/Text/*.md $(patsubst %, $(fdir)/figure%.eps, $(flist))
	cp -R $(fdir) ./
	pandoc -s $(pan_common) -o $@
	rm -r ./Figures

Manuscript/CoverLetter.docx: Manuscript/CoverLetter.md
	pandoc -f markdown $< -o $@

Manuscript/CoverLetter.pdf: Manuscript/CoverLetter.md
	pandoc --pdf-engine=xelatex --template=/Users/asm/.pandoc/letter-templ.tex $< -o $@

autopep:
	autopep8 -i -a --max-line-length 200 ckine/*.py ckine/figures/*.py

clean:
	rm -f ./Manuscript/Manuscript.* Manuscript/CoverLetter.docx Manuscript/CoverLetter.pdf
	rm -f $(fdir)/figure* profile.p* stats.dat .coverage nosetests.xml coverage.xml testResults.xml
	rm -rf html ckine/*.dSYM doxy.log graph_all.svg venv
	find -iname "*.pyc" -delete

spell: Manuscript/Text/*.md
	pandoc --lua-filter common/templates/spell.lua Manuscript/Text/*.md | sort | uniq -ic

test: venv
	. venv/bin/activate && python -m pytest

testcover: venv
	. venv/bin/activate && THEANO_FLAGS='mode=FAST_COMPILE' python -m pytest --junitxml=junit.xml --cov-branch --cov=ckine --cov-report xml:coverage.xml

pylint.log: venv common/pylintrc
	. venv/bin/activate && (pylint --rcfile=./common/pylintrc ckine > pylint.log || echo "pylint3 exited with $?")
