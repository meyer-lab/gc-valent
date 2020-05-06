SHELL := /bin/bash

flist = J1 J2

.PHONY: clean test all testprofile testcover spell

all: output/manuscript.html pylint.log

venv: venv/bin/activate

venv/bin/activate: requirements.txt
	test -d venv || virtualenv --system-site-packages venv
	. venv/bin/activate && pip install -Uqr requirements.txt
	[credential "https://github.com/meyer-lab/gcSolver.jl.git"]
	julia -e 'using Pkg; Pkg.add("PyCall"); Pkg.add(PackageSpec(url="https://borcuttjahns:987faeec2cac1aa066960e9cdc6345a58d3bacb1@github.com/meyer-lab/gcSolver.jl.git")); Pkg.build(); Pkg.precompile()'
	. venv/bin/activate && python3 -c "import julia; julia.install()"
	touch venv/bin/activate

output/figure%.svg: venv genFigures.py ckine/figures/figure%.py
	mkdir -p ./output
	. venv/bin/activate && ./genFigures.py $*

output/manuscript.md: venv manuscript/*.md
	. venv/bin/activate && manubot process --content-directory=./manuscript/ --output-directory=./output --log-level=INFO

output/manuscript.html: venv output/manuscript.md $(patsubst %, output/figure%.svg, $(flist))
	mkdir output/output
	cp output/*.svg output/output/
	. venv/bin/activate && pandoc --verbose \
		--from=markdown --to=html5 --filter=pandoc-fignos --filter=pandoc-eqnos --filter=pandoc-tablenos \
		--bibliography=output/references.json \
		--csl=common/templates/manubot/style.csl \
		--metadata link-citations=true \
		--include-after-body=common/templates/manubot/default.html \
		--include-after-body=common/templates/manubot/plugins/table-scroll.html \
		--include-after-body=common/templates/manubot/plugins/anchors.html \
		--include-after-body=common/templates/manubot/plugins/accordion.html \
		--include-after-body=common/templates/manubot/plugins/tooltips.html \
		--include-after-body=common/templates/manubot/plugins/jump-to-first.html \
		--include-after-body=common/templates/manubot/plugins/link-highlight.html \
		--include-after-body=common/templates/manubot/plugins/table-of-contents.html \
		--include-after-body=common/templates/manubot/plugins/lightbox.html \
		--mathjax \
		--variable math="" \
		--include-after-body=common/templates/manubot/plugins/math.html \
		--include-after-body=common/templates/manubot/plugins/hypothesis.html \
		--output=output/manuscript.html output/manuscript.md

clean:
	mv output/requests-cache.sqlite requests-cache.sqlite || true
	rm -rf prof output coverage.xml .coverage .coverage* junit.xml coverage.xml profile profile.svg pylint.log
	mkdir output
	mv requests-cache.sqlite output/requests-cache.sqlite || true
	rm -f profile.p* stats.dat .coverage nosetests.xml coverage.xml testResults.xml
	rm -rf html doxy.log graph_all.svg venv ./ckine/data/flow
	find -iname "*.pyc" -delete

spell: manuscript/*.md
	pandoc --lua-filter common/templates/spell.lua manuscript/*.md | sort | uniq -ic

download:
	mkdir -p ./ckine/data/flow
	wget -nv -P ./ckine/data/flow/ "https://syno.seas.ucla.edu:9001/gc-cytokines/2019-03-15 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - NK plate.zip"
	wget -nv -P ./ckine/data/flow/ "https://syno.seas.ucla.edu:9001/gc-cytokines/2019-04-18 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - Treg plate - NEW PBMC LOT.zip"
	wget -nv -P ./ckine/data/flow/ "https://syno.seas.ucla.edu:9001/gc-cytokines/2019-04-23 Receptor Quant - Beads.zip"
	wget -nv -P ./ckine/data/flow/ "https://syno.seas.ucla.edu:9001/gc-cytokines/4-23_4-26_Receptor quant.zip"
	unzip -qd ./ckine/data/flow/ './ckine/data/flow/2019-03-15 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - NK plate.zip'
	unzip -qd ./ckine/data/flow/ './ckine/data/flow/2019-04-18 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - Treg plate - NEW PBMC LOT.zip'
	unzip -qd ./ckine/data/flow/ './ckine/data/flow/2019-04-23 Receptor Quant - Beads.zip'
	unzip -qd ./ckine/data/flow/ './ckine/data/flow/4-23_4-26_Receptor quant.zip'

test: venv
	. venv/bin/activate && pytest

testcover: venv
	. venv/bin/activate && pytest --junitxml=junit.xml --cov-branch --cov=ckine --cov-report xml:coverage.xml

pylint.log: venv common/pylintrc
	. venv/bin/activate && (pylint --rcfile=./common/pylintrc ckine > pylint.log || echo "pylint3 exited with $?")
