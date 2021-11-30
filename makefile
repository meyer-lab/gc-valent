SHELL := /bin/bash

flist = C2 C3 C4 S1 S2 S3 S4

.PHONY: clean test all

all: $(patsubst %, output/figure%.svg, $(flist))

output/figure%.svg: ckine/figures/figure%.py
	mkdir -p ./output
	poetry run fbuild $*

clean:
	rm -rf output

test: venv
	poetry run pytest -s -v -x
