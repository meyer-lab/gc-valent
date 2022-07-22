SHELL := /bin/bash

flist = C1 C2 C3 C4 C5 S2 S3 S4 S5 S6

.PHONY: clean test all

all: $(patsubst %, output/figure%.svg, $(flist))

output/figure%.svg: ckine/figures/figure%.py
	mkdir -p ./output
	poetry run fbuild $*

clean:
	rm -rf output

test:
	poetry run pytest -s -v -x
