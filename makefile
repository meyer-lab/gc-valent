SHELL := /bin/bash

flist = C2 C3 C4 S2 S3

.PHONY: clean test all

all: $(patsubst %, output/figure%.svg, $(flist))

output/figure%.svg: ckine/figures/figure%.py
	mkdir -p ./output
	poetry run fbuild $*

clean:
	rm -rf output

test:
	poetry run pytest -s -v -x
