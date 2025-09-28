.PHONY: all

model = AKAP79

all: C/$(model)_gvf.c R/$(model).R


R/%.R: %.tar.gz
	/Users/fedmil/Documents/githubPackages/RPN-derivative/sh/ode.sh -R --maxima $^ > $@

C/%_gvf.c: %.tar.gz
	/Users/fedmil/Documents/githubPackages/RPN-derivative/sh/ode.sh -C --maxima $^ > $@

vfgen/%.vf %.tar.gz %.zip: *.tsv
	/Users/fedmil/Documents/githubPackages/SBtabVFGEN/sbtab_to_vfgen $^


