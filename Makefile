.PHONY: all

model = AKAP79

all: C/$(model)_gvf.c R/$(model).R


R/%.R: %.tar.gz
	ode -R --maxima $^ > $@

C/%_gvf.c: %.tar.gz
	ode -C --maxima $^ > $@

vfgen/%.vf %.tar.gz %.zip: *.tsv
	sbtab_to_vfgen $^


