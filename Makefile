.PHONY: all


model = AKAP79
# turn conservation law analysis on: --cla
OPT = --cla

all: C/$(model)_gvf.c R/$(model).R


C/%_gvf.c R/%.R: *.tsv
	./makeModel.R $(OPT)



