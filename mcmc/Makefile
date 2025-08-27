CC = gcc
CFLAGS = -march=native -O2 -fPIC -Wfatal-errors -Wall `pkg-config --cflags gsl`
LIBS = `pkg-config --libs gsl`

.PHONY: all


all: AKAP79.so


%.so: ../C/%_gvf.c
	$(CC) -shared $(CFLAGS) -o $@ $^ $(LIBS)
