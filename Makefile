CC=gcc
CFLAGS=-std=c99 -O2

OBJECTS = depore

all: $(OBJECTS)

depore: src/main.c
	$(CC) $(CFLAGS) src/main.c src/aln.c -o depore -lz -lm

.PHONY: clean
clean:
	-rm $(OBJECTS)
