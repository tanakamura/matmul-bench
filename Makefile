all: matmul-bench
CFLAGS=-std=gnu99 -Wall -O2 -fopenmp -ffast-math -g # -mavx
LDFLAGS=-g

matmul-bench: matmul-bench.c

clean:
	rm -f matmul-bench
