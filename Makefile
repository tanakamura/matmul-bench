all: matmul-bench
CFLAGS=-std=gnu99 -Wall -O2 -fopenmp -ffast-math -g -mfloat-abi=hard -mfpu=neon -mcpu=cortex-a9 # -mavx
LDFLAGS=-g
CC=arm-linux-gnueabihf-gcc

matmul-bench: matmul-bench.c

clean:
	rm -f matmul-bench
