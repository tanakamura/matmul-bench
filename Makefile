all: matmul-bench
ifeq (${CROSS}, 1)
CFLAGS=-std=gnu99 -Wall -O2 -fopenmp -ffast-math -mfloat-abi=hard -mfpu=neon -mcpu=cortex-a9 -save-temps
CC=arm-linux-gnueabihf-gcc
else
CFLAGS=-std=gnu99 -Wall -O2 -fopenmp -ffast-math -mtune=native -save-temps -march=native
#endif
endif

LDFLAGS=-g $(CFLAGS)

matmul-bench.o: matmul-bench.c
matmul-bench: matmul-bench.o

clean:
	rm -f matmul-bench *.o *.s
