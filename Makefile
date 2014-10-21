all: matmul-bench
ifeq (${CROSS}, 1)
CFLAGS=-std=gnu99 -Wall -O2 -ffast-math -mfloat-abi=hard -fopenmp -mfpu=neon -mtune=cortex-a9 -mcpu=cortex-a9 -save-temps
CC=arm-linux-gcc
else ifeq (${ANDROID}, 1)
CFLAGS=-std=gnu99 -Wall -O2 -ffast-math -mfloat-abi=softfp -mfpu=neon -save-temps
CC=arm-linux-androideabi-gcc
else
CFLAGS=-std=gnu99 -Wall -O2 -fopenmp -ffast-math -mtune=native -save-temps -march=atom
endif

LDFLAGS=-g $(CFLAGS)

matmul-bench.o: matmul-bench.c
#matmul-bench.o: matmul-bench.s
#	$(CC) -o $@ $< $(CFLAGS) -c

matmul-bench: matmul-bench.o

clean:
	rm -f matmul-bench *.o *.s *.i *~
