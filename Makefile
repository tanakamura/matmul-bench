ANDROID_TOOLCHAIN=${HOME}/a/android-ndk-r10c/toolchains/arm-linux-androideabi-4.9/prebuilt/linux-x86_64/bin
ANDROID_PLATFORM=${HOME}/a/android-ndk-r10c/platforms/android-21/arch-arm/usr

all: matmul-bench
ifeq (${CROSS}, 1)
CFLAGS=-std=gnu99 -Wall -O2 -ffast-math -mfloat-abi=hard -fopenmp -mfpu=neon -mtune=cortex-a9 -mcpu=cortex-a9 -falign-loops=16 -save-temps
CC=arm-linux-gnueabihf-gcc
else ifeq (${ANDROID}, 1)
CFLAGS=-std=gnu99 -Wall -O2 -ffast-math -mfloat-abi=softfp -save-temps -fopenmp -I${ANDROID_PLATFORM}/include -L${ANDROID_PLATFORM}/lib --sysroot ${ANDROID_PLATFORM}
CC=${ANDROID_TOOLCHAIN}/arm-linux-androideabi-gcc
else
CFLAGS=-std=gnu99 -Wall -O2 -fopenmp -ffast-math -mtune=native -save-temps -march=native
endif

LDFLAGS=-g $(CFLAGS)

matmul-bench.o: matmul-bench.c avxfunc.h
#matmul-bench.o: matmul-bench.s
#	$(CC) -o $@ $< $(CFLAGS) -c

matmul-bench: matmul-bench.o

clean:
	rm -f matmul-bench *.o *.s *.i *~
