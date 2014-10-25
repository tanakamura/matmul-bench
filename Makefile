#ANDROID_TOOLCHAIN=${HOME}/a/android-ndk-r10c/toolchains/arm-linux-androideabi-4.9/prebuilt/linux-x86_64/bin
#ANDROID_PLATFORM=${HOME}/a/android-ndk-r10c/platforms/android-21/arch-arm/usr
#
#all: matmul-bench
#ifeq (${CROSS}, 1)
#CFLAGS=-std=gnu99 -Wall -O2 -ffast-math -mfloat-abi=hard -fopenmp -mfpu=neon -mtune=cortex-a9 -mcpu=cortex-a9 -falign-loops=16 -save-temps
#CC=arm-linux-gnueabihf-gcc
#else ifeq (${ANDROID}, 1)
#CFLAGS=-std=gnu99 -Wall -O2 -ffast-math -mfloat-abi=softfp -mfpu=neon-vfpv4 -mtune=cortex-a15 -mcpu=cortex-a15 -save-temps -fopenmp -I${ANDROID_PLATFORM}/include -L${ANDROID_PLATFORM}/lib --sysroot ${ANDROID_PLATFORM}
#CC=${ANDROID_TOOLCHAIN}/arm-linux-androideabi-gcc
#else
#CFLAGS=-std=gnu99 -Wall -O2 -fopenmp -ffast-math -mtune=native -save-temps -march=native
#endif
#
#LDFLAGS=-g $(CFLAGS)
#
#matmul-bench.o: matmul-bench.c avxfunc.h
##matmul-bench.o: matmul-bench.s
##	$(CC) -o $@ $< $(CFLAGS) -c
#
#matmul-bench: matmul-bench.o
#
#clean:
#	rm -f matmul-bench *.o *.s *.i *~


ANDROID_TOOLCHAIN=${HOME}/a/android-ndk-r10c/toolchains/arm-linux-androideabi-4.9/prebuilt/linux-x86_64/bin
ANDROID_PLATFORM=${HOME}/a/android-ndk-r10c/platforms/android-21/arch-arm/usr

CFLAGS_COMMON=-std=gnu99 -Wall -O2 -fopenmp -ffast-math -MD -fvisibility=hidden -I$(PWD)

ifdef SAVE_TEMPS
	CFLAGS_COMMON+=-save-temps
endif

ARM_ANDROID_GCC=${ANDROID_TOOLCHAIN}/arm-linux-androideabi-gcc
ARM_LINUX_GCC=arm-linux-gnueabihf-gcc
X86_64_LINUX_GCC=x86_64-linux-gnu-gcc

ALL_TARGET=
TARGET_ARCHS=arm-android arm-linux x86_64

ifneq ("$(shell which ${X86_64_LINUX_GCC})","")
HAVE_X86_64_LINUX=1
ALL_TARGET+=matmul-bench-x86_64-linux
TARGET_ARCHS+=x86_64
endif

ifneq ("$(shell which ${ARM_LINUX_GCC})","")
HAVE_ARM_LINUX=1
ALL_TARGET+=matmul-bench-arm-linux
TARGET_ARCHS+=arm-linux
endif

all: ${ALL_TARGET}

NPR_SRCS=varray.c mempool-c.c
LIBBENCH_SRCS=matmul-bench-simple-c.c matmul-bench.c $(NPR_SRCS)
BENCH_SRCS=matmul-bench-main.c  $(LIBBENCH_SRCS)

X86_OBJS=$(patsubst %.c,obj/x86_64/%.o,${BENCH_SRCS})
ARM_LINUX_OBJS=$(patsubst %.c,obj/arm-linux/%.o,${BENCH_SRCS})
ARM_ANDROID_OBJS=$(patsubst %.c,obj/arm-android/%.o,${BENCH_SRCS})

ALL_OBJS=$(X86_OBJS) $(ARM_LINUX_OBJS) $(ARM_ANDROID_OBJS)
ALL_ASMS=$(ALL_OBJS:.o=.s)
ALL_DEPS=$(ALL_OBJS:.o=.d)
ALL_PPS=$(ALL_OBJS:.o=.i)



matmul-bench-x86_64-linux: $(X86_OBJS)
	${X86_64_LINUX_GCC} ${CFLAGS_COMMON} -o $@ $^

matmul-bench-arm-linux: $(ARM_LINUX_OBJS)
	${ARM_LINUX_GCC} ${CFLAGS_COMMON} -o $@ $^


obj/x86_64/%.o: npr/%.c
	${X86_64_LINUX_GCC} ${CFLAGS_COMMON} -c -o $@ $<

obj/arm-linux/%.o: npr/%.c
	${ARM_LINUX_GCC} ${CFLAGS_COMMON} -c -o $@ $<

obj/x86_64/%.o: %.c
	${X86_64_LINUX_GCC} ${CFLAGS_COMMON} -c -o $@ $<

obj/arm-linux/%.o: %.c
	${ARM_LINUX_GCC} ${CFLAGS_COMMON} -c -o $@ $<

clean:
	rm -f $(ALL_OBJS) $(ALL_ASMS) $(ALL_DEPS) $(ALL_PPS)

-include $(ALL_OBJS:.o=.d)
