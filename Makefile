ANDROID_TOOLCHAIN=${HOME}/a/android-ndk-r10c/toolchains/arm-linux-androideabi-4.9/prebuilt/linux-x86_64/bin
ANDROID_PLATFORM=${HOME}/a/android-ndk-r10c/platforms/android-21/arch-arm/usr

CFLAGS_COMMON=-std=gnu99 -Wall -O2 -fopenmp -ffast-math -falign-loops -MD -fvisibility=hidden -I$(PWD)

ifdef SAVE_TEMPS
	CFLAGS_COMMON+=-save-temps
endif

ARM_ANDROID_GCC=${ANDROID_TOOLCHAIN}/arm-linux-androideabi-gcc
ARM_LINUX_GCC=arm-linux-gnueabihf-gcc
X86_64_LINUX_GCC=x86_64-linux-gnu-gcc

ALL_TARGET=

ifneq ("$(shell which ${X86_64_LINUX_GCC})","")
ALL_TARGET+=matmul-bench-x86_64-linux
endif

ifneq ("$(shell which ${ARM_LINUX_GCC})","")
ALL_TARGET+=matmul-bench-arm-linux
endif

ifneq ("$(shell which ${ARM_ANDROID_GCC})","")
ALL_TARGET+=matmul-bench-arm-android
endif

all: ${ALL_TARGET}

NPR_SRCS=varray.c mempool-c.c
LIBBENCH_SRCS= \
	matmul-bench-simple-c.c \
	matmul-bench-opt-c.c \
	matmul-bench.c \
	$(NPR_SRCS)
BENCH_SRCS=matmul-bench-main.c  $(LIBBENCH_SRCS)

X86_SRCS_BASE=${BENCH_SRCS} matmul-bench-sse.c matmul-bench-avx.c matmul-bench-fma.c
X86_SRCS=$(patsubst %.c,$(PWD)/%.c,$(X86_SRCS_BASE))

X86_OBJS=$(patsubst %.c,$(PWD)/obj/x86_64/%.o,${X86_SRCS_BASE})

ARM_SRCS_BASE=${BENCH_SRCS} matmul-bench-neon.c matmul-bench-vfpv4.c
ARM_SRCS=$(patsubst %.c,$(PWD)/%.c,$(ARM_SRCS_BASE=$))

ARM_LINUX_OBJS=$(patsubst %.c,$(PWD)/obj/arm-linux/%.o,${ARM_SRCS_BASE})
ARM_ANDROID_OBJS=$(patsubst %.c,$(PWD)/obj/arm-android/%.o,${ARM_SRCS_BASE})

ALL_OBJS=$(X86_OBJS) $(ARM_LINUX_OBJS) $(ARM_ANDROID_OBJS)
ALL_SRCS=$(X86_SRCS) $(ARM_SRCS)

ALL_ASMS=$(ALL_SRCS:.c=.s)
ALL_PPS=$(ALL_SRCS:.c=.i)
ALL_DEPS=$(ALL_OBJS:.o=.d)

CFLAGS_ANDROID=$(CFLAGS_COMMON) -I${ANDROID_PLATFORM}/include -L${ANDROID_PLATFORM}/lib --sysroot ${ANDROID_PLATFORM}

matmul-bench-x86_64-linux: $(X86_OBJS)
	${X86_64_LINUX_GCC} ${CFLAGS_COMMON} -o $@ $^

matmul-bench-arm-linux: $(ARM_LINUX_OBJS)
	${ARM_LINUX_GCC} ${CFLAGS_COMMON} -o $@ $^

matmul-bench-arm-android: $(ARM_ANDROID_OBJS)
	${ARM_ANDROID_GCC} ${CFLAGS_ANDROID} -o $@ $^


$(PWD)/obj/x86_64/%.o: $(PWD)/npr/%.c
	cd obj/x86_64; ${X86_64_LINUX_GCC} ${CFLAGS_COMMON} -c -o $@ $<
$(PWD)/obj/arm-linux/%.o: $(PWD)/npr/%.c
	cd obj/arm-linux; ${ARM_LINUX_GCC} ${CFLAGS_COMMON} -c -o $@ $<
$(PWD)/obj/arm-android/%.o: $(PWD)/npr/%.c
	cd obj/arm-android; ${ARM_ANDROID_GCC} ${CFLAGS_ANDROID} -c -o $@ $<


$(PWD)/obj/x86_64/%.o: $(PWD)/%.c
	cd obj/x86_64; ${X86_64_LINUX_GCC} ${CFLAGS_COMMON} -c -o $@ $<
$(PWD)/obj/arm-linux/%.o: $(PWD)/%.c
	cd obj/arm-linux; ${ARM_LINUX_GCC} ${CFLAGS_COMMON} -c -o $@ $<
$(PWD)/obj/arm-android/%.o: $(PWD)/%.c
	cd obj/arm-android; ${ARM_ANDROID_GCC} ${CFLAGS_ANDROID} -c -o $@ $<


$(PWD)/obj/x86_64/matmul-bench-avx.o: $(PWD)/matmul-bench-avx.c
	cd obj/x86_64; ${X86_64_LINUX_GCC} ${CFLAGS_COMMON} -mavx -c -o $@ $<
$(PWD)/obj/x86_64/matmul-bench-fma.o: $(PWD)/matmul-bench-fma.c
	cd obj/x86_64; ${X86_64_LINUX_GCC} ${CFLAGS_COMMON} -mfma -c -o $@ $<

$(PWD)/obj/arm-linux/matmul-bench-neon.o: $(PWD)/matmul-bench-neon.c
	cd obj/arm-linux; ${ARM_LINUX_GCC} ${CFLAGS_COMMON} -mfloat-abi=hard -mfpu=neon -c -o $@ $<
$(PWD)/obj/arm-android/matmul-bench-neon.o: $(PWD)/matmul-bench-neon.c
	cd obj/arm-android; ${ARM_ANDROID_GCC} ${CFLAGS_ANDROID} -mfloat-abi=softfp -mfpu=neon -c -o $@ $<

$(PWD)/obj/arm-linux/matmul-bench-vfpv4.o: $(PWD)/matmul-bench-vfpv4.c
	cd obj/arm-linux; ${ARM_LINUX_GCC} ${CFLAGS_COMMON} -mfloat-abi=hard -mfpu=neon-vfpv4 -c -o $@ $<
$(PWD)/obj/arm-android/matmul-bench-vfpv4.o: $(PWD)/matmul-bench-vfpv4.c
	cd obj/arm-android; ${ARM_ANDROID_GCC} ${CFLAGS_ANDROID} -mfloat-abi=softfp -mfpu=neon-vfpv4 -c -o $@ $<


clean:
	rm -f $(ALL_OBJS) $(ALL_ASMS) $(ALL_DEPS) $(ALL_PPS) ${ALL_TARGET}

-include $(ALL_OBJS:.o=.d)
