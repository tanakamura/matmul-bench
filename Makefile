all: all2

ANDROID_TOOLCHAIN=${HOME}/a/android-ndk-r10c/toolchains/arm-linux-androideabi-4.9/prebuilt/linux-x86_64/bin
ANDROID_PLATFORM=${HOME}/a/android-ndk-r10c/platforms/android-21/arch-arm/usr

WARN_FLAGS=-Wall -Werror=missing-prototypes -Werror=implicit-function-declaration # -Werror=unknown-pragmas

ifeq ($(OS),Windows_NT)
CFLAGS_COMMON= -std=gnu99 $(WARN_FLAGS) -O2 -ffast-math -falign-loops -MD -fvisibility=hidden -D MATMUL_BENCH_BUILD_LIB -I$(CURDIR)

WHICH=$(shell which2 which2)

ifeq ($(WHICH),)

WHICH := $(shell gcc -o which2.exe which.c -lshlwapi)

endif

rm2.exe: rm.c
	gcc -o $@ $<

RM=rm2.exe
WHICH=which2.exe

else
CFLAGS_COMMON= -pthread -std=gnu99 $(WARN_FLAGS) -O2 -ffast-math -falign-loops -MD -fvisibility=hidden -D MATMUL_BENCH_BUILD_LIB -fPIC -I$(CURDIR)

WHICH=which
RM=/bin/rm

endif

SAVE_TEMPS=1
ifdef SAVE_TEMPS
	CFLAGS_COMMON+=-save-temps=obj
endif

ARM_ANDROID_GCC=${ANDROID_TOOLCHAIN}/arm-linux-androideabi-gcc
ARM_LINUX_GCC=arm-linux-gnueabihf-gcc
X86_64_LINUX_GCC=x86_64-linux-gnu-gcc
X86_64_MINGW64_GCC=x86_64-w64-mingw32-gcc
X86_MINGW32_GCC=i686-w64-mingw32-gcc
X86_64_MINGW64_DLLTOOL=x86_64-w64-mingw32-dlltool
X86_MINGW32_DLLTOOL=i686-w64-mingw32-dlltool

ALL_TARGET=

ifneq ("$(shell $(WHICH) ${X86_64_LINUX_GCC})","")
ALL_TARGET+=matmul-bench-x86_64-linux dll/x86_64/libmatmul-bench.so
endif

ifneq ("$(shell $(WHICH) ${X86_64_MINGW64_GCC})","")
ALL_TARGET+=matmul-bench-w64.exe dll/w64/matmul-bench.dll
endif

ifneq ("$(shell $(WHICH) ${X86_MINGW32_GCC})","")
ALL_TARGET+=matmul-bench-w32.exe dll/w32/matmul-bench.dll
endif

ifneq ("$(shell $(WHICH) ${ARM_LINUX_GCC})","")
ALL_TARGET+=matmul-bench-arm-linux dll/arm-linux/libmatmul-bench.so
endif

ifneq ("$(shell $(WHICH) ${ARM_ANDROID_GCC})","")
ALL_TARGET+=matmul-bench-arm-android dll/arm-android/libmatmul-bench.so
endif

all2: ${ALL_TARGET}

NPR_SRCS=varray.c mempool-c.c
LIBBENCH_SRCS= \
	matmul-bench-simple-c.c \
	matmul-bench-opt-c.c \
	matmul-bench.c \
	$(NPR_SRCS)
BENCH_SRCS=$(LIBBENCH_SRCS)

X86_SRCS_BASE=${BENCH_SRCS} matmul-bench-sse.c matmul-bench-avx.c matmul-bench-fma.c
X86_EXE_SRCS=$(patsubst %.c,$(CURDIR)/%.c,$(X86_SRCS_BASE) matmul-bench-main.c)
X86_LIB_SRCS=$(patsubst %.c,$(CURDIR)/%.c,$(X86_SRCS_BASE))

X86_EXE_OBJS=$(patsubst %.c,$(CURDIR)/obj/x86_64/%.o,${X86_SRCS_BASE} matmul-bench-main.c)
X86_LIB_OBJS=$(patsubst %.c,$(CURDIR)/obj/x86_64/%.o,${X86_SRCS_BASE})

LINUX_LDLIBS=-lrt

W64_EXE_OBJS=$(patsubst %.c,$(CURDIR)/obj/w64/%.o,${X86_SRCS_BASE} matmul-bench-main.c)
W64_LIB_OBJS=$(patsubst %.c,$(CURDIR)/obj/w64/%.o,${X86_SRCS_BASE})

W32_EXE_OBJS=$(patsubst %.c,$(CURDIR)/obj/w32/%.o,${X86_SRCS_BASE} matmul-bench-main.c)
W32_LIB_OBJS=$(patsubst %.c,$(CURDIR)/obj/w32/%.o,${X86_SRCS_BASE})

ARM_SRCS_BASE=${BENCH_SRCS} matmul-bench-neon.c matmul-bench-vfpv4.c
ARM_EXE_SRCS=$(patsubst %.c,$(CURDIR)/%.c,$(ARM_SRCS_BASE) matmul-bench-main.c)
ARM_LIB_SRCS=$(patsubst %.c,$(CURDIR)/%.c,$(ARM_SRCS_BASE))

ARM_LINUX_EXE_OBJS=$(patsubst %.c,$(CURDIR)/obj/arm-linux/%.o,$(ARM_SRCS_BASE) matmul-bench-main.c)
ARM_ANDROID_EXE_OBJS=$(patsubst %.c,$(CURDIR)/obj/arm-android/%.o,$(ARM_SRCS_BASE) matmul-bench-main.c)
ARM_LINUX_LIB_OBJS=$(patsubst %.c,$(CURDIR)/obj/arm-linux/%.o,$(ARM_SRCS_BASE))
ARM_ANDROID_LIB_OBJS=$(patsubst %.c,$(CURDIR)/obj/arm-android/%.o,$(ARM_SRCS_BASE))

ALL_OBJS=$(X86_EXE_OBJS) $(W64_EXE_OBJS) $(W32_EXE_OBJS) $(ARM_LINUX_EXE_OBJS) $(ARM_ANDROID_EXE_OBJS)
ALL_SRCS=$(X86_EXE_SRCS) $(ARM_EXE_SRCS)

ALL_ASMS=$(ALL_SRCS:.c=.s)
ALL_PPS=$(ALL_SRCS:.c=.i)
ALL_DEPS=$(ALL_OBJS:.o=.d)

CFLAGS_ANDROID=$(CFLAGS_COMMON) -I${ANDROID_PLATFORM}/include -L${ANDROID_PLATFORM}/lib --sysroot ${ANDROID_PLATFORM}
CFLAGS_W32=-msse2 $(CFLAGS_COMMON) -static 


matmul-bench-x86_64-linux: $(X86_EXE_OBJS)
	${X86_64_LINUX_GCC} ${LINUX_LDLIBS} ${CFLAGS_COMMON} -o $@ $^
dll/x86_64/libmatmul-bench.so: $(X86_LIB_OBJS)
	${X86_64_LINUX_GCC} ${LINUX_LDLIBS} -shared ${CFLAGS_COMMON} -o $@ $^

matmul-bench-w64.exe: $(W64_EXE_OBJS)
	${X86_64_MINGW64_GCC} ${CFLAGS_COMMON} -static -o $@ $^
dll/w64/matmul-bench.dll: $(W64_LIB_OBJS)
	${X86_64_MINGW64_GCC} -Wl,--out-implib=dll/w64/matmul-bench.lib -static -s -shared ${CFLAGS_COMMON} -o $@ $^

matmul-bench-w32.exe: $(W32_EXE_OBJS)
	${X86_MINGW32_GCC} ${CFLAGS_W32} -o $@ $^
dll/w32/matmul-bench.dll: $(W32_LIB_OBJS)
	${X86_MINGW32_GCC} -Wl,--out-implib=dll/w32/matmul-bench.lib -static -s -shared ${CFLAGS_W32} -o $@ $^

matmul-bench-arm-linux: $(ARM_LINUX_EXE_OBJS)
	${ARM_LINUX_GCC} ${LINUX_LDLIBS} ${CFLAGS_COMMON} -o $@ $^
dll/arm-linux/libmatmul-bench.so: $(ARM_LINUX_LIB_OBJS)
	${ARM_LINUX_GCC} ${LINUX_LDLIBS} ${CFLAGS_COMMON} -shared -o $@ $^

matmul-bench-arm-android: $(ARM_ANDROID_EXE_OBJS)
	${ARM_ANDROID_GCC} ${CFLAGS_ANDROID} -o $@ $^
dll/arm-android/libmatmul-bench.so: $(ARM_ANDROID_LIB_OBJS)
	${ARM_ANDROID_GCC} ${CFLAGS_ANDROID} -shared -o $@ $^


$(CURDIR)/obj/x86_64/%.o: $(CURDIR)/npr/%.c
	${X86_64_LINUX_GCC} ${CFLAGS_COMMON} -c -o $@ $<
$(CURDIR)/obj/w64/%.o: $(CURDIR)/npr/%.c
	${X86_64_MINGW64_GCC} ${CFLAGS_COMMON} -c -o $@ $<
$(CURDIR)/obj/w32/%.o: $(CURDIR)/npr/%.c
	${X86_MINGW32_GCC} ${CFLAGS_W32} -c -o $@ $<
$(CURDIR)/obj/arm-linux/%.o: $(CURDIR)/npr/%.c
	${ARM_LINUX_GCC} ${CFLAGS_COMMON} -c -o $@ $<
$(CURDIR)/obj/arm-android/%.o: $(CURDIR)/npr/%.c
	${ARM_ANDROID_GCC} ${CFLAGS_ANDROID} -c -o $@ $<


$(CURDIR)/obj/x86_64/%.o: $(CURDIR)/%.c
	${X86_64_LINUX_GCC} ${CFLAGS_COMMON} -c -o $@ $<
$(CURDIR)/obj/w64/%.o: $(CURDIR)/%.c
	${X86_64_MINGW64_GCC} ${CFLAGS_COMMON} -c -o $@ $<
$(CURDIR)/obj/w32/%.o: $(CURDIR)/%.c
	${X86_MINGW32_GCC} ${CFLAGS_W32} -c -o $@ $<
$(CURDIR)/obj/arm-linux/%.o: $(CURDIR)/%.c
	${ARM_LINUX_GCC} ${CFLAGS_COMMON} -c -o $@ $<
$(CURDIR)/obj/arm-android/%.o: $(CURDIR)/%.c
	${ARM_ANDROID_GCC} ${CFLAGS_ANDROID} -c -o $@ $<


$(CURDIR)/obj/x86_64/matmul-bench-avx.o: $(CURDIR)/matmul-bench-avx.c
	${X86_64_LINUX_GCC} ${CFLAGS_COMMON} -mavx -c -o $@ $<
$(CURDIR)/obj/x86_64/matmul-bench-fma.o: $(CURDIR)/matmul-bench-fma.c
	${X86_64_LINUX_GCC} ${CFLAGS_COMMON} -march=native -mtune=native -mfma -c -o $@ $<
$(CURDIR)/obj/w64/matmul-bench-avx.o: $(CURDIR)/matmul-bench-avx.c
	${X86_64_MINGW64_GCC} ${CFLAGS_COMMON} -mavx -c -o $@ $<
$(CURDIR)/obj/w64/matmul-bench-fma.o: $(CURDIR)/matmul-bench-fma.c
	${X86_64_MINGW64_GCC} ${CFLAGS_COMMON} -march=native -mtune=native -mfma -c -o $@ $<
$(CURDIR)/obj/w32/matmul-bench-avx.o: $(CURDIR)/matmul-bench-avx.c
	${X86_MINGW32_GCC} ${CFLAGS_W32} -mavx -c -o $@ $<
$(CURDIR)/obj/w32/matmul-bench-fma.o: $(CURDIR)/matmul-bench-fma.c
	${X86_MINGW32_GCC} ${CFLAGS_W32} -march=native -mtune=native -mfma -c -o $@ $<

$(CURDIR)/obj/arm-linux/matmul-bench-neon.o: $(CURDIR)/matmul-bench-neon.c
	${ARM_LINUX_GCC} ${CFLAGS_COMMON} -mfloat-abi=hard -mfpu=neon -c -o $@ $<
$(CURDIR)/obj/arm-android/matmul-bench-neon.o: $(CURDIR)/matmul-bench-neon.c
	${ARM_ANDROID_GCC} ${CFLAGS_ANDROID} -mfloat-abi=softfp -mfpu=neon -c -o $@ $<

$(CURDIR)/obj/arm-linux/matmul-bench-vfpv4.o: $(CURDIR)/matmul-bench-vfpv4.c
	${ARM_LINUX_GCC} ${CFLAGS_COMMON} -mfloat-abi=hard -mfpu=neon-vfpv4 -c -o $@ $<
$(CURDIR)/obj/arm-android/matmul-bench-vfpv4.o: $(CURDIR)/matmul-bench-vfpv4.c
	${ARM_ANDROID_GCC} ${CFLAGS_ANDROID} -mfloat-abi=softfp -mfpu=neon-vfpv4 -c -o $@ $<


clean: $(RM)
	$(RM) -f $(ALL_OBJS) $(ALL_ASMS) $(ALL_DEPS) $(ALL_PPS) $(ALL_TARGET)

-include $(ALL_OBJS:.o=.d)
