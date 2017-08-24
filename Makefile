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

find_executable = $(firstword $(foreach a, $(1), $(shell $(WHICH) ${a} 2>/dev/null)))

SAVE_TEMPS=1
ifdef SAVE_TEMPS
	CFLAGS_COMMON+=-save-temps=obj
endif

ARM_ANDROID_GCC=${ANDROID_TOOLCHAIN}/arm-linux-androideabi-gcc
ARM_LINUX_GCC=arm-linux-gnueabihf-gcc
AARCH64_LINUX_GCC=$(call find_executable, aarch64-linux-gnu-gcc aarch64-unknown-linux-gnu-gcc)
X86_64_LINUX_GCC=x86_64-linux-gnu-gcc
X86_64_MINGW64_GCC=x86_64-w64-mingw32-gcc
X86_MINGW32_GCC=i686-w64-mingw32-gcc
X86_64_MINGW64_DLLTOOL=x86_64-w64-mingw32-dlltool
X86_MINGW32_DLLTOOL=i686-w64-mingw32-dlltool

ALL_TARGET=

NPR_SRCS=varray.c mempool-c.c
LIBBENCH_SRCS= \
	matmul-bench-simple-c.c \
	matmul-bench-opt-c.c \
	matmul-bench.c \
	$(NPR_SRCS)
BENCH_SRCS=$(LIBBENCH_SRCS)

X86_SIMD_SRCS=matmul-bench-sse.c matmul-bench-avx.c matmul-bench-fma.c
ARM_SIMD_SRCS=matmul-bench-neon.c matmul-bench-vfpv4.c

LINUX_LDLIBS=-lrt

# CFLAGS_ANDROID=$(CFLAGS_COMMON) -I${ANDROID_PLATFORM}/include -L${ANDROID_PLATFORM}/lib --sysroot ${ANDROID_PLATFORM}

ALL_SRCS=${BENCH_SRCS} matmul-bench-main.c

define genarch # $(1):arch-name  $(2):cross-compiler $(3):additional-src $(4):additional-ldlibs

ifneq ("$(shell $(WHICH) ${2} 2>/dev/null)","")

$(1)_EXE_OBJS=$(patsubst %.c,$(CURDIR)/obj/$(1)/%.o,${BENCH_SRCS} $(3) matmul-bench-main.c)
$(1)_LIB_OBJS=$(patsubst %.c,$(CURDIR)/obj/$(1)/%.o,${BENCH_SRCS} $(3))

$(CURDIR)/obj/$(1)/keep:
	mkdir -p $(CURDIR)/obj/$(1); touch $$@

$(CURDIR)/dll/$(1)/keep:
	mkdir -p $(CURDIR)/dll/$(1); touch $$@

$(1)_DIRS=$(CURDIR)/obj/$(1) $(CURDIR)/dll/$(1)
$(1)_DIR_DEPS=$(CURDIR)/obj/$(1)/keep $(CURDIR)/dll/$(1)/keep

$(CURDIR)/obj/$(1)/%.o: $(CURDIR)/npr/%.c $${$(1)_DIR_DEPS}
	$(2) ${CFLAGS_COMMON} $${$${*}_CFLAGS} -c -o $$@ $$<

$(CURDIR)/obj/$(1)/%.o: $(CURDIR)/%.c $${$(1)_DIR_DEPS}
	$(2) ${CFLAGS_COMMON} $${$${*}_CFLAGS} -c -o $$@ $$<

matmul-bench-$(1): $$($(1)_EXE_OBJS) ${$(1)_DIR_DEPS}
	$(2) $(4) ${CFLAGS_COMMON} -o $$@ $$($(1)_EXE_OBJS)

dll/$(1)/libmatmul-bench.so: $${$(1)_LIB_OBJS} ${$(1)_DIR_DEPS}
	$(2) $(4) -shared ${CFLAGS_COMMON} -o $$@ $${$(1)_LIB_OBJS}

ALL_OBJS+=$${$(1)_EXE_OBJS}
CLEAN_TARGETS+=$(CURDIR)/obj/$(1) $(CURDIR)/dll/$(1)
ALL_TARGET+=matmul-bench-$(1) dll/$(1)/libmatmul-bench.so
ALL_SRCS+=$(3)

endif

endef

$(eval $(call genarch,x86_64-linux,$(X86_64_LINUX_GCC), $(X86_SIMD_SRCS), $(LINUX_LDLIBS)))
$(eval $(call genarch,w64,${X86_64_MINGW64_GCC}, $(X86_SIMD_SRCS), ))
$(eval $(call genarch,w32,${X86_MINGW32_GCC}, $(X86_SIMD_SRCS), ))
$(eval $(call genarch,arm-linux,${ARM_LINUX_GCC}, $(ARM_SIMD_SRCS), ))
$(eval $(call genarch,aarch64-linux,$(AARCH64_LINUX_GCC), , $(LINUX_LDLIBS)))

all2: ${ALL_TARGET}

matmul-bench-avx_CFLAGS=-mavx
matmul-bench-fma_CFLAGS=-mfma
matmul-bench-sse_CFLAGS=-msse
matmul-bench-neon_CFLAGS=-mfloat-abi=hard -mfpu=neon
matmul-bench-vfpv4_CFLAGS=-mfloat-abi=hard -mfpu=neon-vfpv4 

# xx
#$(CURDIR)/obj/arm-android/matmul-bench-neon.o: $(CURDIR)/matmul-bench-neon.c
#	${ARM_ANDROID_GCC} ${CFLAGS_ANDROID} -mfloat-abi=softfp -mfpu=neon -c -o $@ $<
#$(CURDIR)/obj/arm-android/matmul-bench-vfpv4.o: $(CURDIR)/matmul-bench-vfpv4.c
#	${ARM_ANDROID_GCC} ${CFLAGS_ANDROID} -mfloat-abi=softfp -mfpu=neon-vfpv4 -c -o $@ $<


clean: $(RM)
	$(RM) -rf obj dll $(ALL_TARGET)

-include $(ALL_OBJS:.o=.d)

