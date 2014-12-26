all: w32/matmul-bench.exe w64/matmul-bench.exe w32/matmul-bench.dll w64/matmul-bench.dll

CFLAGS_COMMON=-I.. -MD  -Wall -std=c99 -Werror=implicit-function-declaration
CFLAGS_W32=$(CFLAGS_COMMON) -L../dll/w32
CFLAGS_W64=$(CFLAGS_COMMON) -L../dll/w64

SRCS=main.c strbuf.c mempool-c.c
OBJS_32=$(patsubst %.c,$(CURDIR)/w32/%.o,$(SRCS))
OBJS_64=$(patsubst %.c,$(CURDIR)/w64/%.o,$(SRCS))

LDLIBS=-lmatmul-bench -luser32 -lgdi32 -lcomdlg32
LDFLAGS= # -mwindows

w32/matmul-bench.exe: $(OBJS_32) w32/resource.res
	${X86_MINGW32_GCC} ${CFLAGS_W32} -o $@ $^ $(LDLIBS) $(LDFLAGS)

w64/matmul-bench.exe: $(OBJS_64) w64/resource.res
	${X86_64_MINGW64_GCC} ${CFLAGS_W64} -o $@ $^ $(LDLIBS) $(LDFLAGS)

w32/matmul-bench.dll: ..\dll\w32\matmul-bench.dll
	copy ..\dll\w32\matmul-bench.dll w32

w64/matmul-bench.dll: ..\dll\w64\matmul-bench.dll
	copy ..\dll\w64\matmul-bench.dll w64

w32/resource.res: resource.rc
	windres -O coff -i $< -o $@ -F pe-i386

w64/resource.res: resource.rc
	windres -O coff -i $< -o $@ -F pe-x86-64

X86_64_MINGW64_GCC=x86_64-w64-mingw32-gcc
X86_MINGW32_GCC=i686-w64-mingw32-gcc

$(CURDIR)/w32/%.o: $(CURDIR)/%.c
	${X86_MINGW32_GCC} ${CFLAGS_W32} -c -o $@ $<
$(CURDIR)/w64/%.o: $(CURDIR)/%.c
	${X86_64_MINGW64_GCC} ${CFLAGS_W64} -c -o $@ $<

$(CURDIR)/w32/%.o: $(CURDIR)/../npr/%.c
	${X86_MINGW32_GCC} ${CFLAGS_W32} -c -o $@ $<
$(CURDIR)/w64/%.o: $(CURDIR)/../npr/%.c
	${X86_64_MINGW64_GCC} ${CFLAGS_W64} -c -o $@ $<

ALLOBJS=$(OBJS_32) $(OBJS_64)

clean:
	../rm2 $(ALLOBJS)  w32/resource.res w64/resource.res w32/matmul-bench.exe w64/matmul-bench.exe

-include $(ALL_OBJS:.o=.d)
