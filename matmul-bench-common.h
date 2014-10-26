#ifndef MATMUL_BENCH_COMMON_H
#define MATMUL_BENCH_COMMON_H

#if defined __x86_64__  || defined __x86__
#define ARCH_X86
#endif

#if defined __GNUC__ || defined __clang__
#define HAVE_VEC_EXT
#endif

#include "matmul-bench.h"
#include "npr/varray.h"

struct MatmulBench;

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __GNUC__
#define NOINLINE __attribute__((noinline,noclone))
#else
#define NOINLINE __declspec(noinline)
#endif

#ifdef ARCH_X86
#ifdef _WIN32
#include <intrin.h>
#else
#include <x86intrin.h>
#endif
#endif

#ifdef __ARM_NEON__
#include <arm_neon.h>
#endif

#ifdef _WIN32

#else
#include <malloc.h>
#define _aligned_malloc(sz,a) memalign(a,sz)
#define _aligned_free(p) free(p)
#endif

#define STRINGIZE_(a) #a
#define STRINGIZE(a) STRINGIZE_(a)

#define CONCAT_(a,b) a ## b
#define CONCAT(a,b) CONCAT_(a,b)



#define M8(M)                                   \
    M(0) M(1) M(2) M(3) M(4) M(5) M(6) M(7)

#define M16(M)                                  \
    M(0) M(1) M(2) M(3) M(4) M(5) M(6) M(7)     \
    M(8) M(9) M(10) M(11) M(12) M(13) M(14) M(15)

double matmul_bench_sec(void);

void matmulbench_init_simple_c(struct MatmulBench *b, struct npr_varray *test_set);
void matmulbench_init_opt_c(struct MatmulBench *b, struct npr_varray *test_set);
void matmulbench_init_sse(struct MatmulBench *b, struct npr_varray *test_set);
void matmulbench_init_avx(struct MatmulBench *b, struct npr_varray *test_set);
void matmulbench_init_fma(struct MatmulBench *b, struct npr_varray *test_set);
void matmulbench_init_neon(struct MatmulBench *b, struct npr_varray *test_set);
void matmulbench_init_vfpv4(struct MatmulBench *b, struct npr_varray *test_set);

#define MATMULBENCH_TEST_INITIALIZER(name,run,size_step) {name, run, size_step}
#define CEIL_DIV(a,b) (((a)+((b)-1))/(b))

#endif
