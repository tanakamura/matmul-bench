#ifndef MATMUL_BENCH_COMMON_H
#define MATMUL_BENCH_COMMON_H

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

#ifdef __SSE__
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

void matmulbench_init_simple_c(struct npr_varray *test_set);
void matmulbench_init_opt_c(struct npr_varray *test_set);
void matmulbench_init_sse(struct npr_varray *test_set);
void matmulbench_init_avx(struct npr_varray *test_set);
void matmulbench_init_fma(struct npr_varray *test_set);
void matmulbench_init_neon(struct npr_varray *test_set);

#endif