#include "matmul-bench-common.h"

#ifdef __FMA__

#define AVX_OP(I,J)                                                  \
    vout##I##_##J = _mm256_fmadd_ps(lik##I##_8, vr##J, vout##I##_##J); \

#define AVX_FUNC_NAME fma_run
#include "matmul-bench-avxfunc.h"

#endif

static const struct MatmulBenchTest fma = MATMULBENCH_TEST_INITIALIZER("fma", fma_run, 192);

void
matmulbench_init_fma(struct MatmulBench *b, struct npr_varray *test_set)
{
    VA_PUSH(struct MatmulBenchTest, test_set, fma);
}

