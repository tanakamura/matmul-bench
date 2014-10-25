#include "matmul-bench-common.h"

#ifdef __AVX__
#define AVX_OP(I,J)                                                     \
    vout##I##_##J = _mm256_add_ps(vout##I##_##J, _mm256_mul_ps(lik##I##_8, vr##J)); \

#define AVX_FUNC_NAME avx_run
#include "matmul-bench-avxfunc.h"
#endif


static const struct MatmulBenchTest avx = MATMULBENCH_TEST_INITIALIZER("avx", avx_run, 192);

void
matmulbench_init_avx(struct MatmulBench *b, struct npr_varray *test_set)
{
    VA_PUSH(struct MatmulBenchTest, test_set, avx);
}
