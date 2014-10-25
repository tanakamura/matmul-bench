
#ifdef __AVX__
#define AVX_OP(I,J)                                                     \
    vout##I##_##J = _mm256_add_ps(vout##I##_##J, _mm256_mul_ps(lik##I##_8, vr##J)); \

#define AVX_FUNC_NAME matmul_block_outer_avx_omp
#include "matmul-bench-avxfunc.h"

#endif


#ifdef __FMA__

#define AVX_OP(I,J)                                                  \
    vout##I##_##J = _mm256_fmadd_ps(lik##I##_8, vr##J, vout##I##_##J); \

#define AVX_FUNC_NAME matmul_x86_fma
#include "matmul-bench-avxfunc.h"

#endif


