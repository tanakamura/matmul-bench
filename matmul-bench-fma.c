#if 0

#ifdef __FMA__

#define AVX_OP(I,J)                                                  \
    vout##I##_##J = _mm256_fmadd_ps(lik##I##_8, vr##J, vout##I##_##J); \

#define AVX_FUNC_NAME matmul_x86_fma
#include "matmul-bench-avxfunc.h"

#endif

#endif