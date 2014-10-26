
#include <emmintrin.h>

#include "matmul-bench-common.h"

static NOINLINE void
sse_(unsigned long i0,
     unsigned long j0,
     unsigned long k0,
     unsigned long bi,
     float * __restrict out,
     const float* __restrict inL,
     const float* __restrict inR,
     unsigned long n)
{
    unsigned int block_size = 32;
    int i = i0+bi;
    float *outp = &out[i*n+j0];
    __m128 *outp4 = (__m128*)outp;

    __m128 vout0 = _mm_setzero_ps();
    __m128 vout1 = _mm_setzero_ps();
    __m128 vout2 = _mm_setzero_ps();
    __m128 vout3 = _mm_setzero_ps();
    __m128 vout4 = _mm_setzero_ps();
    __m128 vout5 = _mm_setzero_ps();
    __m128 vout6 = _mm_setzero_ps();
    __m128 vout7 = _mm_setzero_ps();

    _mm_prefetch((const char*)(outp + n) ,_MM_HINT_T0);

    for (long bk=0; bk<block_size; bk++) {
        long k = k0+bk;

        const float *inRp = &inR[k*n+j0];

        _mm_prefetch((const char*)(inRp + n),_MM_HINT_T0);

        float lik = inL[i*n+k];
        __m128 lik4 = _mm_set1_ps(lik);

#define OUTER_SSE_J(J)                                                  \
        long j_##J = (J*4);                                             \
            __m128 vr##J = _mm_load_ps(&inRp[j_##J]);                   \
                vout##J = _mm_add_ps(vout##J, _mm_mul_ps(lik4, vr##J)); \

        OUTER_SSE_J(0);
        OUTER_SSE_J(1);
        OUTER_SSE_J(2);
        OUTER_SSE_J(3);

        OUTER_SSE_J(4);
        OUTER_SSE_J(5);
        OUTER_SSE_J(6);
        OUTER_SSE_J(7);
    }

    outp4[0] = _mm_add_ps(outp4[0], vout0);
    outp4[1] = _mm_add_ps(outp4[1], vout1);
    outp4[2] = _mm_add_ps(outp4[2], vout2);
    outp4[3] = _mm_add_ps(outp4[3], vout3);

    outp4[4] = _mm_add_ps(outp4[4], vout4);
    outp4[5] = _mm_add_ps(outp4[5], vout5);
    outp4[6] = _mm_add_ps(outp4[6], vout6);
    outp4[7] = _mm_add_ps(outp4[7], vout7);
    
}


static void
sse_run(float * __restrict out,
        const float * __restrict inL,
        const float * __restrict inR,
        const float * __restrict inL_plus1line,
        const float * __restrict inR_plus1line,
        unsigned int n,
        unsigned int pitch_byte)
{
    const unsigned int block_size = 32;
    int i0;

#pragma omp parallel for schedule(dynamic)
    for (i0=0; i0<n; i0+=block_size) {
        for (int j0=0; j0<n; j0+=block_size) {
            for (int bi=0; bi<block_size; bi++) {
                for (int bj=0; bj<block_size; bj+=4) {
                    int i = i0+bi;
                    int j = j0+bj;
                    _mm_store_ps(&out[i*n+j], _mm_setzero_ps());
                }
            }

            for (int k0=0; k0<n; k0+=block_size) {
                for (int bi=0; bi<block_size; bi++) {
                    sse_(i0, j0, k0, bi, out, inL, inR, n);
                }
            }
        }
    }
}

static const struct MatmulBenchTest sse = MATMULBENCH_TEST_INITIALIZER("sse", sse_run, 32);

void
matmulbench_init_sse(struct MatmulBench *b, struct npr_varray *test_set)
{
    VA_PUSH(struct MatmulBenchTest, test_set, sse);
}
