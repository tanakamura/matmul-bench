#include <immintrin.h>
#include "matmul-bench-common.h"

#define MAT_4x2

#define AVX_OP(I,J)                                                  \
    vout##I##_##J = _mm256_fmadd_ps(lik##I##_8, vr##J, vout##I##_##J); \

#define AVX_FUNC_NAME fma_run
#include "matmul-bench-avxfunc.h"

static void
fma_tmp_run(struct MatmulBenchParam *p)
{
    unsigned long n = p->n;

    float * __restrict out = p->out;
    const float * __restrict inL = p->inL;
    const float * __restrict inR = p->inR;

    float *tR_base = _aligned_malloc(8*3*n*sizeof(float), 64);

    /*     r0  r1  r2
     *    +-----------
     * l0 |o00 o01 o02
     * l1 |o10 o11 o12  4 elem
     * l2 |o20 o21 o22
     * l3 |o30 o31 o32
     *       24 elem
     */

#define O_BLKSIZE_HOR (8U * 3U)
#define O_BLKSIZE_VER 4U

    unsigned long o_bhi, o_bvi, r_vi, l_hi;
    unsigned long num_oblk_hor = n / O_BLKSIZE_HOR;
    unsigned long num_oblk_ver = n / O_BLKSIZE_VER;

    __m256 l0, l1, l2, l3;
    __m256 r0, r1, r2;

    for (o_bhi=0; o_bhi<num_oblk_hor; o_bhi++) {
        __m256 o00=_mm256_setzero_ps(), o01=_mm256_setzero_ps(), o02=_mm256_setzero_ps();
        __m256 o10=_mm256_setzero_ps(), o11=_mm256_setzero_ps(), o12=_mm256_setzero_ps();
        __m256 o20=_mm256_setzero_ps(), o21=_mm256_setzero_ps(), o22=_mm256_setzero_ps();
        __m256 o30=_mm256_setzero_ps(), o31=_mm256_setzero_ps(), o32=_mm256_setzero_ps();

        float *pl0 = inL;
        const float * __restrict Rline = inR + o_bhi * O_BLKSIZE_HOR;
        __m256 *tR = (__m256*)tR_base;

#define pl1 pl0+n*1
#define pl2 pl0+n*2
#define pl3 pl0+n*3

        for (l_hi=0; l_hi<n; l_hi++) {
            r0 = _mm256_loadu_ps(Rline + 0);
            r1 = _mm256_loadu_ps(Rline + 8);
            r2 = _mm256_loadu_ps(Rline + 16);

            Rline += n;

            l0 = _mm256_broadcast_ss(pl0);
            o00 = _mm256_fmadd_ps(l0, r0, o00);
            o01 = _mm256_fmadd_ps(l0, r1, o01);
            o02 = _mm256_fmadd_ps(l0, r2, o02);

            l1 = _mm256_broadcast_ss(pl1);
            o10 = _mm256_fmadd_ps(l1, r0, o10);
            o11 = _mm256_fmadd_ps(l1, r1, o11);
            o12 = _mm256_fmadd_ps(l1, r2, o12);

            l2 = _mm256_broadcast_ss(pl2);
            o20 = _mm256_fmadd_ps(l2, r0, o20);
            o21 = _mm256_fmadd_ps(l2, r1, o21);
            o22 = _mm256_fmadd_ps(l2, r2, o22);

            l3 = _mm256_broadcast_ss(pl3);
            o30 = _mm256_fmadd_ps(l3, r0, o30);
            o31 = _mm256_fmadd_ps(l3, r1, o31);
            o32 = _mm256_fmadd_ps(l3, r2, o32);

            *(tR++) = r0;
            *(tR++) = r1;
            *(tR++) = r2;

            pl0++;
        }

        float * __restrict out_base = out + o_bhi * O_BLKSIZE_HOR;

        _mm256_storeu_ps(out_base+n*0+0,  o00);
        _mm256_storeu_ps(out_base+n*0+8,  o01);
        _mm256_storeu_ps(out_base+n*0+16, o02);

        _mm256_storeu_ps(out_base+n*1+0,  o10);
        _mm256_storeu_ps(out_base+n*1+8,  o11);
        _mm256_storeu_ps(out_base+n*1+16, o12);

        _mm256_storeu_ps(out_base+n*2+0,  o20);
        _mm256_storeu_ps(out_base+n*2+8,  o21);
        _mm256_storeu_ps(out_base+n*2+16, o22);

        _mm256_storeu_ps(out_base+n*3+0,  o30);
        _mm256_storeu_ps(out_base+n*3+8,  o31);
        _mm256_storeu_ps(out_base+n*3+16, o32);

        for (o_bvi=1; o_bvi<num_oblk_ver; o_bvi++) {
            __m256 o00=_mm256_setzero_ps(), o01=_mm256_setzero_ps(), o02=_mm256_setzero_ps();
            __m256 o10=_mm256_setzero_ps(), o11=_mm256_setzero_ps(), o12=_mm256_setzero_ps();
            __m256 o20=_mm256_setzero_ps(), o21=_mm256_setzero_ps(), o22=_mm256_setzero_ps();
            __m256 o30=_mm256_setzero_ps(), o31=_mm256_setzero_ps(), o32=_mm256_setzero_ps();

            float *pl0 = inL + o_bvi * O_BLKSIZE_VER * n;
            const float * __restrict Rline = tR_base;

            for (l_hi=0; l_hi<n; l_hi++) {
                r0 = _mm256_loadu_ps(Rline + 0);
                r1 = _mm256_loadu_ps(Rline + 8);
                r2 = _mm256_loadu_ps(Rline + 16);

                Rline += 8*3;

                l0 = _mm256_broadcast_ss(pl0);
                o00 = _mm256_fmadd_ps(l0, r0, o00);
                o01 = _mm256_fmadd_ps(l0, r1, o01);
                o02 = _mm256_fmadd_ps(l0, r2, o02);

                l1 = _mm256_broadcast_ss(pl1);
                o10 = _mm256_fmadd_ps(l1, r0, o10);
                o11 = _mm256_fmadd_ps(l1, r1, o11);
                o12 = _mm256_fmadd_ps(l1, r2, o12);

                l2 = _mm256_broadcast_ss(pl2);
                o20 = _mm256_fmadd_ps(l2, r0, o20);
                o21 = _mm256_fmadd_ps(l2, r1, o21);
                o22 = _mm256_fmadd_ps(l2, r2, o22);

                l3 = _mm256_broadcast_ss(pl3);
                o30 = _mm256_fmadd_ps(l3, r0, o30);
                o31 = _mm256_fmadd_ps(l3, r1, o31);
                o32 = _mm256_fmadd_ps(l3, r2, o32);

                pl0++;
            }

            float * __restrict out_base = out + o_bhi * O_BLKSIZE_HOR + (o_bvi * O_BLKSIZE_VER * n);

            _mm256_storeu_ps(out_base+n*0+0,  o00);
            _mm256_storeu_ps(out_base+n*0+8,  o01);
            _mm256_storeu_ps(out_base+n*0+16, o02);

            _mm256_storeu_ps(out_base+n*1+0,  o10);
            _mm256_storeu_ps(out_base+n*1+8,  o11);
            _mm256_storeu_ps(out_base+n*1+16, o12);

            _mm256_storeu_ps(out_base+n*2+0,  o20);
            _mm256_storeu_ps(out_base+n*2+8,  o21);
            _mm256_storeu_ps(out_base+n*2+16, o22);

            _mm256_storeu_ps(out_base+n*3+0,  o30);
            _mm256_storeu_ps(out_base+n*3+8,  o31);
            _mm256_storeu_ps(out_base+n*3+16, o32);
        }
    }

    _aligned_free(tR_base);
}

static const struct MatmulBenchTest fma = MATMULBENCH_TEST_INITIALIZER("fma", fma_run, 128);
static const struct MatmulBenchTest fma_tmp = MATMULBENCH_TEST_INITIALIZER("fma_tmp", fma_tmp_run, 24);

void
matmulbench_init_fma(struct MatmulBench *b, struct npr_varray *test_set)
{
    VA_PUSH(struct MatmulBenchTest, test_set, fma);
    VA_PUSH(struct MatmulBenchTest, test_set, fma_tmp);
}

