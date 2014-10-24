#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

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


#define M8(M)                                   \
    M(0) M(1) M(2) M(3) M(4) M(5) M(6) M(7)

#define M16(M)                                  \
    M(0) M(1) M(2) M(3) M(4) M(5) M(6) M(7)     \
    M(8) M(9) M(10) M(11) M(12) M(13) M(14) M(15)

int n;

static float *in0;
static float *in1;
static float *out_simple;
static float *out_simple_outer_omp;
static float *out_simple_omp;
static float *out_block_omp;
static float *out_block_omp_unroll;
static float *out_block_outer_sse_omp;
static float *out_block_outer_avx_omp;
static float *out_x86_fma;
static float *out_block_outer_neon_omp;
static float *out_neon;

#define STRINGIZE_(a) #a
#define STRINGIZE(a) STRINGIZE_(a)

#define CONCAT_(a,b) a ## b
#define CONCAT(a,b) CONCAT_(a,b)

#ifdef _WIN32
#include <windows.h>

LARGE_INTEGER freq;

void
sec_init()
{
    QueryPerformanceFrequency(&freq);
}

double
sec(void)
{
    LARGE_INTEGER c;
    QueryPerformanceCounter(&c);

    return c.QuadPart/ (double)freq.QuadPart;

}

double
drand(void)
{
    unsigned int v;
    rand_s(&v);

    return v / (double)UINT_MAX;
}

#define srand(v)



#else
#define sec_init()

double
sec(void)
{
    struct timespec ts;

#ifdef CLOCK_MONOTONIC_RAW
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
#else
    clock_gettime(CLOCK_MONOTONIC, &ts);
#endif

    return (ts.tv_sec) + (ts.tv_nsec / (1000.0*1000.0*1000.0));
}

#define drand drand48
#define srand srand48

#define _aligned_malloc(sz,a) memalign(a,sz)
#define _aligned_free(p) free(p)

#endif

static void
matmul_simple(float *__restrict out,
              const float * __restrict inL,
              const float * __restrict inR,
              int n)
{
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            float v = 0;
            for (int k=0; k<n; k++) {
                v += inL[i*n+k] * inR[k*n+j];
            }

            out[i*n+j] = v;
        }
    }

    return;
}

static void
matmul_simple_outer_omp(float *__restrict out,
                        const float * __restrict inL,
                        const float * __restrict inR,
                        int n)
{
    int i;

#pragma omp parallel for
    for (i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            out[i*n+j] = 0;
        }
    }

#pragma omp parallel for
    for (i=0; i<n; i++) {
        for (int k=0; k<n; k++) {
            float lik = inL[i*n+k];

            for (int j=0; j<n; j++) {
                out[i*n+j] += lik * inR[k*n + j];
            }
        }
    }

    return;
}

static void
matmul_simple_omp(float * __restrict out,
                  const float * __restrict inL,
                  const float * __restrict inR,
                  int n)
{
    int i;
#pragma omp parallel for
    for (i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            float v = 0;
            for (int k=0; k<n; k++) {
                v += inL[i*n+k] * inR[k*n+j];
            }
            out[i*n+j] = v;
        }
    }

    return;
}

static void
matmul_block_omp(float * __restrict out,
                 const float* __restrict inL,
                 const float* __restrict inR,
                 unsigned int n)
{
    unsigned int block_size = 32;
    int i0, i;

#pragma omp parallel for
    for (i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            out[i*n+j] = 0;
        }
    }



#pragma omp parallel for
    for (i0=0; i0<n; i0+=block_size) {
        for (int j0=0; j0<n; j0+=block_size) {
            for (int k0=0; k0<n; k0+=block_size) {
                for (int bi=0; bi<block_size; bi++) {
                    int i = i0+bi;

                    for (int bk=0; bk<block_size; bk++) {
                        int k = k0+bk;

                        float lik = inL[i*n + k];

                        for (int bj=0; bj<block_size; bj++) {
                            int j = j0+bj;
                            out[i*n + j] += lik * inR[k*n + j];
                        }
                    }
                }
            }
        }
    }
}

#ifdef __SSE__
static void
matmul_block_outer_sse_omp(float * __restrict out,
                           const float* __restrict inL,
                           const float* __restrict inR,
                           unsigned int n)
{
    unsigned int block_size = 32;
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

#define OUTER_SSE_J(J)                                               \
                        long j_##J = (J*4);                           \
                        __m128 vr##J = _mm_load_ps(&inRp[j_##J]); \
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
            }
        }
    }
}
#endif

#ifdef __AVX__
#define AVX_OP(I,J)                                                     \
    vout##I##_##J = _mm256_add_ps(vout##I##_##J, _mm256_mul_ps(lik##I##_8, vr##J)); \

#define AVX_FUNC_NAME matmul_block_outer_avx_omp
#include "avxfunc.h"

#endif


#ifdef __FMA__

#define AVX_OP(I,J)                                                  \
    vout##I##_##J = _mm256_fmadd_ps(lik##I##_8, vr##J, vout##I##_##J); \

#define AVX_FUNC_NAME matmul_x86_fma
#include "avxfunc.h"

#endif



#ifdef __ARM_NEON__

static void
matmul_neon(float * __restrict out,
            const float* __restrict inL,
            const float* __restrict inR,
            unsigned int n)
{
    /* C=4x4x(2simd) register */
    unsigned int block_size_i = 4;
    unsigned int block_size_j = 8;
    unsigned int block_size_k = 64;
    int i00;

#pragma omp parallel for schedule(dynamic)
    for (i00=0; i00<n; i00+=block_size_i) {
        for (int j0=0; j0<n; j0+=block_size_j) {
            for (int k0=0; k0<n; k0+=block_size_k) {
                for (int bi=0; bi<block_size_i; bi+=4) {
                    int i0 = i00+bi+0;

                    float *outp0 = &out[i0*n+j0];

#define outp1 (outp0+n*1)
#define outp2 (outp0+n*2)
#define outp3 (outp0+n*3)

                    float32x2_t *outp2_0 = (float32x2_t*)outp0;
                    float32x2_t *outp2_1 = (float32x2_t*)outp1;
                    float32x2_t *outp2_2 = (float32x2_t*)outp2;
                    float32x2_t *outp2_3 = (float32x2_t*)outp3;

                    __builtin_prefetch(outp0);
                    __builtin_prefetch(outp1);
                    __builtin_prefetch(outp2);
                    __builtin_prefetch(outp3);

                    /* 4x4x(2simd) reg */
                    float32x2_t vout0_0;
                    float32x2_t vout0_1;
                    float32x2_t vout0_2;
                    float32x2_t vout0_3;

                    float32x2_t vout1_0;
                    float32x2_t vout1_1;
                    float32x2_t vout1_2;
                    float32x2_t vout1_3;

                    float32x2_t vout2_0;
                    float32x2_t vout2_1;
                    float32x2_t vout2_2;
                    float32x2_t vout2_3;

                    float32x2_t vout3_0;
                    float32x2_t vout3_1;
                    float32x2_t vout3_2;
                    float32x2_t vout3_3;

                    vout0_0 = vdup_n_f32(0);
                    vout0_1 = vdup_n_f32(0);
                    vout0_2 = vdup_n_f32(0);
                    vout0_3 = vdup_n_f32(0);
                    vout1_0 = vdup_n_f32(0);
                    vout1_1 = vdup_n_f32(0);
                    vout1_2 = vdup_n_f32(0);
                    vout1_3 = vdup_n_f32(0);
                    vout2_0 = vdup_n_f32(0);
                    vout2_1 = vdup_n_f32(0);
                    vout2_2 = vdup_n_f32(0);
                    vout2_3 = vdup_n_f32(0);
                    vout3_0 = vdup_n_f32(0);
                    vout3_1 = vdup_n_f32(0);
                    vout3_2 = vdup_n_f32(0);
                    vout3_3 = vdup_n_f32(0);

                    const float *inRp1 = (float*)&inR[k0*n+j0];

                    const float32x2_t *inRp;

                    float32x2_t vr0;
                    float32x2_t vr1;
                    float32x2_t vr2;
                    float32x2_t vr3;

                    const float *inL00 = (inL + i0*n + k0);

#define x_vmlaq_n_f32(a, b, c, n)                                       \
    __asm__ __volatile__ ("vmla.f32 %P[A], %P[B], %P[C][" STRINGIZE(n) "]\n\t" \
                          :[A]"+w"(a)                                   \
                          :[B]"w"(b), [C]"x"(c));

#define NEON_K(K)                                                       \
                    {                                                   \
                        inRp = (float32x2_t*)(inRp1);                   \
                        __builtin_prefetch(inRp1 + n*2);                \
                        vr0 = inRp[0];                                  \
                        vr1 = inRp[1];                                  \
                        vr2 = inRp[2];                                  \
                        vr3 = inRp[3];                                  \
                                                                        \
                        lik = vdup_n_f32(inL00[0*n+bk+K]);              \
                        vout0_0 = vmla_f32(vout0_0, vr0, lik);          \
                        vout0_1 = vmla_f32(vout0_1, vr1, lik);          \
                        vout0_2 = vmla_f32(vout0_2, vr2, lik);          \
                        vout0_3 = vmla_f32(vout0_3, vr3, lik);          \
                                                                        \
                        lik = vdup_n_f32(inL00[1*n+bk+K]);              \
                        vout1_0 = vmla_f32(vout1_0, vr0, lik);          \
                        vout1_1 = vmla_f32(vout1_1, vr1, lik);          \
                        vout1_2 = vmla_f32(vout1_2, vr2, lik);          \
                        vout1_3 = vmla_f32(vout1_3, vr3, lik);          \
                                                                        \
                        lik = vdup_n_f32(inL00[2*n+bk+K]);              \
                        vout2_0 = vmla_f32(vout2_0, vr0, lik);          \
                        vout2_1 = vmla_f32(vout2_1, vr1, lik);          \
                        vout2_2 = vmla_f32(vout2_2, vr2, lik);          \
                        vout2_3 = vmla_f32(vout2_3, vr3, lik);          \
                                                                        \
                        lik = vdup_n_f32(inL00[3*n+bk+K]);              \
                        vout3_0 = vmla_f32(vout3_0, vr0, lik);          \
                        vout3_1 = vmla_f32(vout3_1, vr1, lik);          \
                        vout3_2 = vmla_f32(vout3_2, vr2, lik);          \
                        vout3_3 = vmla_f32(vout3_3, vr3, lik);          \
                                                                        \
                        inRp1 += n;                                     \
                    }

                    for (int bk=0; bk<block_size_k; bk+=2) {
                        float32x2_t lik;
                        NEON_K(0);
                        NEON_K(1);
                    }

                    if (k0==0) {
                        outp2_0[0] = vout0_0;
                        outp2_0[1] = vout0_1;
                        outp2_0[2] = vout0_2;
                        outp2_0[3] = vout0_3;

                        outp2_1[0] = vout1_0;
                        outp2_1[1] = vout1_1;
                        outp2_1[2] = vout1_2;
                        outp2_1[3] = vout1_3;

                        outp2_2[0] = vout2_0;
                        outp2_2[1] = vout2_1;
                        outp2_2[2] = vout2_2;
                        outp2_2[3] = vout2_3;

                        outp2_3[0] = vout3_0;
                        outp2_3[1] = vout3_1;
                        outp2_3[2] = vout3_2;
                        outp2_3[3] = vout3_3;
                    } else {
                        outp2_0[0] = vadd_f32(outp2_0[0], vout0_0);
                        outp2_0[1] = vadd_f32(outp2_0[1], vout0_1);
                        outp2_0[2] = vadd_f32(outp2_0[2], vout0_2);
                        outp2_0[3] = vadd_f32(outp2_0[3], vout0_3);

                        outp2_1[0] = vadd_f32(outp2_1[0], vout1_0);
                        outp2_1[1] = vadd_f32(outp2_1[1], vout1_1);
                        outp2_1[2] = vadd_f32(outp2_1[2], vout1_2);
                        outp2_1[3] = vadd_f32(outp2_1[3], vout1_3);

                        outp2_2[0] = vadd_f32(outp2_2[0], vout2_0);
                        outp2_2[1] = vadd_f32(outp2_2[1], vout2_1);
                        outp2_2[2] = vadd_f32(outp2_2[2], vout2_2);
                        outp2_2[3] = vadd_f32(outp2_2[3], vout2_3);

                        outp2_3[0] = vadd_f32(outp2_3[0], vout3_0);
                        outp2_3[1] = vadd_f32(outp2_3[1], vout3_1);
                        outp2_3[2] = vadd_f32(outp2_3[2], vout3_2);
                        outp2_3[3] = vadd_f32(outp2_3[3], vout3_3);
                    }

                }
            }
        }
    }
}
#endif



static void
matmul_block_omp_unroll(float * __restrict out,
                        const float* __restrict inL,
                        const float* __restrict inR,
                        unsigned int n)
{
    unsigned int block_size = 16;
    int i0;

#pragma omp parallel for schedule(dynamic)
    for (i0=0; i0<n; i0+=block_size) {
        for (int j0=0; j0<n; j0+=block_size) {
            for (int bi=0; bi<block_size; bi++) {
                for (int bj=0; bj<block_size; bj++) {
                    int i = i0+bi;
                    int j = j0+bj;
                    out[i*n+j] = 0;
                }
            }

            for (int k0=0; k0<n; k0+=block_size) {
                for (int bi=0; bi<block_size; bi+=4) {
                    int i_0 = i0+bi + 0;
                    int i_1 = i0+bi + 1;
                    int i_2 = i0+bi + 2;
                    int i_3 = i0+bi + 3;

                    for (int bj=0; bj<block_size; bj++) {
                        float v_0 = 0;
                        float v_1 = 0;
                        float v_2 = 0;
                        float v_3 = 0;

                        int j = j0+bj;

                        const float *lp_0 = inL + i_0*n + k0;
                        const float *rp = inR + k0*n + j;

                        for (int bk=0; bk<block_size; bk++) {
                            //int k = k0+bk;

                            v_0 += lp_0[0*n] * *rp;
                            v_1 += lp_0[1*n] * *rp;
                            v_2 += lp_0[2*n] * *rp;
                            v_3 += lp_0[3*n] * *rp;

                            lp_0 += 1;

                            rp += n;
                        }

                        out[i_0*n+j] += v_0;
                        out[i_1*n+j] += v_1;
                        out[i_2*n+j] += v_2;
                        out[i_3*n+j] += v_3;
                    }
                }
            }
        }
    }
}

static void
dump_mat(int n, float *data)
{
    if (n <= 4) {
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                printf("%8.1f, ", data[i*n + j]);
            }

            printf("\n");
        }
    }
}

static void
dump_flops(const char *tag,
           int ni,
           double sec,
           float *data,
           int iter)
{
    double n = ni;
    printf("%d:%-20s: sec=%8.5f, %8.5f[GFLOPS], %8.5f[GB/s]\n",
           iter,
           tag,
           sec,
           n*n*n*2/(sec*1024.0*1024.0*1024.0),
           (n*n*3.0*sizeof(float))/(sec*1024.0*1024.0*1024.0));

    dump_mat(ni, data);

    if (data != out_simple) {
        for (int i=0; i<n*n; i++) {
            float delta = fabs(data[i]-out_simple[i]);
            double ratio = (delta/fabs(data[i]))*100;

            if (ratio > 1e-3) {
                printf("error delta=%e(%e[%%]), simple=%e, opt=%e\n", delta, ratio, data[i], out_simple[i]);
                exit(1);
            }
        }
    }

}

static void
bench(int iter)
{
    double t0, t1;

    if (n < 512) {
        t0 = sec();
        matmul_simple(out_simple, in0, in1, n);
        t1 = sec();

        dump_flops("simple(gold)", n, t1-t0, out_simple, iter);


        t0 = sec();
        matmul_simple_outer_omp(out_simple_outer_omp, in0, in1, n);
        t1 = sec();
        dump_flops("simple_outer_omp", n, t1-t0, out_simple_outer_omp, iter);

        t0 = sec();
        matmul_simple_omp(out_simple_omp, in0, in1, n);
        t1 = sec();

        dump_flops("simple_omp", n, t1-t0, out_simple_omp, iter);
    } else {
        t0 = sec();
        matmul_simple_outer_omp(out_simple, in0, in1, n);
        t1 = sec();

        dump_flops("outer_omp(gold)", n, t1-t0, out_simple, iter);
    }


    if (n > 16) {
        t0 = sec();
        matmul_block_omp(out_block_omp, in0, in1, n);
        t1 = sec();

        dump_flops("block_omp", n, t1-t0, out_block_omp, iter);

        if (n > 64) {
            t0 = sec();
            matmul_block_omp_unroll(out_block_omp_unroll, in0, in1, n);
            t1 = sec();

            dump_flops("block_omp_unroll", n, t1-t0, out_block_omp_unroll, iter);


#ifdef __SSE__
            t0 = sec();
            matmul_block_outer_sse_omp(out_block_outer_sse_omp, in0, in1, n);
            t1 = sec();

            dump_flops("sse", n, t1-t0, out_block_outer_sse_omp, iter);
#endif

#ifdef __AVX__
            t0 = sec();
            matmul_block_outer_avx_omp(out_block_outer_avx_omp, in0, in1, n);
            t1 = sec();

            dump_flops("avx", n, t1-t0, out_block_outer_avx_omp, iter);
#endif



#ifdef __FMA__
            t0 = sec();
            matmul_x86_fma(out_x86_fma, in0, in1, n);
            t1 = sec();

            dump_flops("fma", n, t1-t0, out_x86_fma, iter);
#endif


#ifdef __ARM_NEON__
            t0 = sec();
            matmul_neon(out_neon, in0, in1, n);
            t1 = sec();

            dump_flops("neon", n, t1-t0, out_neon, iter);
#endif
        }
    }
}


int
main(int argc, char **argv)
{
    int iter = 3;
    sec_init();

    n = 512;

    if (argc >= 2) {
        n = atoi(argv[1]);
    }

#ifdef _OPENMP
    if (argc >= 3) {
        omp_set_num_threads(atoi(argv[2]));
    }
#endif

    if (argc >= 4) {
        iter = atoi(argv[3]);
    }

    int align = 64;

    in0 = _aligned_malloc(n*n * sizeof(float), align);
    in1 = _aligned_malloc(n*n * sizeof(float), align);

    out_simple = _aligned_malloc(n*n * sizeof(float), align);
    out_simple_outer_omp = _aligned_malloc(n*n * sizeof(float), align);
    out_simple_omp = _aligned_malloc(n*n * sizeof(float), align);
    out_block_omp = _aligned_malloc(n*n * sizeof(float), align);
    out_block_omp_unroll = _aligned_malloc(n*n * sizeof(float), align);
    out_block_outer_sse_omp = _aligned_malloc(n*n * sizeof(float), align);
    out_block_outer_avx_omp = _aligned_malloc(n*n * sizeof(float), align);
    out_block_outer_neon_omp = _aligned_malloc(n*n * sizeof(float), align);
    out_x86_fma = _aligned_malloc(n*n * sizeof(float), align);
    out_neon = _aligned_malloc(n*n * sizeof(float), align);

    printf("n=%d\n", n);
    printf("total nelem=%d, mat size=%f[MB]\n",
           n*n,
           (n*n)/(1024.0*1024.0)*sizeof(float));

    srand(100);

    for (int i=0; i<n*n; i++) {
        in0[i] = drand()+1.0;
        in1[i] = drand()+1.0;
    }

    dump_mat(n, in0);
    dump_mat(n, in1);

    for (int i=0; i<iter; i++) {
        bench(i);
    }

    _aligned_free(in0);
    _aligned_free(in1);
    _aligned_free(out_simple);
    _aligned_free(out_simple_outer_omp);
    _aligned_free(out_simple_omp);
    _aligned_free(out_block_omp);
    _aligned_free(out_block_omp_unroll);
    _aligned_free(out_block_outer_sse_omp);
    _aligned_free(out_block_outer_avx_omp);
    _aligned_free(out_block_outer_neon_omp);
    _aligned_free(out_neon);
    _aligned_free(out_x86_fma);

}
