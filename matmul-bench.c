#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <omp.h>

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
static float *out_for_a9;

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

static void
matmul_block_omp_unroll_xx(float * __restrict out,
                           const float* __restrict inL,
                           const float* __restrict inR,
                           unsigned int n)
{
    unsigned int block_size = 32;
    int i0, i;

#pragma omp parallel for schedule(dynamic)
    for (i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            out[i*n+j] = 0;
        }
    }



#pragma omp parallel for schedule(dynamic)
    for (i0=0; i0<n; i0+=block_size) {
        for (int j0=0; j0<n; j0+=block_size) {
            for (int k0=0; k0<n; k0+=block_size) {
                for (int bi=0; bi<block_size; bi++) {
                    int i = i0+bi;

                    for (int bk=0; bk<block_size; bk++) {
                        int k = k0+bk;

                        float lik = inL[i*n + k];

                        for (int bj=0; bj<block_size; bj+=4) {
                            int j;

#define BLOCK_OMP_UNROLL_J(J)                                   \
                            j = j0+bj+J;                        \
                            out[i*n + j] += lik * inR[k*n + j]; \

                            BLOCK_OMP_UNROLL_J(0);
                            BLOCK_OMP_UNROLL_J(1);
                            BLOCK_OMP_UNROLL_J(2);
                            BLOCK_OMP_UNROLL_J(3);
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
    int i0, i;

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
static void
matmul_block_outer_avx_omp(float * __restrict out,
                           const float* __restrict inL,
                           const float* __restrict inR,
                           unsigned int n)
{
    unsigned int block_size = 32;
    int i0, i;

#pragma omp parallel for schedule(dynamic)
    for (i0=0; i0<n; i0+=block_size) {
        for (int j0=0; j0<n; j0+=block_size) {
            for (int bi=0; bi<block_size; bi++) {
                for (int bj=0; bj<block_size; bj+=8) {
                    int i = i0+bi;
                    int j = j0+bj;
                    _mm256_store_ps(&out[i*n+j], _mm256_setzero_ps());
                }
            }


            for (int k0=0; k0<n; k0+=block_size) {
                for (int bi=0; bi<block_size; bi++) {
                    int i = i0+bi;
                    float *outp = &out[i*n+j0];
                    __m256 *outp8 = (__m256*)outp;

                    __m256 vout0 = _mm256_setzero_ps();
                    __m256 vout1 = _mm256_setzero_ps();
                    __m256 vout2 = _mm256_setzero_ps();
                    __m256 vout3 = _mm256_setzero_ps();

                    _mm_prefetch((const char*)(outp + n) ,_MM_HINT_T0);

                    for (long bk=0; bk<block_size; bk++) {
                        long k = k0+bk;

                        const float *inRp = &inR[k*n+j0];

                        _mm_prefetch((const char*)(inRp + n), _MM_HINT_T0);

                        float lik = inL[i*n+k];
                        __m256 lik8 = _mm256_set1_ps(lik);

#define OUTER_AVX_J(J)                                               \
                        long j_##J = (J*8);                          \
                        __m256 vr##J = _mm256_load_ps(&inRp[j_##J]); \
                        vout##J = _mm256_add_ps(vout##J, _mm256_mul_ps(lik8, vr##J)); \

                        OUTER_AVX_J(0);
                        OUTER_AVX_J(1);
                        OUTER_AVX_J(2);
                        OUTER_AVX_J(3);
                    }

                    outp8[0] = _mm256_add_ps(outp8[0], vout0);
                    outp8[1] = _mm256_add_ps(outp8[1], vout1);
                    outp8[2] = _mm256_add_ps(outp8[2], vout2);
                    outp8[3] = _mm256_add_ps(outp8[3], vout3);
                }
            }
        }
    }
}
#endif


#ifdef __FMA__
static void
matmul_x86_fma(float * __restrict out,
               const float* __restrict inL,
               const float* __restrict inR,
               unsigned int n)
{
    unsigned int block_size = 32;
    int i0, i;

#pragma omp parallel for schedule(dynamic)
    for (i0=0; i0<n; i0+=block_size) {
        for (int j0=0; j0<n; j0+=block_size) {
            for (int bi=0; bi<block_size; bi++) {
                for (int bj=0; bj<block_size; bj+=8) {
                    int i = i0+bi;
                    int j = j0+bj;
                    _mm256_store_ps(&out[i*n+j], _mm256_setzero_ps());
                }
            }


            for (int k0=0; k0<n; k0+=block_size) {
                for (int bi=0; bi<block_size; bi++) {
                    int i = i0+bi;
                    float *outp = &out[i*n+j0];
                    __m256 *outp8 = (__m256*)outp;

                    __m256 vout0 = _mm256_setzero_ps();
                    __m256 vout1 = _mm256_setzero_ps();
                    __m256 vout2 = _mm256_setzero_ps();
                    __m256 vout3 = _mm256_setzero_ps();

                    _mm_prefetch((const char*)(outp + n) ,_MM_HINT_T0);

                    for (long bk=0; bk<block_size; bk++) {
                        long k = k0+bk;

                        const float *inRp = &inR[k*n+j0];

                        _mm_prefetch((const char*)(inRp + n), _MM_HINT_T0);

                        float lik = inL[i*n+k];
                        __m256 lik8 = _mm256_set1_ps(lik);

#define OUTER_AVX_J(J)                                               \
                        long j_##J = (J*8);                          \
                        __m256 vr##J = _mm256_load_ps(&inRp[j_##J]); \
                        vout##J = _mm256_fmadd_ps(lik8, vr##J, vout##J);     \

                        OUTER_AVX_J(0);
                        OUTER_AVX_J(1);
                        OUTER_AVX_J(2);
                        OUTER_AVX_J(3);
                    }

                    outp8[0] = _mm256_add_ps(outp8[0], vout0);
                    outp8[1] = _mm256_add_ps(outp8[1], vout1);
                    outp8[2] = _mm256_add_ps(outp8[2], vout2);
                    outp8[3] = _mm256_add_ps(outp8[3], vout3);
                }
            }
        }
    }
}
#endif



#ifdef __ARM_NEON__
static void
matmul_for_a9(float * __restrict out,
              const float* __restrict inL,
              const float* __restrict inR,
              unsigned int n)
{
    unsigned int block_size = 16;
    int i0, i;

#pragma omp parallel for
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
                for (int bi=0; bi<block_size; bi++) {
                    int i = i0+bi;

                    float *outp = &out[i*n+j0];
                    float32x2_t *outp2 = (float32x2_t*)outp;

                    __builtin_prefetch(outp + n);

                    float32x2_t vout0 = outp2[0];
                    float32x2_t vout1 = outp2[1];
                    float32x2_t vout2 = outp2[2];
                    float32x2_t vout3 = outp2[3];
                    float32x2_t vout4 = outp2[4];
                    float32x2_t vout5 = outp2[5];
                    float32x2_t vout6 = outp2[6];
                    float32x2_t vout7 = outp2[7];

                    const float *inRp1 = (float*)&inR[k0*n+j0];
                    const float *inL1 = inL + i*n + k0;

                    float32x2_t vr0, vr1, vr2, vr3, vr4, vr5, vr6, vr7;

                    {
                        const float32x2_t *inRp = (float32x2_t*)(inRp1);
                        __builtin_prefetch(inRp1 + n);
                        __builtin_prefetch(inRp1 + n*2);

                        vr0 = inRp[0];
                        vr1 = inRp[1];
                        vr2 = inRp[2];
                        vr3 = inRp[3];
                        vr4 = inRp[4];
                        vr5 = inRp[5];
                        vr6 = inRp[6];
                        vr7 = inRp[7];

                        inRp1 += n;
                    }

                    for (int bk=0; bk<block_size-1; bk++) {
                        float lik = *inL1;
                        float32x2_t lik2 = vdup_n_f32(lik);

                        const float32x2_t *inRp = (float32x2_t*)(inRp1);
                        __builtin_prefetch(inRp1 + n*2);

                        vout0 = vmla_f32(vout0, lik2, vr0);
                        vout1 = vmla_f32(vout1, lik2, vr1);
                        vout2 = vmla_f32(vout2, lik2, vr2);
                        vout3 = vmla_f32(vout3, lik2, vr3);
                        vout4 = vmla_f32(vout4, lik2, vr4);
                        vout5 = vmla_f32(vout5, lik2, vr5);
                        vout6 = vmla_f32(vout6, lik2, vr6);
                        vout7 = vmla_f32(vout7, lik2, vr7);

                        vr0 = inRp[0];
                        vr1 = inRp[1];
                        vr2 = inRp[2];
                        vr3 = inRp[3];
                        vr4 = inRp[4];
                        vr5 = inRp[5];
                        vr6 = inRp[6];
                        vr7 = inRp[7];

                        inRp1 += n;
                        inL1 ++;
                    }

                    {
                        float lik = *inL1;
                        float32x2_t lik2 = vdup_n_f32(lik);
                        vout0 = vmla_f32(vout0, lik2, vr0);
                        vout1 = vmla_f32(vout1, lik2, vr1);
                        vout2 = vmla_f32(vout2, lik2, vr2);
                        vout3 = vmla_f32(vout3, lik2, vr3);
                        vout4 = vmla_f32(vout4, lik2, vr4);
                        vout5 = vmla_f32(vout5, lik2, vr5);
                        vout6 = vmla_f32(vout6, lik2, vr6);
                        vout7 = vmla_f32(vout7, lik2, vr7);
                    }


                    outp2[0] = vout0;
                    outp2[1] = vout1;
                    outp2[2] = vout2;
                    outp2[3] = vout3;
                    outp2[4] = vout4;
                    outp2[5] = vout5;
                    outp2[6] = vout6;
                    outp2[7] = vout7;
                }
            }
        }
    }
}


static void
matmul_neon(float * __restrict out,
            const float* __restrict inL,
            const float* __restrict inR,
            unsigned int n)
{
    unsigned int block_size = 16;
    int i0, i;

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
                for (int bi=0; bi<block_size; bi++) {
                    int i = i0+bi;

                    float *outp = &out[i*n+j0];
                    float32x4_t *outp4 = (float32x4_t*)outp;

                    __builtin_prefetch(outp + n);

                    float32x4_t vout0 = outp4[0];
                    float32x4_t vout1 = outp4[1];
                    float32x4_t vout2 = outp4[2];
                    float32x4_t vout3 = outp4[3];

                    const float *inRp1 = (float*)&inR[k0*n+j0];
                    const float *inL1 = inL + i*n + k0;

                    for (int bk=0; bk<block_size; bk++) {
                        float lik = *inL1;
                        float32x4_t lik4 = vdupq_n_f32(lik);

                        const float32x4_t *inRp = (float32x4_t*)(inRp1);
                        __builtin_prefetch(inRp1 + n*2);

                        vout0 = vmlaq_f32(vout0, lik4, inRp[0]);
                        vout1 = vmlaq_f32(vout1, lik4, inRp[1]);
                        vout2 = vmlaq_f32(vout2, lik4, inRp[2]);
                        vout3 = vmlaq_f32(vout3, lik4, inRp[3]);

                        inRp1 += n;
                        inL1 ++;
                    }


                    outp4[0] = vout0;
                    outp4[1] = vout1;
                    outp4[2] = vout2;
                    outp4[3] = vout3;
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

    t0 = sec();
    matmul_simple(out_simple, in0, in1, n);
    t1 = sec();

    dump_flops("simple", n, t1-t0, out_simple, iter);


    t0 = sec();
    matmul_simple_outer_omp(out_simple_outer_omp, in0, in1, n);
    t1 = sec();

    dump_flops("simple_outer_omp", n, t1-t0, out_simple_outer_omp, iter);

    t0 = sec();
    matmul_simple_omp(out_simple_omp, in0, in1, n);
    t1 = sec();

    dump_flops("simple_omp", n, t1-t0, out_simple_omp, iter);

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
            matmul_for_a9(out_for_a9, in0, in1, n);
            t1 = sec();

            dump_flops("tune for A9", n, t1-t0, out_for_a9, iter);

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
    out_for_a9 = _aligned_malloc(n*n * sizeof(float), align);
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

    bench(0);
    bench(1);
    bench(2);

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
    _aligned_free(out_for_a9);
    _aligned_free(out_neon);
    _aligned_free(out_x86_fma);

}
