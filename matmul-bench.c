#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <omp.h>

#ifdef __SSE__
#include <x86intrin.h>
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
static float *out_block_outer_neon_omp;

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

    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);

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
    unsigned int block_size = 16;
    int i0;

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

                    for (int bj=0; bj<block_size; bj++) {
                        float v = 0;
                        int j = j0+bj;

                        const float *lp = inL + i*n + k0;
                        const float *rp = inR + k0*n + j;

                        for (int bk=0; bk<block_size; bk++) {
                            //int k = k0+bk;
                            v += *lp * *rp;
                            lp += 1;
                            rp += n;
                        }

                        out[i*n+j] += v;
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
    for (i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            out[i*n+j] = 0;
        }
    }

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
                    __m128 *outp4 = (__m128*)outp;

                    __m128 vout0 = outp4[0];
                    __m128 vout1 = outp4[1];
                    __m128 vout2 = outp4[2];
                    __m128 vout3 = outp4[3];

                    __m128 vout4 = outp4[4];
                    __m128 vout5 = outp4[5];
                    __m128 vout6 = outp4[6];
                    __m128 vout7 = outp4[7];
 
                    _mm_prefetch(outp + n ,_MM_HINT_T0);

                    for (int bk=0; bk<block_size; bk++) {
                        int k = k0+bk;

                        const float *inRp = &inR[k*n+j0];

                        _mm_prefetch(inRp + n ,_MM_HINT_T0);

                        float lik = inL[i*n+k];
                        __m128 lik4 = _mm_set1_ps(lik);

#define OUTER_SSE_J(J)                                               \
                        int j_##J = (J*4);                            \
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

                    outp4[0] = vout0;
                    outp4[1] = vout1;
                    outp4[2] = vout2;
                    outp4[3] = vout3;

                    outp4[4] = vout4;
                    outp4[5] = vout5;
                    outp4[6] = vout6;
                    outp4[7] = vout7;
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
    for (i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            out[i*n+j] = 0;
        }
    }

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
                    __m256 *outp8 = (__m256*)outp;

                    __m256 vout0 = outp8[0];
                    __m256 vout1 = outp8[1];
                    __m256 vout2 = outp8[2];
                    __m256 vout3 = outp8[3];

                    _mm_prefetch(outp + n ,_MM_HINT_T0);

                    for (int bk=0; bk<block_size; bk++) {
                        int k = k0+bk;

                        const float *inRp = &inR[k*n+j0];

                        _mm_prefetch(inRp + n ,_MM_HINT_T0);

                        float lik = inL[i*n+k];
                        __m256 lik8 = _mm256_set1_ps(lik);

#define OUTER_AVX_J(J)                                               \
                        int j_##J = (J*8);                            \
                        __m256 vr##J = _mm256_load_ps(&inRp[j_##J]); \
                        vout##J = _mm256_add_ps(vout##J, _mm256_mul_ps(lik8, vr##J)); \

                        OUTER_AVX_J(0);
                        OUTER_AVX_J(1);
                        OUTER_AVX_J(2);
                        OUTER_AVX_J(3);
                    }

                    outp8[0] = vout0;
                    outp8[1] = vout1;
                    outp8[2] = vout2;
                    outp8[3] = vout3;
                }
            }
        }
    }
}
#endif


#ifdef __ARM_NEON__
static void
matmul_block_outer_neon_omp(float * __restrict out,
                            const float* __restrict inL,
                            const float* __restrict inR,
                            unsigned int n)
{
    unsigned int block_size = 16;
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

                    for (int bk=0; bk<block_size; bk++) {
                        int k = k0+bk;

                        const float32x2_t *inRp = (float32x2_t*)&inR[k*n+j0];

                        float lik = inL[i*n+k];
                        float32x2_t lik2 = vdup_n_f32(lik);

                        __builtin_prefetch(inRp + n*3);

#define OUTER_NEON_J0(J)                                               \
                        float32x2_t vr##J = *(inRp++);                 \
                        vout##J = vmla_f32(vout##J, lik2, vr##J);      \

                        OUTER_NEON_J0(0);
                        OUTER_NEON_J0(1);
                        OUTER_NEON_J0(2);
                        OUTER_NEON_J0(3);
                        OUTER_NEON_J0(4);
                        OUTER_NEON_J0(5);
                        OUTER_NEON_J0(6);
                        OUTER_NEON_J0(7);

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
           float *data)
{
    double n = ni;
    printf("%-20s: sec=%8.5f, %8.5f[GFLOPS], %8.5f[GB/s]\n",
           tag,
           sec,
           n*n*n*2/(sec*1024.0*1024.0*1024.0),
           (n*n*3.0*sizeof(float))/(sec*1024.0*1024.0*1024.0));

    dump_mat(ni, data);

    if (data != out_simple) {
        for (int i=0; i<n*n; i++) {
            float delta = fabs(data[i]-out_simple[i]);
            double ratio = delta/(fabs(data[i])*100.0);

            if (ratio > 1e-5) {
                printf("error delta=%e(%e[%%]), simple=%e, opt=%e\n", delta, ratio, data[i], out_simple[i]);
                exit(1);
            }
        }
    }

}

static void
bench(void)
{
    double t0, t1;

    t0 = sec();
    matmul_simple(out_simple, in0, in1, n);
    t1 = sec();

    dump_flops("simple", n, t1-t0, out_simple);


    t0 = sec();
    matmul_simple_outer_omp(out_simple_outer_omp, in0, in1, n);
    t1 = sec();

    dump_flops("simple_outer_omp", n, t1-t0, out_simple_outer_omp);

    t0 = sec();
    matmul_simple_omp(out_simple_omp, in0, in1, n);
    t1 = sec();

    dump_flops("simple_omp", n, t1-t0, out_simple_omp);

    if (n > 16) {
        t0 = sec();
        matmul_block_omp(out_block_omp, in0, in1, n);
        t1 = sec();

        dump_flops("block_omp", n, t1-t0, out_block_omp);

        if (n > 64) {
            t0 = sec();
            matmul_block_omp_unroll(out_block_omp_unroll, in0, in1, n);
            t1 = sec();

            dump_flops("block_omp_unroll", n, t1-t0, out_block_omp_unroll);


#ifdef __SSE__
            t0 = sec();
            matmul_block_outer_sse_omp(out_block_outer_sse_omp, in0, in1, n);
            t1 = sec();

            dump_flops("sse", n, t1-t0, out_block_omp_unroll);
#endif

#ifdef __AVX__
            t0 = sec();
            matmul_block_outer_avx_omp(out_block_outer_avx_omp, in0, in1, n);
            t1 = sec();

            dump_flops("avx", n, t1-t0, out_block_omp_unroll);
#endif


#ifdef __ARM_NEON__
            t0 = sec();
            matmul_block_outer_neon_omp(out_block_outer_neon_omp, in0, in1, n);
            t1 = sec();

            dump_flops("neon", n, t1-t0, out_block_outer_neon_omp);
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

    if (argc >= 3) {
        omp_set_num_threads(atoi(argv[2]));
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

    bench();
    bench();

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

}
