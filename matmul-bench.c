#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <omp.h>
#include <emmintrin.h>

enum bench_code {
    BENCH_SIMPLE_C,
    BENCH_CACHED_C,
#ifdef __x86_64__
    BENCH_AVX,
#endif
};

int n;

static float *in0;
static float *in1;
static float *out_simple;
static float *out_simple_outer_omp;
static float *out_simple_omp;
static float *out_block_omp;
static float *out_block_omp_unroll;

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
    printf("%-20s: sec=%f, %f[GFLOPS], %f[GB/s]\n",
           tag,
           sec,
           n*n*n*2/(sec*1024.0*1024.0*1024.0),
           (n*n*3.0*sizeof(float))/(sec*1024.0*1024.0*1024.0));

    dump_mat(ni, data);

    if (data != out_simple) {
        for (int i=0; i<n*n; i++) {
            float delta = fabs(data[i]-out_simple[i]);
            double ratio = delta/(fabs(data[i])*100.0);

            if (ratio > 1e-3) {
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

    in0 = _aligned_malloc(n*n * sizeof(float), 64);
    in1 = _aligned_malloc(n*n * sizeof(float), 64);

    out_simple = _aligned_malloc(n*n * sizeof(float), 64);
    out_simple_outer_omp = _aligned_malloc(n*n * sizeof(float), 64);
    out_simple_omp = _aligned_malloc(n*n * sizeof(float), 64);
    out_block_omp = _aligned_malloc(n*n * sizeof(float), 64);
    out_block_omp_unroll = _aligned_malloc(n*n * sizeof(float), 64);

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
}
