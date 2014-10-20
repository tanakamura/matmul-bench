#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <omp.h>

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
static float *out_simple_omp;
static float *out_block_omp;
static float *out_block_omp_unroll;



float
sec()
{
    struct timespec ts;

    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);

    return (ts.tv_sec) + (ts.tv_nsec / (1000.0*1000.0*1000.0));
}

static void
matmul_simple(float *restrict out,
              const float * restrict inL,
              const float * restrict inR,
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
matmul_simple_omp(float * restrict out,
                  const float * restrict inL,
                  const float * restrict inR,
                  int n)
{
#pragma omp parallel for
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
matmul_block_omp(float * restrict out,
                 const float* restrict inL,
                 const float* restrict inR,
                 unsigned int n)
{
    unsigned int block_size = 16;

#pragma omp parallel for
    for (int i0=0; i0<n; i0+=block_size) {
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


static void
matmul_block_omp_unroll(float * restrict out,
                        const float* restrict inL,
                        const float* restrict inR,
                        unsigned int n)
{
    unsigned int block_size = 16;

#pragma omp parallel for
    for (int i0=0; i0<n; i0+=block_size) {
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
                        const float *lp_1 = inL + i_1*n + k0;
                        const float *lp_2 = inL + i_2*n + k0;
                        const float *lp_3 = inL + i_3*n + k0;

                        const float *rp = inR + k0*n + j;

                        for (int bk=0; bk<block_size; bk++) {
                            //int k = k0+bk;

                            v_0 += *lp_0 * *rp;
                            v_1 += *lp_1 * *rp;
                            v_2 += *lp_2 * *rp;
                            v_3 += *lp_3 * *rp;

                            lp_0 += 1;
                            lp_1 += 1;
                            lp_2 += 1;
                            lp_3 += 1;

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
    n = atoi(argv[1]);

    if (argc >= 3) {
        omp_set_num_threads(atoi(argv[2]));
    }

    in0 = malloc(n*n * sizeof(float));
    in1 = malloc(n*n * sizeof(float));

    out_simple = malloc(n*n * sizeof(float));
    out_simple_omp = malloc(n*n * sizeof(float));
    out_block_omp = memalign(64, n*n * sizeof(float));
    out_block_omp_unroll = memalign(64, n*n * sizeof(float));

    printf("n=%d\n", n);
    printf("total nelem=%d, mat size=%f[MB]\n",
           n*n,
           (n*n)/(1024.0*1024.0)*sizeof(float));

    srand48(100);

    for (int i=0; i<n*n; i++) {
        in0[i] = drand48()+1.0;
        in1[i] = drand48()+1.0;
    }

    dump_mat(n, in0);
    dump_mat(n, in1);

    bench();
    bench();

    free(in0);
    free(in1);
    free(out_simple);
    free(out_simple_omp);
    free(out_block_omp);
    free(out_block_omp_unroll);
}
