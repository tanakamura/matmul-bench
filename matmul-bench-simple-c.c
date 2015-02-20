#include "matmul-bench-common.h"

static void
simple_run(struct MatmulBenchParam *p)
{
    float * __restrict out = p->out;
    const float * __restrict inL = p->inL;
    const float * __restrict inR = p->inR;

    unsigned int n = p->n;

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            float v = 0;
            for (int k=0; k<n; k++) {
                v += inL[i*n+k] * inR[k*n+j];
            }

            out[i*n+j] = v;
        }
    }
}


static void
thread_func(struct MatmulBenchParam *p,
            unsigned long i_start,
            unsigned long i_end,
            unsigned int thread_id)
{
    float * __restrict out = p->out;
    const float * __restrict inL = p->inL;
    const float * __restrict inR = p->inR;
    unsigned long n = p->n;

    for (unsigned long i=i_start; i<i_end; i++) {
        for (unsigned long j=0; j<n; j++) {
            float v = 0;
            for (unsigned long k=0; k<n; k++) {
                v += inL[i*n+k] * inR[k*n+j];
            }

            out[i*n+j] = v;
        }
    }
}

static void
simple_thread_run(struct MatmulBenchParam *p)
{
    matmul_bench_thread_call(p, p->i_block_size, p->n, thread_func);
}


static void
outer_func(struct MatmulBenchParam *p,
           unsigned long i_start,
           unsigned long i_end,
           unsigned int thread_id)
{
    const float * __restrict inL = p->inL;
    const float * __restrict inR = p->inR;
    float * __restrict out = p->out;
    unsigned long n = p->n;

    for (unsigned long i=i_start; i<i_end; i++) {
        for (int k=0; k<n; k++) {
            float lik = inL[i*n+k];

            for (int j=0; j<n; j++) {
                out[i*n+j] += lik * inR[k*n + j];
            }
        }
    }
}


static void
outer_run(struct MatmulBenchParam *p)
{
    float * __restrict out = p->out;

    unsigned int n = p->n;
    int i;

    for (i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            out[i*n+j] = 0;
        }
    }

    matmul_bench_thread_call(p, p->i_block_size, p->n, outer_func);
}


static const struct MatmulBenchTest simple = MATMULBENCH_TEST_INITIALIZER("simple", simple_run, 1);
static const struct MatmulBenchTest simple_thread = MATMULBENCH_TEST_INITIALIZER("simple_thread", simple_thread_run, 1);
static const struct MatmulBenchTest outer = MATMULBENCH_TEST_INITIALIZER("outer_thread", outer_run, 1);

void
matmulbench_init_simple_c(struct MatmulBench *b, struct npr_varray *test_set)
{
    VA_PUSH(struct MatmulBenchTest, test_set, simple);
    VA_PUSH(struct MatmulBenchTest, test_set, simple_thread);
    VA_PUSH(struct MatmulBenchTest, test_set, outer);
}
