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
simple_omp_run(struct MatmulBenchParam *p)
{
    float * __restrict out = p->out;
    const float * __restrict inL = p->inL;
    const float * __restrict inR = p->inR;

    unsigned int n = p->n;
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
}

static void
outer_run(struct MatmulBenchParam *p)
{
    float * __restrict out = p->out;
    const float * __restrict inL = p->inL;
    const float * __restrict inR = p->inR;

    unsigned int n = p->n;
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
}


static const struct MatmulBenchTest simple = MATMULBENCH_TEST_INITIALIZER("simple", simple_run, 1);
static const struct MatmulBenchTest simple_omp = MATMULBENCH_TEST_INITIALIZER("simple_omp", simple_omp_run, 1);
static const struct MatmulBenchTest outer = MATMULBENCH_TEST_INITIALIZER("outer_omp", outer_run, 1);

void
matmulbench_init_simple_c(struct MatmulBench *b, struct npr_varray *test_set)
{
    VA_PUSH(struct MatmulBenchTest, test_set, simple);
    VA_PUSH(struct MatmulBenchTest, test_set, simple_omp);
    VA_PUSH(struct MatmulBenchTest, test_set, outer);
}
