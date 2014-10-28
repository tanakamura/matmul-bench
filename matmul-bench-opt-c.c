#include "matmul-bench-common.h"

#ifdef HAVE_VEC_EXT
typedef float v4sf __attribute__((vector_size (16)));
typedef float v8sf __attribute__((vector_size (32)));

static void
gccvec4_run(struct MatmulBenchParam *p)
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
            v4sf vlik = {lik,lik,lik,lik};

            for (int j=0; j<n; j+=4) {
                *(v4sf*)&out[i*n+j] += vlik * *(v4sf*)&inR[k*n + j];
            }
        }
    }
}

static void
gccvec8_run(struct MatmulBenchParam *p)
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
            v8sf vlik = {lik,lik,lik,lik,
                         lik,lik,lik,lik};

            for (int j=0; j<n; j+=8) {
                *(v8sf*)&out[i*n+j] += vlik * *(v8sf*)&inR[k*n + j];
            }
        }
    }
}
#endif

static void
block_run(struct MatmulBenchParam *p)
{
    float * __restrict out = p->out;
    const float * __restrict inL = p->inL;
    const float * __restrict inR = p->inR;

    unsigned int n = p->n;

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
block_unroll_run(struct MatmulBenchParam *p)
{
    float * __restrict out = p->out;
    const float * __restrict inL = p->inL;
    const float * __restrict inR = p->inR;

    unsigned int n = p->n;

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


static const struct MatmulBenchTest block = MATMULBENCH_TEST_INITIALIZER("block", block_run, 32);
static const struct MatmulBenchTest block_unroll = MATMULBENCH_TEST_INITIALIZER("block_unroll", block_unroll_run, 16);

#ifdef HAVE_VEC_EXT
static const struct MatmulBenchTest gccvec4 = MATMULBENCH_TEST_INITIALIZER("gccvec4", gccvec4_run, 4);
static const struct MatmulBenchTest gccvec8 = MATMULBENCH_TEST_INITIALIZER("gccvec8", gccvec8_run, 8);
#endif

void
matmulbench_init_opt_c(struct MatmulBench *b, struct npr_varray *test_set)
{
    VA_PUSH(struct MatmulBenchTest, test_set, block);
    VA_PUSH(struct MatmulBenchTest, test_set, block_unroll);

#ifdef HAVE_VEC_EXT
    b->feature_bits |= MATMULBENCH_FEATURE_GCCVEC;
    VA_PUSH(struct MatmulBenchTest, test_set, gccvec4);
    VA_PUSH(struct MatmulBenchTest, test_set, gccvec8);
#endif
}
