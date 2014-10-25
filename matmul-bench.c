

#if 0
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

int n;

static float *in0;
static float *in1;
static float *in0_plus1line;
static float *in1_plus1line;
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
static float *out_gcc_vec4;
static float *out_gcc_vec8;

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

#endif



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
    printf("%d:%-20s: sec=%8.5f, %9.5f[GFLOPS], %8.5f[GB/s]\n",
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

    t0 = sec();
    gcc_vec4(out_gcc_vec4, in0, in1, n);
    t1 = sec();

    dump_flops("gcc_vec4", n, t1-t0, out_gcc_vec4, iter);


    t0 = sec();
    gcc_vec8(out_gcc_vec8, in0, in1, n);
    t1 = sec();

    dump_flops("gcc_vec8", n, t1-t0, out_gcc_vec8, iter);


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
            matmul_block_outer_avx_omp(out_block_outer_avx_omp, in0_plus1line, in1_plus1line, n, n+16);
            t1 = sec();

            dump_flops("avx", n, t1-t0, out_block_outer_avx_omp, iter);
#endif



#ifdef __FMA__
            t0 = sec();
            matmul_x86_fma(out_x86_fma, in0_plus1line, in1_plus1line, n, n+16);
            t1 = sec();

            dump_flops("fma", n, t1-t0, out_x86_fma, iter);
#endif


#ifdef __ARM_NEON__
            t0 = sec();
            matmul_neon(out_neon, in0_plus1line, in1_plus1line, n, n+16);
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
    in0_plus1line = _aligned_malloc((n+64)*n * sizeof(float), align);
    in1_plus1line = _aligned_malloc((n+64)*n * sizeof(float), align);

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
    out_gcc_vec4 = _aligned_malloc(n*n * sizeof(float), align);
    out_gcc_vec8 = _aligned_malloc(n*n * sizeof(float), align);

    printf("n=%d\n", n);
    printf("total nelem=%d, mat size=%f[MB]\n",
           n*n,
           (n*n)/(1024.0*1024.0)*sizeof(float));

    srand(100);

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            in0_plus1line[i*(n+16)+j] = in0[i*n+j] = drand()+1.0;
            in1_plus1line[i*(n+16)+j] = in1[i*n+j] = drand()+1.0;
        }
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
    _aligned_free(out_gcc_vec4);
    _aligned_free(out_gcc_vec8);

}
#endif