#define _WIN32_WINNT 0x0600
#define _CRT_RAND_S

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "matmul-bench.h"
#include "matmul-bench-common.h"
#include "npr/varray.h"

#ifdef ARCH_X86
#include <cpuid.h>
#endif


#ifdef _WIN32
#include <windows.h>

static LARGE_INTEGER freq;
static INIT_ONCE g_InitOnce = INIT_ONCE_STATIC_INIT;

static BOOL CALLBACK
init1(PINIT_ONCE InitOnce,
      PVOID Parameter,
      LPVOID *lpContext)
{
    QueryPerformanceFrequency(&freq);
    return TRUE;
}

double
matmul_bench_sec(void)
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


double
matmul_bench_sec(void)
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

struct MatmulBench *
matmul_bench_init(void)
{
    struct MatmulBench *ret = malloc(sizeof(struct MatmulBench));

#ifdef _WIN32
    InitOnceExecuteOnce(&g_InitOnce, init1,
                        NULL, NULL);
#endif


    ret->feature_bits = 0;

    struct npr_varray test_set;

    npr_varray_init(&test_set, 16, sizeof(struct MatmulBenchTest));

    matmulbench_init_simple_c(ret, &test_set);
    matmulbench_init_opt_c(ret, &test_set);

#ifdef ARCH_X86
    unsigned int eax=0, ebx=0, ecx=0, edx=0;
    __get_cpuid(1, &eax, &ebx, &ecx, &edx);

    if (edx & (1<<25)) {
        matmulbench_init_sse(ret, &test_set);
        ret->feature_bits |= MATMULBENCH_FEATURE_SSE;
    }

    if ((ecx & 0x18000000) == 0x18000000) {
        matmulbench_init_avx(ret, &test_set);
        ret->feature_bits |= MATMULBENCH_FEATURE_AVX;
    }

    if ((ecx & (1<<12))) {
        matmulbench_init_fma(ret, &test_set);
        ret->feature_bits |= MATMULBENCH_FEATURE_FMA;
    }
#endif


#ifdef __arm__

    {
        FILE *fp = fopen("/proc/cpuinfo", "rb");

        if (fp) {
            while (1) {
                char line[1024];
                if (fgets(line, 1024, fp) == NULL) {
                    break;
                }
                if (strncmp(line, "Features", 8) == 0) {
                    int cur = 0;
                    size_t len = strlen(line);
                    for (cur=0; cur<len; cur++) {
                        if (line[cur]==':') {
                            cur++;
                            break;
                        }
                    }


                    while (cur < len) {
                        if (strncmp(line+cur, "vfpv4", 5) == 0) {
                            ret->feature_bits |= MATMULBENCH_FEATURE_VFPV4;
                        }
                        if (strncmp(line+cur, "neon", 4) == 0) {
                            ret->feature_bits |= MATMULBENCH_FEATURE_NEON;
                        }
                        while (cur<len) {
                            if (line[cur] == ' ') {
                                cur++;
                                break;
                            }
                            cur++;
                        }
                    }
                }
            }
            fclose(fp);
        }
    }

    if (ret->feature_bits & MATMULBENCH_FEATURE_NEON) {
        matmulbench_init_neon(ret, &test_set);
    }
    if (ret->feature_bits & MATMULBENCH_FEATURE_VFPV4) {
        matmulbench_init_vfpv4(ret, &test_set);
    }

#endif

    ret->num_test = test_set.nelem;
    ret->test_set = npr_varray_malloc_close(&test_set);

    return ret;
}

void
matmul_bench_fini(struct MatmulBench *mb)
{
    free(mb->test_set);
    free(mb);
}


struct MatmulBenchConfig *
matmul_bench_config_init(struct MatmulBench *mb)
{
    int i;
    struct MatmulBenchConfig *c = malloc(sizeof(*c));

    c->iter = 3;
    c->enable = malloc(sizeof(int) * mb->num_test);

    for (i=0; i<mb->num_test; i++) {
        c->enable[i] = 1;
    }

    c->mat_size = 0;
    c->size_min = 1;
    c->size_step = 1;
    c->max_time_sec = 2.0;

    return c;
}

void
matmul_bench_config_fini(struct MatmulBench *mb,
                         struct MatmulBenchConfig *mbc)
{
    free(mbc->enable);
    free(mbc);
}

int
matmul_bench_config_enable_test(struct MatmulBench *mb,
                                struct MatmulBenchConfig *config,
                                const char *test_name)
{
    int i;
    for (i=0; i<mb->num_test; i++) {
        if (strcmp(mb->test_set[i].name, test_name) == 0) {
            config->enable[i] = 1;
            return 0;
        }
    }

    return -1;
}

static unsigned long
gcd(unsigned long m, unsigned long n)
{
    if (m < n) {
        unsigned long t;
        t = n;
        n = m;
        m = t;
    }

    while (1) {
        if (n == 0) {
            break;
        }

        unsigned long mod = m%n;

        m = n;
        n = mod;
    }

    return m;
}

static unsigned long
lcm(unsigned long m, unsigned long n)
{
    unsigned long nm = m*n;

    return nm / gcd(m,n);
}


struct MatmulBenchResult *
matmul_bench_run(struct MatmulBench *b,
                 struct MatmulBenchConfig *c,
                 matmul_bench_finish_callback_t callback)
{
    struct MatmulBenchResult *r = malloc(sizeof(*r));
    int num_enable = 0;
    int i;
    for (i=0; i<b->num_test; i++) {
        if (c->enable[i]) {
            num_enable++;
        }
    }

    int *test_map = malloc(sizeof(int) * num_enable);

    int *timeout = malloc(sizeof(int) * num_enable);
    int *timeout_size = malloc(sizeof(int) * num_enable);

    int ti = 0;
    unsigned long test_size_step_lcm = 1;

    for (i=0; i<b->num_test; i++) {
        if (c->enable[i]) {
            timeout[ti] = 0;
            test_map[ti] = i;
            ti++;

            test_size_step_lcm = lcm(test_size_step_lcm, b->test_set[i].size_step);
        }
    }

    unsigned long run_size_step;
    unsigned long run_size_min;

    run_size_step = c->size_step;

    if (run_size_step < test_size_step_lcm) {
        run_size_step = test_size_step_lcm;
    } else {
        run_size_step = CEIL_DIV(run_size_step,test_size_step_lcm) * test_size_step_lcm;
    }     

    if (c->mat_size) {
        run_size_min = CEIL_DIV(c->mat_size, run_size_step) * run_size_step;
    } else {
        run_size_min = CEIL_DIV(c->size_min, run_size_step) * run_size_step;
    }

    unsigned long cur_size = run_size_min;
    float **out_set = malloc(sizeof(float*) * num_enable);

    while (1) {
        unsigned long n = cur_size;
        int align = 64;
        int all_timeout = 1;

        float *in0_plus1line = _aligned_malloc((n+64)*n * sizeof(float), align);
        float *in1_plus1line = _aligned_malloc((n+64)*n * sizeof(float), align);
        float *in0 = _aligned_malloc(n*n * sizeof(float), align);
        float *in1 = _aligned_malloc(n*n * sizeof(float), align);

        for (int i=0; i<num_enable; i++) {
            out_set[i] = _aligned_malloc(n*n * sizeof(float), align);
        }

        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                in0_plus1line[i*(n+16)+j] = in0[i*n+j] = drand()+1.0;
                in1_plus1line[i*(n+16)+j] = in1[i*n+j] = drand()+1.0;
            }
        }

        for (int iter_i=0; iter_i<c->iter; iter_i++) {
            for (int ti=0; ti<num_enable; ti++) {
                float *out = out_set[ti];
                struct MatmulBenchTest *t = &b->test_set[test_map[ti]];

                double tb = matmul_bench_sec();

                t->run(out,
                       in0, in1,
                       in0_plus1line, in1_plus1line,
                       n,
                       (n+16)*4);

                double te = matmul_bench_sec();

                callback(t, te-tb, iter_i, n);

                if (ti != 0) {
                    float *out0 = out_set[0];

                    for (int i=0; i<n*n; i++) {
                        float delta = fabs(out[i]-out0[i]);
                        double ratio = (delta/fabs(out0[i]))*100;

                        if (ratio > 1e-3) {
                            printf("error delta=%e(%e[%%]), simple=%e, opt=%e\n", delta, ratio, out[i], out0[i]);
                            exit(1);
                        }
                    }
                }
            }
        }

        for (int i=0; i<num_enable; i++) {
            free(out_set[i]);
        }

        free(in0);
        free(in1);
        free(in0_plus1line);
        free(in1_plus1line);

        if (all_timeout || c->mat_size) {
            break;
        }

        cur_size += run_size_step;
    }

    free(out_set);
    free(timeout_size);
    free(timeout);

    r->test_map = test_map;

    return r;
}

void
matmul_bench_result_fini(struct MatmulBench *b,
                         struct MatmulBenchResult *r)
{
    free(r->test_map);
    free(r);
}       


#if 0
#include <stdio.h>
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
