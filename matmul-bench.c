#define _WIN32_WINNT 0x0600
#define _CRT_RAND_S

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "matmul-bench.h"
#include "matmul-bench-common.h"
#include "npr/varray.h"

#ifdef ARCH_X86
#include <cpuid.h>
#define rmb() __asm__ __volatile__ ("":::"memory")
#define wmb() __asm__ __volatile__ ("":::"memory")
#endif

#ifdef __arm__
#define rmb() __asm__ __volatile__ ("dsb":::"memory")
#define wmb() __asm__ __volatile__ ("dsb st":::"memory")
#endif

#ifdef linux
#include <sys/eventfd.h>
#endif

#ifdef _WIN32
#include <windows.h>
#include <process.h>

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

static double
drand(void)
{
    unsigned int v;
    rand_s(&v);

    return v / (double)UINT_MAX;
}

#define srand(v)


static void
notify_event(HANDLE ev)
{
    wmb();
    SetEvent(ev);
}

static void
wait_event(HANDLE ev)
{
    WaitForSingleObject(ev, INFINITE);
    rmb();
}


#else
#include <unistd.h>


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

static void
notify_event(int fd)
{
    uint64_t ev_val = 1;
    wmb();
    ssize_t s = write(fd, &ev_val, sizeof(ev_val));
    if (s != sizeof(uint64_t)) {
        perror("write");    /* ?? */
    }
}

static void
wait_event(int fd)
{
    uint64_t ev_val;

    ssize_t s = read(fd, &ev_val, sizeof(ev_val));
    if (s != sizeof(uint64_t)) {
        perror("read");    /* ?? */
    }

    rmb();
}


#endif

static void *
thread_func(void *ap)
{
    struct MatmulBenchThreadArg *a = (struct MatmulBenchThreadArg*)ap;
    struct MatmulBench *bench = a->b;
    struct MatmulBenchThreadPool *pool = bench->threads;
    unsigned int *current_i = pool->current_i;

    while (1) {
        wait_event(a->from_master_ev);

        if (pool->fini) {
            break;
        }

        if (pool->fini) {
            break;
        }

        unsigned int max_i = pool->max_i;
        unsigned int i_block_size = pool->i_block_size;
        matmul_bench_thread_func_t func = pool->func;

        struct MatmulBenchParam *param = pool->param;;

        while (1) {
            unsigned long i;

            i = __sync_fetch_and_add(current_i, i_block_size);

            if (i >= max_i) {
                break;
            }

            unsigned long end = i+i_block_size;
            end = (end<max_i)?end:max_i;

            func(param, i, end);
        }

        a->fini = 1;

        notify_event(pool->to_master_ev);
    }

    return NULL;
}

#ifdef _WIN32
static unsigned __stdcall
thread_func_w32(void *ap)
{
    thread_func(ap);
    return 0;
}
#endif

struct MatmulBench *
matmul_bench_init(unsigned int num_thread)
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

#ifdef _WIN32
    if (num_thread == 0) {
        SYSTEM_INFO si;

        GetSystemInfo(&si);

        num_thread = si.dwNumberOfProcessors;
    }
#else
    if (num_thread == 0) {
        num_thread = sysconf(_SC_NPROCESSORS_ONLN);
    }
#endif

    struct MatmulBenchThreadPool *pl = malloc(sizeof(*pl));
    ret->threads = pl;
    pl->num_thread = num_thread;
    pl->current_i = _aligned_malloc(64, 64);
    pl->fini = 0;
    pl->args = (struct MatmulBenchThreadArg*)_aligned_malloc(sizeof(struct MatmulBenchThreadArg) * num_thread, 64);

    for (int ti=0; ti<num_thread; ti++) {
        struct MatmulBenchThreadArg *a = &pl->args[ti];
        a->b = ret;

#ifdef _WIN32
        unsigned int threadID;

        a->from_master_ev = CreateEvent(NULL, FALSE, FALSE, NULL);
        a->thread = (HANDLE)_beginthreadex(NULL, 0, thread_func_w32, a, 0, &threadID);
#else
        a->from_master_ev = eventfd(0,0);

        pthread_create(&a->thread, NULL,
                       thread_func, a);
#endif
    }


#ifdef _WIN32
    pl->to_master_ev = CreateEvent(NULL, FALSE, FALSE, NULL);
#else
    pl->to_master_ev = eventfd(0,0);
#endif

    ret->num_test = test_set.nelem;
    ret->test_set = npr_varray_malloc_close(&test_set);

    return ret;
}

void
matmul_bench_fini(struct MatmulBench *mb)
{
    struct MatmulBenchThreadPool *pl = mb->threads;
    int num_thread = pl->num_thread;
    uint64_t ev_val = 1;

    pl->fini = 1;

    for (int ti=0; ti<num_thread; ti++) {
        notify_event(pl->args[ti].from_master_ev);
    }

    for (int ti=0; ti<num_thread; ti++) {
#ifdef _WIN32
        WaitForSingleObject(pl->args[ti].thread, INFINITE);
        CloseHandle(pl->args[ti].thread);
        CloseHandle(pl->args[ti].from_master_ev);
#else
        pthread_join(pl->args[ti].thread, NULL);
        close(pl->args[ti].from_master_ev);
#endif
    }

    free(pl->args);
    _aligned_free(pl->current_i);

    free(mb->test_set);
    free(mb);
}
void
matmul_bench_thread_call(struct MatmulBenchParam *param,
                         unsigned int i_block_size,
                         unsigned int max_i,
                         matmul_bench_thread_func_t func)
{
    struct MatmulBench *mb = param->mb;
    struct MatmulBenchThreadPool *pl = mb->threads;
    int num_thread = pl->num_thread;

    *(pl->current_i) = 0;
    pl->fini = 0;
    pl->i_block_size = i_block_size;
    pl->max_i = max_i;
    pl->func = func;
    pl->param = param;

    for (int ti=0; ti<num_thread; ti++) {
        pl->args[ti].fini = 0;

        notify_event(pl->args[ti].from_master_ev);
    }

    while (1) {
        wait_event(pl->to_master_ev);

        int all_fini = 1;

        for (int ti=0; ti<num_thread; ti++) {
            all_fini &= pl->args[ti].fini;
        }

        if (all_fini) {
            break;
        }
    }
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
    c->i_block_size = 1;

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

struct RunTestState {
    struct npr_varray time_list;
    double *sec_buffer;
};

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
    r->test_map = test_map;
    r->num_test = num_enable;

    int *timeout = malloc(sizeof(int) * num_enable);
    int *timeout_cur_size = malloc(sizeof(int) * num_enable);
    struct RunTestState *run_st;
    run_st = malloc(sizeof(*run_st) * num_enable);


    int ti = 0;
    unsigned long test_size_step_lcm = 1;

    for (i=0; i<b->num_test; i++) {
        if (c->enable[i]) {
            timeout[ti] = 0;
            test_map[ti] = i;

            npr_varray_init(&run_st[ti].time_list, 16, sizeof(double*));

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
    r->run_size_min = run_size_min;
    r->run_size_step = run_size_step;

    struct MatmulBenchParam *run_param = _aligned_malloc(sizeof(*run_param), 64);
    run_param->mb = b;
    run_param->i_block_size = c->i_block_size;

    while (1) {
        unsigned long n = cur_size;
        int align = 64;

        float *in0_plus1line = _aligned_malloc((n+64)*n * sizeof(float), align);
        float *in1_plus1line = _aligned_malloc((n+64)*n * sizeof(float), align);
        float *in0 = _aligned_malloc(n*n * sizeof(float), align);
        float *in1 = _aligned_malloc(n*n * sizeof(float), align);

        for (int i=0; i<num_enable; i++) {
            out_set[i] = _aligned_malloc(n*n * sizeof(float), align);
            timeout_cur_size[i] = 0;
            if (!timeout[i]) {
                run_st[i].sec_buffer = malloc(c->iter * sizeof(double));
            }
        }

        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                in0_plus1line[i*(n+16)+j] = in0[i*n+j] = drand()+1.0;
                in1_plus1line[i*(n+16)+j] = in1[i*n+j] = drand()+1.0;
            }
        }

        run_param->inL = in0;
        run_param->inR = in1;
        run_param->inL_plus1line = in0_plus1line;
        run_param->inR_plus1line = in1_plus1line;
        run_param->n = n;
        run_param->pitch_byte = (n+16)*4;

        for (int iter_i=0; iter_i<c->iter; iter_i++) {
            int first_run = -1;

            for (int ti=0; ti<num_enable; ti++) {
                if (!timeout[ti]) {
                    float *out = out_set[ti];
                    struct MatmulBenchTest *t = &b->test_set[test_map[ti]];

                    double tb = matmul_bench_sec();

                    run_param->out = out;

                    t->run(run_param);

                    double te = matmul_bench_sec();
                    double dt = te-tb;

                    callback(t, dt, iter_i, n);

                    if (dt > c->max_time_sec) {
                        timeout_cur_size[ti] = 1;
                    }

                    run_st[ti].sec_buffer[iter_i] = dt;

                    if (first_run == -1) {
                        first_run = ti;
                    } else {
                        float *out0 = out_set[first_run];

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
        }

        int all_timeout = 1;
        for (int i=0; i<num_enable; i++) {
            if (!timeout[i]) {
                VA_PUSH(double *, &run_st[i].time_list, run_st[i].sec_buffer);
                run_st[i].sec_buffer = NULL;
            }

            _aligned_free(out_set[i]);
            timeout[i] |= timeout_cur_size[i];
            all_timeout &= timeout[i];
        }

        _aligned_free(in0);
        _aligned_free(in1);
        _aligned_free(in0_plus1line);
        _aligned_free(in1_plus1line);

        if (all_timeout || c->mat_size) {
            break;
        }

        cur_size += run_size_step;
    }

    _aligned_free(run_param);

    struct MatmulBenchTestResult *test_results = malloc(sizeof(*test_results) * num_enable);

    int num_run_max = 0;

    for (int i=0; i<num_enable; i++) {
        int num_run = run_st[i].time_list.nelem;
        test_results[i].num_run = num_run;
        test_results[i].sec = npr_varray_malloc_close(&run_st[i].time_list);

        if (num_run > num_run_max) {
            num_run_max = num_run;
        }
    }

    r->num_run_max = num_run_max;
    r->results = test_results;

    free(run_st);
    free(out_set);
    free(timeout);
    free(timeout_cur_size);

    return r;
}

void
matmul_bench_result_fini(struct MatmulBench *b,
                         struct MatmulBenchResult *r)
{
    for (int ti=0; ti<r->num_test; ti++) {
        struct MatmulBenchTestResult *tr = &r->results[ti];
        for (int ri=0; ri<tr->num_run; ri++) {
            free(tr->sec[ri]);
        }

        free(tr->sec);
    }

    free(r->results);
    free(r->test_map);
    free(r);
}

