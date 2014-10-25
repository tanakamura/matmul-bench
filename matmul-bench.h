#ifndef MATMUL_BENCH_H
#define MATMUL_BENCH_H

#ifdef _WIN32
#include <windows.h>
#endif

struct MatmulBenchTest;
struct MatmulBenchRunTestConfig;
struct MatmulBenchRunConfig;
struct MatmulBench;

struct MatmulBenchRunTestConfig {
    float *out;
    float *inL;
    float *inR;

    double *time_sec;           // [iter]
};

struct MatmulBenchRunConfig {
    unsigned int iter;
    unsigned int num_test;

    unsigned long mat_size;
    unsigned int *test_map_table;

    struct MatmulBenchRunTestConfig *test_config_list;
};

typedef void (*matmul_bench_test_run_t)(struct MatmulBench *mb,
                                        struct MatmulBenchRunConfig *run_config,
                                        struct MatmulBenchRunTestConfig *test_config,
                                        int iter);

struct MatmulBenchTest {
    char *name;
    matmul_bench_test_run_t run;

    unsigned long size_step;
};

struct MatmulBench {
    int num_test;
    struct MatmulBenchTest *t;

    const char *info;
#ifdef _WIN32
    LARGE_INTEGER freq;
#endif
};

struct MatmulBenchResult;

struct MatmulBench *matmul_bench_init(void);
void matmul_bench_fini(struct MatmulBench *mb);

#endif