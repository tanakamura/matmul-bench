#ifndef MATMUL_BENCH_H
#define MATMUL_BENCH_H

#ifdef _WIN32
#include <windows.h>

#ifdef MATMULBENCH_BUILD_LIB
#define MATMULBENCH_EXPORT __declspec(dllexport)
#else
#define MATMULBENCH_EXPORT __declspec(dllimport)
#endif

#else
#define MATMULBENCH_EXPORT __attribute__((visivility("default")))
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

#define MATMULBENCH_FEATURE_SSE (1<<0)
#define MATMULBENCH_FEATURE_AVX (1<<1)
#define MATMULBENCH_FEATURE_FMA (1<<2)
#define MATMULBENCH_FEATURE_NEON (1<<3)

struct MatmulBench {
    int num_test;
    struct MatmulBenchTest *t;

    int feature_bits;

    const char *info;
};

struct MatmulBenchResult;

struct MatmulBench *matmul_bench_init(void);
void matmul_bench_fini(struct MatmulBench *mb);

#endif