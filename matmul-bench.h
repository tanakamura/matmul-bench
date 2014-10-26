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
struct MatmulBenchConfig;
struct MatmulBench;

struct MatmulBenchConfig {
    int iter;                   /* テスト回数 */

    int *enable;                /* test_set と対応 0:やらない 0以外:やる */
    unsigned long mat_size;     /* 計測サイズ。0で自動(下みっつのパラメータでテストする) */

    unsigned long size_min;     /* テスト開始サイズ (これより大きくて、size_stepの倍数が実際のsize_minになる) */
    unsigned long size_step;    /* サイズ増加 (これより大きく、かつ、全テストのsize_stepの最小公倍数が実際のstepになる) */
    double max_time_sec;        /* 処理時間がこれを超えたらやめる */
};

struct MatmulBenchTestResult {
    int num_run;
    double **sec;               /* sec[iter][num_run] */
};

struct MatmulBenchResult {
    int num_test;
    int *test_map;              /* 結果とMatmulBench::test_setの対応 */

    unsigned int num_run;       /* 最大run数 */
    unsigned long run_size_step; /* config::size_stepを全テストの最小公倍数になるようにテストした値 */
    unsigned long run_size_min; /* config::size_minをrun_size_stepの倍になるように調整した値 */

    struct MatmulBenchTestResult *results;
};

typedef void (*matmul_bench_test_run_t)(float * __restrict out,
                                        const float * __restrict inL,
                                        const float * __restrict inR,

                                        const float * __restrict inL_plus1line,
                                        const float * __restrict inR_plus1line,
                                        
                                        unsigned int n,
                                        unsigned int pitch_byte);

struct MatmulBenchTest {
    const char *name;
    matmul_bench_test_run_t run;

    unsigned long size_step;
};

typedef void (*matmul_bench_finish_callback_t)(const struct MatmulBenchTest *test,
                                               double sec,
                                               unsigned int iter,
                                               unsigned long mat_size);

#define MATMULBENCH_FEATURE_SSE (1<<0)
#define MATMULBENCH_FEATURE_AVX (1<<1)
#define MATMULBENCH_FEATURE_FMA (1<<2)
#define MATMULBENCH_FEATURE_NEON (1<<3)
#define MATMULBENCH_FEATURE_VFPV4 (1<<4)
#define MATMULBENCH_FEATURE_GCCVEC (1<<5)

struct MatmulBench {
    int num_test;
    struct MatmulBenchTest *test_set;
    int feature_bits;
};

struct MatmulBenchResult;

struct MatmulBench *matmul_bench_init(void);
void matmul_bench_fini(struct MatmulBench *mb);

/*
 * iter = 3, 全テスト, サイズ自動でパラメータ設定
 */
struct MatmulBenchConfig *matmul_bench_config_init(struct MatmulBench *mb);
void matmul_bench_config_fini(struct MatmulBench *mb, struct MatmulBenchConfig *c);

/* 名前が見つからなかったら-1を返す */
int matmul_bench_config_enable_test(struct MatmulBench *mb,
                                    struct MatmulBenchConfig *config,
                                    const char *test_name);
int matmul_bench_config_disable_test(struct MatmulBench *mb,
                                     struct MatmulBenchConfig *config,
                                     const char *test_name);

struct MatmulBenchResult *matmul_bench_run(struct MatmulBench *mb,
                                           struct MatmulBenchConfig *config,
                                           matmul_bench_finish_callback_t callback);
void matmul_bench_result_fini(struct MatmulBench *mb,
                              struct MatmulBenchResult *param);

#endif