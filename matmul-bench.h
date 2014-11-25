#ifndef MATMUL_BENCH_H
#define MATMUL_BENCH_H

#ifdef _WIN32
#include <windows.h>

#ifdef MATMUL_BENCH_BUILD_LIB
#define MATMUL_BENCH_EXPORT __declspec(dllexport)
#else
#define MATMUL_BENCH_EXPORT __declspec(dllimport)
#endif

#else

#define MATMUL_BENCH_EXPORT __attribute__((visibility("default")))

#endif



struct MatmulBenchTest;
struct MatmulBenchConfig;
struct MatmulBench;

struct MatmulBenchConfig {
    int iter;                   /* テスト回数 */

    int *enable;                /* test_set と対応 0:やらない 0以外:やる */
    unsigned int mat_size;     /* 計測サイズ。0で自動(下みっつのパラメータでテストする) */

    unsigned int size_min;     /* テスト開始サイズ (これより大きくて、size_stepの倍数が実際のsize_minになる) */
    unsigned int i_block_size; /* ループのブロックサイズ */
    unsigned int size_step;    /* サイズ増加 (これより大きく、かつ、全テストのsize_stepの最小公倍数が実際のstepになる) */
    double max_time_sec;        /* 処理時間がこれを超えたらやめる */
};

struct MatmulBenchTestResult {
    int num_run;
    double **sec;               /* sec[num_run][iter] */
};

struct MatmulBenchResult {
    int num_test;
    int *test_map;              /* 結果とMatmulBench::test_setの対応 */

    int num_run_max;

    unsigned int num_run;       /* 最大run数 */
    unsigned int run_size_step; /* config::size_stepを全テストの最小公倍数になるようにテストした値 */
    unsigned int run_size_min; /* config::size_minをrun_size_stepの倍になるように調整した値 */

    struct MatmulBenchTestResult *results;
};

__attribute__((aligned(64))) struct MatmulBenchParam {
    struct MatmulBench *mb;
    float * out;
    const float *inL, *inR;
    const float *inL_plus1line, *inR_plus1line;
    unsigned int n;
    unsigned int pitch_byte;
    unsigned int i_block_size;
};

typedef void (*matmul_bench_test_run_t)(struct MatmulBenchParam *p);

struct MatmulBenchTest {
    const char *name;
    matmul_bench_test_run_t run;

    unsigned int size_step;
};

typedef void (*matmul_bench_finish_callback_t)(const struct MatmulBenchTest *test,
                                               double sec,
                                               unsigned int iter,
                                               unsigned long mat_size,
                                               void *ptr);

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
    struct MatmulBenchThreadPool *threads;
};

struct MatmulBenchResult;

/* 0 でシステムのスレッド数 */
MATMUL_BENCH_EXPORT struct MatmulBench *matmul_bench_init(unsigned int num_thread);
MATMUL_BENCH_EXPORT void matmul_bench_fini(struct MatmulBench *mb);

/*
 * iter = 3, 全テスト, サイズ自動でパラメータ設定
 */
MATMUL_BENCH_EXPORT struct MatmulBenchConfig *matmul_bench_config_init(struct MatmulBench *mb);
MATMUL_BENCH_EXPORT void matmul_bench_config_fini(struct MatmulBench *mb, struct MatmulBenchConfig *c);

/* 名前が見つからなかったら-1を返す */
MATMUL_BENCH_EXPORT int matmul_bench_config_enable_test(struct MatmulBench *mb,
                                                        struct MatmulBenchConfig *config,
                                                        const char *test_name);
MATMUL_BENCH_EXPORT int matmul_bench_config_disable_test(struct MatmulBench *mb,
                                                         struct MatmulBenchConfig *config,
                                                         const char *test_name);

MATMUL_BENCH_EXPORT struct MatmulBenchResult *matmul_bench_run(struct MatmulBench *mb,
                                                               struct MatmulBenchConfig *config,
                                                               matmul_bench_finish_callback_t callback,
                                                               void *callback_data);
MATMUL_BENCH_EXPORT void matmul_bench_result_fini(struct MatmulBench *mb,
                                                  struct MatmulBenchResult *param);

#endif
