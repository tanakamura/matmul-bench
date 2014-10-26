#include <stdio.h>
#include <stdlib.h>
#include "matmul-bench.h"

enum run_type {
    RUN_ALL,
    RUN_SIZE1,
    DISPLAY_TEST_INFO,
};

static void
usage(void)
{
    puts("matmul-bench [-n size] [-i] [-t test1,test2,...] [-i iter]");
    puts(" -i : display test info");
    puts(" -n : matrix size (default auto)");
    puts(" -t : set test list (default all)");
    puts(" -i : num iter (default 3)");
}

int
main(int argc, char **argv)
{
    int ai;
    enum run_type run_type = RUN_ALL;
    unsigned long size1 = 512;
    const char *test_list;
    int iter = 3;

    for (ai=1; ai<argc; ai++) {
        if (argv[ai][0] == '-') {
            switch (argv[ai][1]) {
            case 'i':
                run_type = DISPLAY_TEST_INFO;
                break;

            case 'n':
                run_type = RUN_SIZE1;
                if (ai == argc-1) {
                    usage();
                    exit(1);
                }

                size1 = atoi(argv[ai+1]);
                break;

            case 't':
                if (ai == argc-1) {
                    usage();
                    exit(1);
                }

                test_list = argv[ai+1];
                break;

            case 'h':
                usage();
                exit(0);

            default:
                usage();
                exit(1);
            }
        }
    }

    if (argc != ai) {
        usage();
        exit(0);
    }

    struct MatmulBench *b = matmul_bench_init();
    if (run_type == DISPLAY_TEST_INFO) {
        int i;
        printf("<cpu feature =");

        if (b->feature_bits & MATMULBENCH_FEATURE_SSE) {
            printf(" sse");
        }
        if (b->feature_bits & MATMULBENCH_FEATURE_AVX) {
            printf(" avx");
        }
        if (b->feature_bits & MATMULBENCH_FEATURE_FMA) {
            printf(" fma");
        }
        if (b->feature_bits & MATMULBENCH_FEATURE_NEON) {
            printf(" neon");
        }
        if (b->feature_bits & MATMULBENCH_FEATURE_VFPV4) {
            printf(" vfpv4");
        }
        if (b->feature_bits & MATMULBENCH_FEATURE_GCCVEC) {
            printf(" gccvec");
        }
        printf(">\n");

        for (i=0; i<b->num_test; i++) {
            printf("%15s: size_step=%ld\n", b->test_set[i].name, b->test_set[i].size_step);
        }

        matmul_bench_fini(b);
        exit(0);
    }

    struct MatmulBenchConfig *config = matmul_bench_config_init(b);

    if (run_type == RUN_SIZE1) {
        config->mat_size = size1;
    }

    config->iter = iter;

    matmul_bench_config_fini(b, config);
    matmul_bench_fini(b);

    return 0;

}
