#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "matmul-bench.h"

enum run_type {
    RUN_ALL,
    RUN_SIZE1,
    DISPLAY_TEST_INFO,
};

static void
usage(void)
{
    puts("matmul-bench [-n size] [-I] [-t test1,test2,...] [-i iter] [-m <min size>] [-s <size step>] [-T <time limit>] [-o out.csv] [-O theradnum] ");
    puts(" -I : display test info");
    puts(" -n : matrix size (default : auto)");
    puts(" -t : set test list (default : all)");
    puts(" -i : num iter (default 3)");
    puts(" -m : min size (default 64)");
    puts(" -s : size step (default 64)");
    puts(" -T : time limit [float sec] (default 0.5)");
    puts(" -o : out csv");
    puts(" -O : omp thread num");
}

static void
disp_flops(const struct MatmulBenchTest *test,
           double sec,
           unsigned int iter,
           unsigned long mat_size)
{
    double n = mat_size;
    printf("(%5d-%2d):%-20s: sec=%8.5f, %9.5f[GFLOPS], %8.5f[GB/s]\n",
           (int)mat_size,
           (int)iter,
           test->name,
           sec,
           n*n*n*2/(sec*1024.0*1024.0*1024.0),
           (n*n*3.0*sizeof(float))/(sec*1024.0*1024.0*1024.0));
}

int
main(int argc, char **argv)
{
    int ai;
    enum run_type run_type = RUN_ALL;
    unsigned long size1 = 512;
    const char *test_list = NULL;
    const char *out_csv = NULL;
    int iter = 3;
    unsigned long size_min = 64;
    unsigned long size_step = 64;
    double timeout_sec = 0.5;

    for (ai=1; ai<argc; ai++) {
        if (argv[ai][0] == '-') {
            switch (argv[ai][1]) {
            case 'I':
                run_type = DISPLAY_TEST_INFO;
                break;

            case 'i':
                if (ai == argc-1) {
                    usage();
                    exit(1);
                }

                iter = atoi(argv[ai+1]);
                ai++;
                break;

            case 'T':
                if (ai == argc-1) {
                    usage();
                    exit(1);
                }

                timeout_sec = atof(argv[ai+1]);
                ai++;
                break;

            case 'n':
                run_type = RUN_SIZE1;
                if (ai == argc-1) {
                    usage();
                    exit(1);
                }

                size1 = atoi(argv[ai+1]);
                ai++;
                break;

            case 't':
                if (ai == argc-1) {
                    usage();
                    exit(1);
                }

                test_list = argv[ai+1];
                ai++;
                break;

            case 's':
                if (ai == argc-1) {
                    usage();
                    exit(1);
                }

                size_step = atoi(argv[ai+1]);
                ai++;
                break;

            case 'm':
                if (ai == argc-1) {
                    usage();
                    exit(1);
                }

                size_min = atoi(argv[ai+1]);
                ai++;
                break;

            case 'h':
                usage();
                exit(0);

            case 'o':
                if (ai == argc-1) {
                    usage();
                    exit(1);
                }
                out_csv = argv[ai+1];
                ai++;
                break;

            case 'O':
                if (ai == argc-1) {
                    usage();
                    exit(1);
                }
#ifdef _OPENMP
                omp_set_num_threads(atoi(argv[ai+1]));
#endif
                ai++;
                break;
                
                

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
    config->size_min = size_min;
    config->size_step = size_step;
    config->max_time_sec = timeout_sec;

    if (test_list) {
        int i;
        for (i=0; i<b->num_test; i++) {
            config->enable[i] = 0;
        }

        char name[512];
        size_t cur = 0;
        size_t len = strlen(test_list);

        while (cur < len) {
            for (i=0; i<511; i++, cur++) {
                if (test_list[cur] == ',' || test_list[cur] == '\0') {
                    cur++;
                    name[i] = '\0';
                    break;
                }
                name[i] = test_list[cur];
            }

            int r = matmul_bench_config_enable_test(b, config, name);
            if (r < 0) {
                fprintf(stderr, "test '%s' is not defined.\n", name);
                exit(1);
            }
        }
    }

    struct MatmulBenchResult *result = matmul_bench_run(b, config, disp_flops);

    if (out_csv) {
        FILE *fp = fopen(out_csv, "wb");
        if (fp == NULL) {
            perror(out_csv);
        } else {
            int ri, ti, ii;
            fprintf(fp, ",");
            for (ri=0; ri<result->num_run_max; ri++) {
                fprintf(fp, "%d,", (int)(result->run_size_min + ri*result->run_size_step));
            }
            fprintf(fp, "\n");

            for (ti=0; ti<result->num_test; ti++) {
                struct MatmulBenchTest *t = &b->test_set[result->test_map[ti]];
                struct MatmulBenchTestResult *tr = &result->results[ti];

                fprintf(fp, "%s,", t->name);

                for (ri=0; ri<tr->num_run; ri++) {
                    unsigned long mat_size = result->run_size_min + ri*result->run_size_step;
                    double min = DBL_MAX;

                    for (ii=0; ii<iter; ii++) {
                        if (tr->sec[ri][ii] < min) {
                            min = tr->sec[ri][ii];
                        }
                    }

                    double flops = mat_size*(double)mat_size*mat_size*2/(min*1024.0*1024.0*1024.0);

                    fprintf(fp, "%f,", flops);
                }

                fprintf(fp, "\n");
            }
        }
    }


    matmul_bench_result_fini(b, result);
    matmul_bench_config_fini(b, config);
    matmul_bench_fini(b);

    return 0;

}
