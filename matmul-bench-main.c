#include <stdio.h>
#include <stdlib.h>
#include "matmul-bench.h"

enum run_type {
    RUN_ALL,
    RUN_SIZE1,
    DISPLAY_LIST,
};

static void
usage(void)
{
    puts("matmul-bench [-n size] [-l]");
    puts(" -n : matrix size");
    puts(" -l : display test info");
}

int
main(int argc, char **argv)
{
    int ai;
    enum run_type run_type = RUN_ALL;
    unsigned long size1 = 512;

    for (ai=1; ai<argc; ai++) {
        if (argv[ai][0] == '-') {
            switch (argv[ai][1]) {
            case 'l':
                run_type = DISPLAY_LIST;
                break;

            case 'n':
                run_type = RUN_SIZE1;
                if (ai == argc-1) {
                    usage();
                    exit(1);
                }

                size1 = atoi(argv[ai+1]);
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

}
