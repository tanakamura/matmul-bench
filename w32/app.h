#ifndef MATMUL_BENCH_W32_APP_H
#define MATMUL_BENCH_W32_APP_H

#include "matmul-bench.h"
#include <windows.h>

#ifdef __cplusplus
extern "C" {
#endif

struct bench_test {
    HWND check_box;
    HWND result_box;
};

enum app_state {
    APP_IDLE,
    APP_RUNNING
};

struct app {
    enum app_state st;

    HINSTANCE hInst;

    struct MatmulBench *bench;
    struct MatmulBenchConfig *config;

    HWND start_button;
    HWND console;

    struct bench_test *tests;

    HFONT control_font;
    HFONT console_font;

    HANDLE run_thread;
    HANDLE finish_event;

    size_t console_text_size;
    char *console_text;
};

#ifdef __cplusplus
}
#endif

#endif
