#ifndef MATMUL_BENCH_W32_APP_H
#define MATMUL_BENCH_W32_APP_H

#include "matmul-bench.h"
#include "npr/strbuf.h"
#include <windows.h>

#ifdef __cplusplus
extern "C" {
#endif

struct bench_test {
    HWND check_box;
    HWND result_box;
    double max_flops;
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
    struct MatmulBenchResult *last_result;

    HWND main_win;
    HWND start_button;
    HWND console;
    HWND iter_in;
    HWND sec_in;

    struct bench_test *tests;

    HFONT control_font;
    HFONT console_font;

    HANDLE run_thread;
    HANDLE finish_event;

    struct npr_strbuf console_text;
};

#ifdef __cplusplus
}
#endif

#endif
