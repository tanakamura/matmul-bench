#include <windows.h>
#include <windowsx.h>
#include <process.h>
#include <stdio.h>
#include "resource.h"
#include "app.h"

#define GET_APP()  struct app *app = (struct app*)GetWindowLongPtr(hWnd, GWLP_USERDATA)

static void
on_destroy(HWND hWnd)
{
    PostQuitMessage(0);
}

static int
on_create(HWND hWnd, LPCREATESTRUCT cs)
{
    struct app *app = (struct app*)cs->lpCreateParams;
    app->main_win = hWnd;

    SetWindowLongPtr(hWnd, GWLP_USERDATA, (LONG_PTR)app);

    npr_strbuf_init(&app->console_text);

    app->bench = matmul_bench_init(0);
    app->config = matmul_bench_config_init(app->bench);
    app->last_result = NULL;

    app->control_font = CreateFont(18, 0, 0, 0, FW_DONTCARE,
                                   FALSE, FALSE, FALSE,
                                   DEFAULT_CHARSET,
                                   OUT_DEFAULT_PRECIS,
                                   CLIP_DEFAULT_PRECIS,
                                   DEFAULT_QUALITY,
                                   DEFAULT_PITCH,
                                   "ƒƒCƒŠƒI");

    app->console_font = CreateFont(18, 0, 0, 0, FW_DONTCARE,
                                   FALSE, FALSE, FALSE,
                                   DEFAULT_CHARSET,
                                   OUT_DEFAULT_PRECIS,
                                   CLIP_DEFAULT_PRECIS,
                                   DEFAULT_QUALITY,
                                   FIXED_PITCH,
                                   "Courier New");

    app->start_button = CreateWindow("BUTTON",
                                     "Start",
                                     WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
                                     20,
                                     20,
                                     150,
                                     30,
                                     hWnd,
                                     (HMENU)IDC_START,
                                     app->hInst,
                                     NULL);

    app->console = CreateWindowEx(WS_EX_CLIENTEDGE,
                                  "EDIT",
                                  "Welcome!\n",
                                  WS_TABSTOP | WS_VISIBLE | WS_CHILD | ES_LEFT | ES_READONLY
                                  | ES_AUTOVSCROLL | ES_AUTOHSCROLL | ES_MULTILINE | WS_VSCROLL | WS_HSCROLL,
                                  14,
                                  600,
                                  984,
                                  90,
                                  hWnd,
                                  NULL,
                                  app->hInst,
                                  NULL);

    SetWindowFont(app->start_button, app->control_font, TRUE);
    SetWindowFont(app->console, app->console_font, TRUE);

    struct MatmulBenchTest *ts = app->bench->test_set;
    int nt = app->bench->num_test;

    app->tests = (struct bench_test *)malloc(sizeof(struct bench_test) * nt);

    HWND l0 = CreateWindowEx(0,
                             "STATIC",
                             "max sec",
                             WS_VISIBLE | WS_CHILD | SS_RIGHT,
                             40,
                             20 + 40,
                             150,
                             25,
                             hWnd,
                             NULL,
                             app->hInst,
                             NULL);

    app->sec_in = CreateWindowEx(WS_EX_CLIENTEDGE,
                                 "EDIT",
                                 "1.5",
                                 WS_VISIBLE | WS_CHILD | ES_LEFT,
                                 240,
                                 20 + 40,
                                 150,
                                 25,
                                 hWnd,
                                 NULL,
                                 app->hInst,
                                 NULL);

    HWND l1 = CreateWindowEx(0,
                             "STATIC",
                             "iter",
                             WS_VISIBLE | WS_CHILD | SS_RIGHT,
                             40,
                             20 + 40 + 24,
                             150,
                             25,
                             hWnd,
                             NULL,
                             app->hInst,
                             NULL);
    app->iter_in = CreateWindowEx(WS_EX_CLIENTEDGE,
                                  "EDIT",
                                  "3",
                                  WS_VISIBLE | WS_CHILD | ES_LEFT,
                                  240,
                                  20 + 40 + 24,
                                  150,
                                  25,
                                  hWnd,
                                  NULL,
                                  app->hInst,
                                  NULL);

    SetWindowFont(l0, app->control_font, TRUE);
    SetWindowFont(l1, app->control_font, TRUE);
    SetWindowFont(app->sec_in, app->control_font, TRUE);
    SetWindowFont(app->iter_in, app->control_font, TRUE);

    for (int i=0; i<nt; i++) {
        DWORD s = WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON | BS_CHECKBOX | BS_AUTOCHECKBOX;
        HWND b = CreateWindow("BUTTON",
                              ts[i].name,
                              s,
                              20,
                              20 + 100 + 30*i,
                              150,
                              25,
                              hWnd,
                              NULL,
                              app->hInst,
                              NULL);

        app->tests[i].check_box = b;
        SetWindowFont(b, app->control_font, TRUE);

        Button_SetCheck(b, 0);

        s = WS_VISIBLE | WS_CHILD | ES_LEFT | ES_READONLY;
        HWND e = CreateWindowEx(WS_EX_CLIENTEDGE,
                                "EDIT",
                                "n/a",
                                s,
                                240,
                                20 + 100 + 30*i,
                                150,
                                25,
                                hWnd,
                                NULL,
                                app->hInst,
                                NULL);

        app->tests[i].result_box = e;
        app->tests[i].max_flops = 0;
        SetWindowFont(e, app->control_font, TRUE);
    }

    return 1;
}

static int
on_size(HWND hWnd, UINT state, int cx, int cy)
{
    return 0;
}

static void
update_state(struct app *app,
             enum app_state st)
{
    app->st = st;

    if (st == APP_RUNNING) {
        Button_Enable(app->start_button, 0);
        Button_SetText(app->start_button, "‚¯‚¢‚»‚­‚¿‚ã..");
    } else {
        Button_Enable(app->start_button, 1);
        Button_SetText(app->start_button, "Start");
    }
}

static void
disp_flops(const struct MatmulBenchTest *test,
           double sec,
           unsigned int iter,
           unsigned long mat_size,
           void *ap)
{
    struct app *app = (struct app*)ap;
    double n = mat_size;

    npr_strbuf_printf(&app->console_text,
                      "\r\n(%5d-%2d):%-20s: sec=%8.5f, %9.5f[GFLOPS], %8.5f[GB/s]",
                      (int)mat_size,
                      (int)iter,
                      test->name,
                      sec,
                      n*n*n*2/(sec*1000.0*1000.0*1000.0),
                      (n*n*3.0*sizeof(float))/(sec*1000.0*1000.0*1000.0));

    char *con = npr_strbuf_c_str(&app->console_text);
    Edit_SetText(app->console, con);
    Edit_SetSel(app->console, app->console_text.cur, app->console_text.cur);
    Edit_ScrollCaret(app->console);
}


static unsigned __stdcall
bench_func(void *ap)
{
    struct app *app = (struct app*)ap;
    struct MatmulBenchResult *r;
    int rti;

    r = matmul_bench_run(app->bench,
                         app->config,
                         disp_flops,
                         app);

    int size_step = r->run_size_step;
    int size_min = r->run_size_min;

    for (rti=0; rti<r->num_test; rti++) {
        int ti = r->test_map[rti];

        double max = app->tests[ti].max_flops;
        struct MatmulBenchTestResult *tr = &r->results[rti];
        int config_iter = app->config->iter;
        int num_run = tr->num_run;
        int ri, ii;
        int max_size = 0;

        for (ri=0; ri<num_run; ri++) {
            for (ii=0; ii<config_iter; ii++) {
                int n = size_min + size_step * ri;
                double sec = tr->sec[ri][ii];
                double flops = matmul_bench_calc_gflops(n, sec);

                if (flops > max) {
                    max = flops;
                    max_size = n;
                }
            }
        }

        char flops_text[64];
        snprintf(flops_text, 64, "%.3f@%d^3", max, max_size);

        Edit_SetText(app->tests[ti].result_box, flops_text);
    }

    if (app->last_result) {
        matmul_bench_result_fini(app->bench, app->last_result);
    }
    app->last_result = r;

    SetEvent(app->finish_event);

    return 0;
}


static void
start_benchmark(HWND hWnd)
{
    GET_APP();

    int num_test = app->bench->num_test;
    struct MatmulBenchConfig *config = app->config;

    for (int ti=0; ti<num_test; ti++) {
        int enable = Button_GetState(app->tests[ti].check_box);
        config->enable[ti] = enable;
    }

    config->size_step = 64;
    config->size_min = 64;

    char text[4096];

    Edit_GetText(app->iter_in, text, 4096);
    int iter = atoi(text);

    if (iter <= 0 || iter > 10) {
        iter = 3;
        Edit_SetText(app->iter_in, "3");
    }
    config->iter = iter;

    Edit_GetText(app->sec_in, text, 4096);
    double sec = atof(text);
    if (sec <= 0 || sec > 120) {
        sec = 1.5;
        Edit_SetText(app->sec_in, "1.5");
    }

    config->max_time_sec = sec;

    app->console_text.cur = 0;

    update_state(app, APP_RUNNING);

    unsigned int threadID;
    app->run_thread = (HANDLE)_beginthreadex(NULL, 0, bench_func, app, 0, &threadID);
}

static void
on_command(HWND hWnd, int id, HWND com_wnd, UINT code)
{
    GET_APP();

    switch (id) {
    case IDM_ABOUT:
        MessageBox(NULL, "‚Ü‚Æ‚Ü‚é‚¿‚á‚ñ\nver 1.0", "‚Ü‚Æ‚Ü‚é‚¿‚á‚ñ", MB_OK);
        break;

    case IDM_EXPORT: {
        OPENFILENAME ofn;
        char szFile[1024];

        ZeroMemory(&ofn, sizeof(ofn));

        ofn.lStructSize = sizeof(ofn);

        ofn.hwndOwner = app->main_win;
        ofn.lpstrFile = szFile;
        ofn.lpstrFile[0] = '\0';
        ofn.nMaxFile = sizeof(szFile);
        ofn.lpstrFilter = "csv(*.csv)\0*.csv\0All(*.*)\0*.*\0\0";
        ofn.nFilterIndex = 1;
        ofn.lpstrFileTitle = NULL;
        ofn.nMaxFileTitle = 0;
        ofn.lpstrInitialDir = NULL;
        ofn.Flags = 0;
        ofn.lpstrDefExt = "csv";

        if (GetOpenFileName(&ofn)) {
            
        }
    }
        break;


    case IDC_START:
        start_benchmark(hWnd);
        break;
    }
}

static HBRUSH
on_ctl_color_static(HWND hwnd, HDC hdc, HWND ctl, int message)
{
    return GetStockObject(NULL_BRUSH);
}

static void
on_paint(HWND hwnd)
{
    PAINTSTRUCT ps;

    HDC dc = BeginPaint(hwnd, &ps);

    FillRect(dc, &ps.rcPaint, (HBRUSH)(COLOR_WINDOW+1));

    EndPaint(hwnd, &ps);
}

static LRESULT WINAPI
wndproc( HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam )
{
    switch (msg) {
        HANDLE_MSG(hWnd, WM_DESTROY, on_destroy);
        HANDLE_MSG(hWnd, WM_CREATE, on_create);
        HANDLE_MSG(hWnd, WM_SIZE, on_size);
        HANDLE_MSG(hWnd, WM_COMMAND, on_command);
        HANDLE_MSG(hWnd, WM_CTLCOLORSTATIC, on_ctl_color_static);
        HANDLE_MSG(hWnd, WM_PAINT, on_paint);

    default:
        return DefWindowProc( hWnd, msg, wParam, lParam );
    }
}


static int
main2(HINSTANCE hInst, int nCmdShow)
{
    {
        WNDCLASSEX wc = {
            sizeof( WNDCLASSEX ), CS_CLASSDC, wndproc, 0L, 0L,
            hInst, NULL, LoadCursor(NULL,IDC_ARROW), NULL, MAKEINTRESOURCE(IDC_MAINMENU),
            TEXT("matmulbench"), NULL
        };
        RegisterClassEx(&wc);
    }

    struct app *app = (struct app*)malloc(sizeof(struct app));
    app->st = APP_IDLE;
    app->hInst = hInst;
    app->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
    HWND hWnd = CreateWindow( TEXT("matmulbench"),
			      TEXT("matmulbench"),
                              WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, 1024, 768,
                              NULL, NULL, hInst, app );

    ShowWindow(hWnd, nCmdShow);
    UpdateWindow(hWnd);

    MSG msg;
    int run = 1;

    while (run) {
        DWORD r = MsgWaitForMultipleObjects(1,
                                            &app->finish_event,
                                            FALSE,
                                            INFINITE,
                                            QS_ALLINPUT);
        if (r == WAIT_OBJECT_0) {
            WaitForSingleObject(app->run_thread, INFINITE);
            CloseHandle(app->run_thread);

            update_state(app, APP_IDLE);
        } else {
            while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))  {
                if (msg.message == WM_QUIT) {
                    run = 0;
                } else {
                    TranslateMessage(&msg);
                    DispatchMessage(&msg);
                }
            }
        }
    }

    if (app->last_result) {
        matmul_bench_result_fini(app->bench, app->last_result);
    }
    matmul_bench_fini(app->bench);

    free(app);

    return msg.wParam;
}




int WINAPI
_tWinMain(HINSTANCE hInst, HINSTANCE hPrev,
	  LPTSTR cmdline,
	  INT nCmdShow)
{
    return main2(hInst, nCmdShow);
}

int
main()
{
    return main2(GetModuleHandle(NULL), SW_SHOW);
}
