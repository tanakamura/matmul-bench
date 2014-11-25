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
    SetWindowLongPtr(hWnd, GWLP_USERDATA, (LONG_PTR)app);

    app->console_text = malloc(128);
    app->console_text_size = 128;

    app->bench = matmul_bench_init(0);
    app->config = matmul_bench_config_init(app->bench);

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

    for (int i=0; i<nt; i++) {
        DWORD s = WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON | BS_CHECKBOX | BS_AUTOCHECKBOX;
        HWND b = CreateWindow("BUTTON",
                              ts[i].name,
                              s,
                              20,
                              20 + 40 + 30*i,
                              150,
                              25,
                              hWnd,
                              NULL,
                              app->hInst,
                              NULL);

        app->tests[i].check_box = b;
        SetWindowFont(b, app->control_font, TRUE);

        Button_SetCheck(b, 1);

        s = WS_VISIBLE | WS_CHILD | ES_LEFT | ES_READONLY;
        HWND e = CreateWindowEx(WS_EX_CLIENTEDGE,
                                "EDIT",
                                "n/a",
                                s,
                                180,
                                20 + 40 + 30*i,
                                150,
                                25,
                                hWnd,
                                NULL,
                                app->hInst,
                                NULL);

        app->tests[i].result_box = e;
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
    char buffer0[4096];

    snprintf(buffer0, 4096,
             "(%5d-%2d):%-20s: sec=%8.5f, %9.5f[GFLOPS], %8.5f[GB/s]\n",
             (int)mat_size,
             (int)iter,
             test->name,
             sec,
             n*n*n*2/(sec*1000.0*1000.0*1000.0),
             (n*n*3.0*sizeof(float))/(sec*1000.0*1000.0*1000.0));

    size_t len0 = strlen(buffer0);
    size_t edit_len = Edit_GetTextLength(app->console);
    size_t len = edit_len + len0 + 1;
    char *buffer = malloc(len);

    Edit_GetText(app->console, buffer, edit_len+1);
    memcpy(buffer + edit_len, buffer0, len0);
    buffer[edit_len + len0] = '\0';
    printf("%d %d %s\n", edit_len, len0, buffer);
    puts("==");
    Edit_SetText(app->console, buffer0);

    free(buffer);
}


static unsigned __stdcall
bench_func(void *ap)
{
    struct app *app = (struct app*)ap;

    matmul_bench_run(app->bench,
                     app->config,
                     disp_flops,
                     app);

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

    update_state(app, APP_RUNNING);

    unsigned int threadID;
    app->run_thread = (HANDLE)_beginthreadex(NULL, 0, bench_func, app, 0, &threadID);
}

static void
on_command(HWND hWnd, int id, HWND com_wnd, UINT code)
{
    switch (id) {
    case IDM_ABOUT:
        MessageBox(NULL, "‚Ü‚Æ‚Ü‚é‚¿‚á‚ñ\nver 1.0", "‚Ü‚Æ‚Ü‚é‚¿‚á‚ñ", MB_OK);
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

static LRESULT WINAPI
wndproc( HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam )
{
    switch (msg) {
        HANDLE_MSG(hWnd, WM_DESTROY, on_destroy);
        HANDLE_MSG(hWnd, WM_CREATE, on_create);
        HANDLE_MSG(hWnd, WM_SIZE, on_size);
        HANDLE_MSG(hWnd, WM_COMMAND, on_command);
        HANDLE_MSG(hWnd, WM_CTLCOLORSTATIC, on_ctl_color_static);

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

    while (1) {
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
            int r = GetMessage(&msg, NULL, 0, 0);
            if (!r) {
                break;
            }
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
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
