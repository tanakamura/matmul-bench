#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <sys/unistd.h>
#include <sys/mman.h>
#endif

#include "xstdint.h"
#include "npr/coro.h"
#include "npr/align.h"
#ifdef HAVE_VALGRIND
#include <valgrind/valgrind.h>
#endif

void npr_coro_done(void);
void npr_coro_start_x8664(void);

int npr_coro_init(struct npr_coro *coro,
                  npr_coro_entry_t entry,
                  void *entry_arg,
                  int stack_size,
                  int flags)
{
    int page_size;
    int num_page;
    int r;
    void *page;
    uintptr_t page_val;
    uintptr_t stack_bottom;
#if (defined __x86_64__) || (defined _AMD64_)
    uint64_t *stack_start;
#else
    uint32_t *stack_start;    
#endif


#ifdef _WIN32
    DWORD old_protect;
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    page_size = si.dwPageSize;
#else
    page_size = sysconf(_SC_PAGESIZE);
#endif
    num_page = NPR_CEIL_DIV(stack_size, page_size);
    num_page++;

#ifdef _WIN32
    page = VirtualAlloc(NULL,
                        num_page * page_size,
                        MEM_COMMIT,
                        PAGE_NOACCESS);

    if (page == NULL) {
        return -1;
    }

    page_val = (uintptr_t)page;

    r = VirtualProtect((void*)(page_val+page_size),
                       (num_page-1)*page_size,
                       PAGE_READWRITE,
                       &old_protect);
    if (r == 0) {
        VirtualFree(page, 0, MEM_RELEASE);
        return -1;
    }
#else
    page = mmap(NULL,
                num_page * page_size,
                PROT_NONE,
                MAP_ANONYMOUS|MAP_PRIVATE,
                -1, 0);

    if (page == MAP_FAILED) {
        return -1;
    }

    page_val = (uintptr_t)page;

    r = mprotect((void*)(page_val+page_size),
                 (num_page-1)*page_size,
                 PROT_READ|PROT_WRITE);
    if (r == -1) {
        munmap(page, num_page*page_size);
        return -1;
    }
#endif

#if (defined HAVE_VALGRIND) && (defined DEBUG)
    VALGRIND_MALLOCLIKE_BLOCK(page_val, num_page*page_size, 0, 0);
    coro->stackid = VALGRIND_STACK_REGISTER(page_val+page_size,
                                            page_val+num_page*page_size);
    if (RUNNING_ON_VALGRIND) {
        printf("stackid = %d\n", coro->stackid);
    }
#endif

    stack_bottom = page_val + num_page*page_size;
    coro->stack_top = (void*)page_val;
    coro->stack_alloc_size = num_page*page_size;
    coro->entry = (void*)entry;

#if (defined __x86_64__) || (defined _AMD64_)

    /*
     * x86-64 (before npr_coro_start_x8664)
     *
     *        32
     *             -------
     *              coro_start_x8664 <- sp
     *        16   ------- <- should be aligned to 16
     *              self
     *              arg
     * base + size ------- <- aligned to 16
     *
     *    regparam : undefined
     *
     * x86-64 (after npr_coro_start_x8664)
     *
     *        32
     *             -------
     *              
     *              ret addr <- sp
     *        16   ------- <- should be aligned to 16
     *              self
     *              arg
     * base + size ------- <- aligned to 16
     *
     *                        linux  | win64
     *              self   <- rsi    |  rcx
     *              arg    <- rdi    |  rdx
     */
    stack_start = (uint64_t*)(stack_bottom - (3*8));
    coro->sp = stack_start;
    stack_start[0] = (uint64_t)npr_coro_start_x8664;
    stack_start[1] = (uint64_t)coro;
    stack_start[2] = (uint64_t)entry_arg;
#elif (defined __i386__) || (defined _WIN32)
    /*
     * x86-32:
     *
     *             entry    <- sp
     *             ret addr 
     *             -------- <- should be aligned to 16
     *               self
     *               arg
     *               empty
     *               empty
     * base + size -------- <- aligned to 16
     */

    stack_start = (uint32_t*)(stack_bottom - (6*4));
    coro->sp = stack_start;
    stack_start[0] = (uint32_t)entry;
    stack_start[1] = (uint32_t)npr_coro_done;
    stack_start[2] = (uint32_t)coro;
    stack_start[3] = (uint32_t)entry_arg;
#endif

    return 0;
}

void
npr_coro_fini(struct npr_coro *coro)
{
#ifdef _WIN32
    VirtualFree(coro->stack_top, 0, MEM_RELEASE);
#else
#if (defined HAVE_VALGRIND) && (defined DEBUG)
    VALGRIND_STACK_DEREGISTER(coro->stackid);
    VALGRIND_FREELIKE_BLOCK(coro->stack_top, 0);
#endif
    munmap(coro->stack_top, coro->stack_alloc_size);
#endif
}
