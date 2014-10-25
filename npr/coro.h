#ifndef NPR_CORO_H
#define NPR_CORO_H

#ifdef __cplusplus
extern "C" {
#endif


struct npr_coro {
    void *sp;
    void *save_sp;

    void *entry;
    void *stack_top;
    int stack_alloc_size;

#if (defined HAVE_VALGRIND) && (defined DEBUG)
    int stackid;
#endif
};

typedef void * (*npr_coro_entry_t)(struct npr_coro *self,
                                   void *);

static const int NPR_CORO_DISABLE_STACK_GUARD = (1<<0);

int npr_coro_init(struct npr_coro *coro,
                  npr_coro_entry_t entry,
                  void *entry_arg,
                  int stack_size,
                  int flags);

void npr_coro_fini(struct npr_coro *coro);

void *npr_coro_yield(struct npr_coro *self, void *arg);
void *npr_coro_cont(struct npr_coro *coro, void *arg);

#ifdef __cplusplus
}
#endif


#endif
