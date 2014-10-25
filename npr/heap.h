#ifndef NPR_FIRSTFIT_HEAP_H
#define NPR_FIRSTFIT_HEAP_H

#include "xstdint.h"
#include <stddef.h>
#include <stdio.h>

enum {
    NPR_LARGE_CHUNK = (1<<0)
};

struct npr_heap_page_header {
    struct npr_heap_page_header *next, *prev; /* 16 */
    int flags;                   /* 20 */
    int page_size;               /* 24 to free large chunk */
};
struct npr_free_chunk_header {
    size_t sz;
    struct npr_free_chunk_header *next, *prev; /* 24 */
};

struct npr_heap {
    int flags;
    int page_size;
    int chunk_size;
    int num_ptr_per_page;
    uintptr_t pfn_mask;
    struct npr_free_chunk_header free_head, free_tail;
    struct npr_heap_page_header page_head, page_tail;
};

void npr_heap_init(struct npr_heap *h, int is_exec);
void npr_heap_fini(struct npr_heap *h);
void *npr_heap_alloc(struct npr_heap *h, size_t n);
void npr_heap_free(struct npr_heap *h, void *p, size_t n);
void npr_heap_dump(FILE *fp, struct npr_heap *h);

#endif
