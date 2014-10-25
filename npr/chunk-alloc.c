#include "npr/chunk-alloc.h"
#include <stdlib.h>

struct npr_chunk_allocator_elem {
    struct npr_chunk_allocator_elem *chain;    
};

struct npr_chunk_allocator_chunk {
    void *data;
    struct npr_chunk_allocator_chunk *chain;
};

static struct npr_chunk_allocator_chunk *
alloc_chunk(int elem_size,
            int chunk_count)
{
    char *buf = malloc(elem_size * chunk_count);
    struct npr_chunk_allocator_chunk *c;

    int i;
    typedef struct npr_chunk_allocator_elem elem_t;
    elem_t *last;
    for (i=0; i<chunk_count-1; i++) {
        elem_t *e = (elem_t *)(buf + i*elem_size);
        elem_t *next = (elem_t *)(buf + (i+1)*elem_size);
        e->chain = next;
    }

    last = (elem_t*)(buf + ((chunk_count-1)*elem_size));
    last->chain = NULL;

    c = malloc(sizeof(*c));
    c->chain = NULL;
    c->data = buf;

    return c;
}

void
npr_chunk_allocator_init(struct npr_chunk_allocator *a,
                         int elem_size,
                         int chunk_count)
{
    struct npr_chunk_allocator_chunk *c;

    c = alloc_chunk(elem_size, chunk_count);

    a->top = c->data;
    a->chunks = c;
    a->elem_size = elem_size;
    a->chunk_count = chunk_count;
}

void *
npr_chunk_allocator_alloc(struct npr_chunk_allocator *a) {
    if (a->top) {
        struct npr_chunk_allocator_elem *ret = a->top, *next;
        next = ret->chain;
        a->top = next;
        return ret;
    } else {
        struct npr_chunk_allocator_chunk *c;
        struct npr_chunk_allocator_elem *e;
        c = alloc_chunk(a->elem_size, a->chunk_count);

        c->chain = a->chunks;
        a->chunks = c;

        e = c->data;
        a->top = e->chain;

        return e;
    }
}

void
npr_chunk_allocator_free(struct npr_chunk_allocator *a,
                         void *data)
{
    struct npr_chunk_allocator_elem *e = data;
    e->chain = a->top;
    a->top = e;
}

void
npr_chunk_allocator_fini(struct npr_chunk_allocator *a)
{
    struct npr_chunk_allocator_chunk *c, *n;
    c = a->chunks;
    while (c) {
        n = c->chain;
        free(c->data);
        free(c);
        c = n;
    }
}
                         
void
npr_chunk_allocator_stat(FILE *out,
                         struct npr_chunk_allocator *a)
{
    int chunk_count = 0;
    int free_count = 0;
    int total;

    struct npr_chunk_allocator_chunk *c;
    struct npr_chunk_allocator_elem *e;

    c = a->chunks;
    while (c) {
        c = c->chain;
        chunk_count ++;
    }

    total = chunk_count * a->chunk_count;

    e = a->top;
    while (e) {
        e = e->chain;
        free_count++;
    }

    fprintf(out,
            "total: %d(chunk=%d), used:%d, free:%d\n",
            total,
            chunk_count,
            total-free_count,
            free_count);
}
