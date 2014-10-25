#include "npr/varray.h"
#include "npr/mempool.h"
#include <stdlib.h>
#include <string.h>

void
npr_varray_init(struct npr_varray *a,
		   size_t n,
		   size_t elem_size)
{
    size_t sz;
    a->nelem = 0;
    a->size = n;
    a->elem_size = elem_size;
    sz = elem_size * n;

    a->elements = malloc(sz);
}
void
npr_varray_init_pool(struct npr_varray *a,
			size_t n,
			size_t elem_size,
			struct npr_mempool *p)
{
    size_t sz;
    a->nelem = 0;
    a->size = n;
    a->elem_size = elem_size;
    sz = elem_size * n;

    a->elements = npr_mempool_alloc(p, sz);
}

void
npr_varray_realloc(struct npr_varray *a)
{
    a->size *= 2;
    a->elements = realloc(a->elements, a->size * a->elem_size);
}

void
npr_varray_realloc_pool(struct npr_varray *a, struct npr_mempool *p)
{
    void *old = a->elements;
    size_t sz = a->size * a->elem_size;
    a->size *= 2;
    a->elements = npr_mempool_alloc(p, a->size * a->elem_size);
    memcpy(a->elements, old, sz);
}

void
npr_varray_resize(struct npr_varray *a, int n)
{
    a->nelem = n;
    if (n <= a->size) {
        return;
    }
    a->size = n*2;
    a->elements = realloc(a->elements, a->size * a->elem_size);
}

void
npr_varray_resize_pool(struct npr_varray *a, int n, struct npr_mempool *p)
{
    void *old;
    size_t sz;

    a->nelem = n;

    if (n <= a->size) {
        return;
    }

    old = a->elements;
    sz = a->size * a->elem_size;

    a->size = n*2;
    a->elements = npr_mempool_alloc(p, a->size * a->elem_size);
    memcpy(a->elements, old, sz);
}


void *
npr_varray_copy(struct npr_varray *a,
                struct npr_mempool *p)
{
    size_t sz = a->nelem * a->elem_size;
    void *ret;

    ret = npr_mempool_alloc(p, sz);
    memcpy(ret, a->elements, sz);

    return ret;
}

void *
npr_varray_malloc_copy(struct npr_varray *a)
{
    size_t sz = a->nelem * a->elem_size;
    void *ret;

    ret = malloc(sz);
    memcpy(ret, a->elements, sz);

    return ret;
}

void *
npr_varray_close(struct npr_varray *a,
		    struct npr_mempool *p)
{
    void *ret = npr_varray_copy(a, p);
    free(a->elements);
    return ret;
}

void *
npr_varray_malloc_close(struct npr_varray *a)
{
    void *ret = npr_varray_malloc_copy(a);
    free(a->elements);
    return ret;
}

void
npr_varray_discard(struct npr_varray *a)
{
    free(a->elements);
}
		    
