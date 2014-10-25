#include <stdlib.h>
#include "npr/free-chain.h"

struct dtor {
    void (*f)(void*);
    void *data;
};

void
npr_free_chain_init(struct npr_free_chain *fc)
{
    npr_varray_init(&fc->frees, 16, sizeof(void*));
    npr_varray_init(&fc->dtors, 16, sizeof(struct dtor));
}
void
npr_free_chain_append(struct npr_free_chain *fc,
                      void *ptr)
{
    VA_PUSH(void *, &fc->frees, ptr);
}
void *
npr_free_chain_append_malloc(struct npr_free_chain *fc,
                             size_t sz)
{
    void *ret = malloc(sz);
    npr_free_chain_append(fc, ret);
    return ret;
}


void
npr_free_chain_append_dtor(struct npr_free_chain *fc,
                           void (*dtor)(void*),
                           void *ptr)
{
    struct dtor *d;
    VA_NEWELEM_LASTPTR(struct dtor, &fc->dtors, d);
    d->data = ptr;
    d->f = dtor;
}

void
npr_free_chain_free_all(struct npr_free_chain *fc)
{
    int n = fc->frees.nelem, i;
    for (i=0; i<n; i++) {
        free(VA_ELEM(void *, &fc->frees, i));
    }
    npr_varray_discard(&fc->frees);

    n = fc->dtors.nelem;
    for (i=0; i<n; i++) {
        struct dtor *d = &VA_ELEM(struct dtor, &fc->dtors, i);
        d->f(d->data);
    }
    npr_varray_discard(&fc->dtors);
}
void
npr_free_chain_close(struct npr_free_chain *fc)
{
    npr_varray_discard(&fc->frees);
    npr_varray_discard(&fc->dtors);
}
