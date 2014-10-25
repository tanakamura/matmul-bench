#ifndef NPR_VARRAY_H
#define NPR_VARRAY_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct npr_mempool;

struct npr_varray {
    size_t nelem;               /**< number of elements */
    size_t size;                /**< size of buffer */
    size_t elem_size;           /**< size of a element */
    void *elements;             /**< elements buffer */
};

#define VA_ELEM(t,a,n) (((t*)((a)->elements))[n])
#define VA_TOP(t,a) VA_ELEM(t,a,(a)->nelem-1)
#define VA_ELEM_PTR(t,a,n) (&(((t*)((a)->elements))[n]))
#define VA_LAST_PTR(t,a) (&(((t*)((a)->elements))[(a)->nelem-1]))
#define VA_PUSH(t,a,e) do { if ((a)->nelem>=(a)->size) npr_varray_realloc(a); VA_ELEM(t,a,(a)->nelem++) = (e); } while(0)
#define VA_NEWELEM(a) do { if ((a)->nelem>=(a)->size) npr_varray_realloc(a); (a)->nelem++; } while(0)
#define VA_NEWELEM_LASTPTR(t,a,s) do { VA_NEWELEM(a); s = VA_LAST_PTR(t,a); } while (0)
#define VA_NEWELEM_P(a,p) do { if ((a)->nelem>=(a)->size) npr_varray_realloc_pool(a,p); (a)->nelem++; } while(0)
#define VA_POP(t,a) VA_ELEM(t, (a), --((a)->nelem))
#define VA_PUSH_P(t,a,e,p) do { if ((a)->nelem>=(a)->size) npr_varray_realloc_pool(a,p); VA_ELEM(t,a,(a)->nelem++) = (e); } while(0)
void npr_varray_init(struct npr_varray *a, size_t n, size_t es);
void npr_varray_init_pool(struct npr_varray *a, size_t n, size_t es,
                                struct npr_mempool *pool);
void npr_varray_realloc(struct npr_varray *a);
void npr_varray_realloc_pool(struct npr_varray *a,
                             struct npr_mempool *pool);
void npr_varray_resize(struct npr_varray *a, int nelem);
void npr_varray_resize_pool(struct npr_varray *a,
                            int nelem,
                            struct npr_mempool *pool);
void *npr_varray_close(struct npr_varray *a,
                       struct npr_mempool *p);
void *npr_varray_malloc_close(struct npr_varray *a);
void *npr_varray_copy(struct npr_varray *a,
                      struct npr_mempool *p);
void *npr_varray_malloc_copy(struct npr_varray *a);
void npr_varray_discard(struct npr_varray *a);

#ifdef __cplusplus
}
#endif

#endif
