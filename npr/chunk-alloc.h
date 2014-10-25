#ifndef NPR_CHUNK_ALLOC_H
#define NPR_CHUNK_ALLOC_H
#include <stdio.h>
#ifdef __cplusplus
extern "C" {
#endif

struct npr_chunk_allocator_chunk;

struct npr_chunk_allocator {
    void *top;

    int elem_size;
    int chunk_count;

    struct npr_chunk_allocator_chunk *chunks;
};

void npr_chunk_allocator_init(struct npr_chunk_allocator *a,
                              int elem_size,
                              int chunk_count);
void npr_chunk_allocator_fini(struct npr_chunk_allocator *a);

void *npr_chunk_allocator_alloc(struct npr_chunk_allocator *a);
void npr_chunk_allocator_free(struct npr_chunk_allocator *a, void *data);

void npr_chunk_allocator_stat(FILE *out, struct npr_chunk_allocator *a);
                              

#ifdef __cplusplus
}
#endif

#endif
