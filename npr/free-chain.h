#ifndef NPR_FREE_CHAIN_H
#define NPR_FREE_CHAIN_H

#include "npr/varray.h"

#ifdef __cplusplus
extern "C" {
#endif

struct npr_free_chain {
    struct npr_varray frees;
    struct npr_varray dtors;
};

void npr_free_chain_init(struct npr_free_chain *fc);
void npr_free_chain_append(struct npr_free_chain *fc,
                           void *ptr);
void *npr_free_chain_append_malloc(struct npr_free_chain *fc,
                                   size_t sz);
void npr_free_chain_append_dtor(struct npr_free_chain *fc,
                                void (*dtor)(void*),
                                void *ptr);
void npr_free_chain_free_all(struct npr_free_chain *fc);
void npr_free_chain_close(struct npr_free_chain *fc);

#ifdef __cplusplus
}
#endif

#endif
