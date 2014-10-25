#include <stdlib.h>
#include <string.h>
#include "xstdint.h"
#include "npr/mempool.h"
#include "npr/bits.h"

#define MT_ARG MEMPOOL_MT_ARG
#ifdef GATHER_STAT
#define PASS_MT , mt
#define MT_OR_ZERO mt
#else
#define PASS_MT
#define MT_OR_ZERO 0
#endif

const char *memtype_name[] = {
#define DEFINE_MEMTYPE(name) #name,
#include "npr/memtype.def"
#undef DEFINE_MEMTYPE
};

void
npr_mempool_init(struct npr_mempool *pool,
             unsigned int size_hint)
{
    pool->entry_num = 16;
    pool->data_entry = (unsigned char**)malloc(sizeof(unsigned char*)*pool->entry_num);
    pool->entry_index = 0;
    pool->entry_byte_remain = pool->entry_byte_size = size_hint;
    pool->entry_byte_pos = 0;
    pool->data_entry[0] = (unsigned char*)malloc(size_hint);

    pool->large_num = 16;
    pool->large_index = 0;
    pool->large = (unsigned char **)malloc(sizeof(void*)*pool->large_num);

    pool->alloc_small = 0;
#ifdef GATHER_STAT
    pool->alloc_total = 0;
    pool->unused_mem = 0;
    pool->align_padding = 0;
    memset(pool->size_each_type, 0, sizeof(pool->size_each_type));
#endif
}

#ifdef GATHER_STAT
#define INC_STAT(name,val) pool->name += val
#else
#define INC_STAT(name,val) 
#endif

void
npr_mempool_destroy(struct npr_mempool *pool)
{
    int i;
    int n = pool->entry_index;
    for (i=0; i<=n; i++) {
        free(pool->data_entry[i]);
    }

    free(pool->data_entry);

    n = pool->large_index;
    for (i=0; i<n; i++) {
        free(pool->large[i]);
    }
    free(pool->large);
}

void
npr_mempool_clear(struct npr_mempool *p)
{
    if (p->entry_index == 0 && p->large_index == 0) {
        p->entry_byte_pos = 0;
        return;
    }

    npr_mempool_destroy(p);
    npr_mempool_init(p, roundup2(p->alloc_small));
}

static void *
alloc_large(struct npr_mempool *pool,
            unsigned int size
            MT_ARG
)
{
    void *ret;
    if (pool->large_index+1 >= pool->large_num) {
        pool->large_num *= 2;
        pool->large = (unsigned char**)realloc(pool->large, pool->large_num * sizeof(char*));
    }

    ret = malloc(size);
    pool->large[pool->large_index] = (unsigned char*)ret;
    pool->large_index++;

    INC_STAT(alloc_total, size);
    INC_STAT(size_each_type[mt], size);

    return ret;
}

void *
npr_mempool_alloc_align_stat(struct npr_mempool *pool,
                             unsigned int align_shift, unsigned int size
                             MT_ARG
)
{
    uintptr_t ptr, ptr_aligned;
    int ptr_diff, remain;
    unsigned char *alloc_ptr;
    unsigned int alloc_size;
    unsigned int align;
    uintptr_t align_mask;
    unsigned int aligned_size;


    if (size >= 1024) {
        return alloc_large(pool, size PASS_MT);
    }

    align = 1<<align_shift;
    align_mask = ~((uintptr_t)(1<<align_shift) - 1);
    aligned_size = (size + (align-1))&align_mask;

retry:
    ptr = (uintptr_t)(pool->data_entry[pool->entry_index] + pool->entry_byte_pos);
    ptr_aligned = (ptr+(align-1))&align_mask;
    ptr_diff = ptr_aligned - ptr;
    remain = pool->entry_byte_remain - ptr_diff;
    alloc_ptr = (unsigned char*)ptr_aligned;

    if ((int)size > remain) {
        unsigned int next_entry_index;
        INC_STAT(unused_mem, pool->entry_byte_remain);

        /* realloc entry */
        next_entry_index = pool->entry_index + 1;
        alloc_size = pool->entry_byte_size*2;

        if ((size+align) > alloc_size)
            alloc_size = (size+align)*2;

        if (next_entry_index >= pool->entry_num) {
            /* realloc entry pointer */
            int next_entry_num = pool->entry_num * 2;
            pool->data_entry = (unsigned char**)realloc(pool->data_entry, next_entry_num * sizeof(unsigned char*));
            pool->entry_num = next_entry_num;
        }

        pool->entry_byte_remain = pool->entry_byte_size = alloc_size;
        pool->entry_index = next_entry_index;
        pool->entry_byte_pos = 0;

        pool->data_entry[next_entry_index] = (unsigned char*)malloc(alloc_size);

        goto retry;
    }
    INC_STAT(align_padding, ptr_diff);

    alloc_size = (ptr_diff + aligned_size);

    remain -= alloc_size;
    pool->entry_byte_remain = remain;
    pool->entry_byte_pos += alloc_size;

    pool->alloc_small += alloc_size;
    INC_STAT(alloc_total, alloc_size);
    INC_STAT(size_each_type[mt], alloc_size);

    /* memset(alloc_ptr, 0x33, size); */

    return alloc_ptr;
}

char *
npr_mempool_strdup(struct npr_mempool *p,
                   const char *str)
{
    int len = strlen(str);
    char *mem = (char*)npr_mempool_alloc(p, len+1);
    memcpy(mem, str, len);
    mem[len] = '\0';
    return mem;
}

char *
npr_mempool_strndup(struct npr_mempool *p,
                    const char *str,
                    int len)
{
    char *mem = (char*)npr_mempool_alloc(p, len+1);
    memcpy(mem, str, len);
    mem[len] = '\0';
    return mem;
}

void *
npr_mempool_copy(struct npr_mempool *p,
                 void *data,
                 int len
                 MT_ARG)
{
    void *mem = npr_mempool_alloc_memtype(p, len, MT_OR_ZERO);
    memcpy(mem, data, len);
    return mem;
}

