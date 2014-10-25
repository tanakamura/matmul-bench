#ifndef NPR_MEMPOOL_H
#define NPR_MEMPOOL_H

#ifdef __cplusplus
extern "C" {
#endif

enum npr_mem_type {
#define DEFINE_MEMTYPE(name) NPR_MEM_##name,
#include "npr/memtype.def"
#undef DEFINE_MEMTYPE
    NUM_NPR_MEM_TYPE
};

extern const char *memtype_name[];

#ifdef GATHER_STAT
#define MEMPOOL_MT_ARG     , enum npr_mem_type mt
#else
#define MEMPOOL_MT_ARG
#endif


/**
 * @file
 * @brief メモリアロケータ
 *
 * 関連するオブジェクトをまとめて解放したい場合に使うアロケータ
 */

/**
 * @brief メモリアロケータ
 */
struct npr_mempool {
    unsigned int entry_num;
    unsigned int entry_index;
    unsigned int entry_byte_size;
    unsigned int entry_byte_remain;
    unsigned int entry_byte_pos;
    unsigned char **data_entry;

    unsigned int large_num;
    unsigned int large_index;
    unsigned char **large;

    int alloc_small;

#ifdef GATHER_STAT
    int alloc_total;
    int unused_mem;
    int align_padding;
    int size_each_type[NUM_NPR_MEM_TYPE];
#endif
};

/**
 * @brief create mempool
 */
extern void npr_mempool_init(struct npr_mempool *p,
                         unsigned int size_hint);


/**
 * @brief delete mempool
 */
extern void npr_mempool_destroy(struct npr_mempool *p);

extern void npr_mempool_clear(struct npr_mempool *p);

/**
 * @brief allocate memory from pool.
 *        memory is aligend to `2^align'.
 * 
 * @param p pool object
 * @param align_shift align
 * @param size object size
 * @return allocated memory
 */
extern void *npr_mempool_alloc_align_stat(struct npr_mempool *p,
                                          unsigned int align_shift,
                                          unsigned int size
                                          MEMPOOL_MT_ARG
);

#ifdef GATHER_STAT
#define npr_mempool_alloc_align(p,a,s,t) npr_mempool_alloc_align_stat(p,a,s,t)
#else
#define npr_mempool_alloc_align(p,a,s,t) ((void)(t),npr_mempool_alloc_align_stat(p,a,s))
#endif

/**
 * @brief allocate memory from pool.
 *        memory is aligned to 8.
 * @param p pool
 * @param s size
 * @return allocated memory
 */
#define npr_mempool_alloc(p,s) (npr_mempool_alloc_align(p, 3, s, NPR_MEM_OTHER))

#define npr_mempool_alloc_memtype(p,s,t) (npr_mempool_alloc_align(p, 3, s, t))

extern char *npr_mempool_strdup(struct npr_mempool *p,
                                const char *str
);
extern char *npr_mempool_strndup(struct npr_mempool *p,
                                 const char *str,
                                 int n
);

extern void *npr_mempool_copy(struct npr_mempool *p,
                              void *data,
                              int data_len
                              MEMPOOL_MT_ARG);


#ifdef __cplusplus
}
#endif

#endif
