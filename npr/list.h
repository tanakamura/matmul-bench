#ifndef NPR_LIST_H
#define NPR_LIST_H

#include "npr/chunk-alloc.h"

#ifdef __cplusplus
extern "C" {
#endif

struct npr_singly_list_elem {
    struct npr_singly_list_elem *chain;
    void *value;
};

struct npr_singly_list {
    struct npr_singly_list_elem *head;
    struct npr_singly_list_elem **tail_chain;
};

void npr_singly_list_init(struct npr_singly_list *l);
void npr_singly_list_fini(struct npr_singly_list *l, int free_elem);
void npr_singly_list_push_head(struct npr_singly_list *l, void *ptr);
void npr_singly_list_push_tail(struct npr_singly_list *l, void *ptr);

#define NPR_SINGLY_LIST_FOR_EACH(NSlist,NStype,NSiter)  \
{                                                       \
    NStype NSiter;                                      \
    struct npr_singly_list_elem *NSle;                  \
    NSle = (NSlist)->head;                              \
    while (NSle) {                                      \
        NSiter = (NStype)NSle->value;

#define NPR_SINGLY_LIST_END_FOR_EACH()                  \
        NSle = NSle->chain;                             \
    }                                                   \
}



struct npr_dlist_elem {
    struct npr_dlist_elem *next;
    struct npr_dlist_elem *prev;
    void *val;
};

struct npr_dlist {
    struct npr_dlist_elem head, tail;
    struct npr_chunk_allocator a;
};

void npr_dlist_init(struct npr_dlist *l);
void npr_dlist_fini(struct npr_dlist *l);

struct npr_dlist_elem *npr_dlist_push_back(struct npr_dlist *l,
                                         void *val);

void npr_dlist_remove(struct npr_dlist *l,
                     struct npr_dlist_elem *e);

#define NPR_DLIST_FOR_EACH(NSlist,NStype,NSiter)           \
    {                                                     \
    struct npr_dlist_elem *npr_cur = (NSlist)->head.next;     \
    struct npr_dlist_elem *npr_end = &(NSlist)->tail;         \
    NStype NSiter;                                          \
    for (;npr_cur!=npr_end; npr_cur = npr_cur->next) {      \
    NSiter = (NStype)npr_cur->val;

#define NPR_DLIST_END_FOR_EACH()                \
    }}


#ifdef __cplusplus
}
#endif

#endif
