#include <stdlib.h>
#include "npr/list.h"

void
npr_singly_list_init(struct npr_singly_list *l)
{
    l->head = NULL;
    l->tail_chain = &l->head;
}
void
npr_singly_list_fini(struct npr_singly_list *l, int free_elem)
{
    struct npr_singly_list_elem *e, *n;
    e = l->head;
    while (e) {
        n = e->chain;
        if (free_elem) {
            free(e->value);
        }
        free(e);
        e = n;
    }
}

void
npr_singly_list_push_head(struct npr_singly_list *l, void *ptr)
{
    struct npr_singly_list_elem *e = malloc(sizeof(*e));

    if (l->head == NULL) {
        l->tail_chain = &e->chain;
    }

    e->value = ptr;
    e->chain = l->head;
    l->head = e;
}

void npr_singly_list_push_tail(struct npr_singly_list *l, void *ptr)
{
    struct npr_singly_list_elem *e = malloc(sizeof(*e));
    *l->tail_chain = e;
    l->tail_chain = &e->chain;
    e->value = ptr;
    e->chain = NULL;
}

void
npr_dlist_init(struct npr_dlist *l)
{
    npr_chunk_allocator_init(&l->a, sizeof(struct npr_dlist_elem), 8);
    l->head.next = &l->tail;
    l->tail.prev = &l->head;
}

void
npr_dlist_fini(struct npr_dlist *l)
{
    npr_chunk_allocator_fini(&l->a);
}

struct npr_dlist_elem *
npr_dlist_push_back(struct npr_dlist *l,
                    void *val)
{
    struct npr_dlist_elem *e = npr_chunk_allocator_alloc(&l->a);
    e->next = &l->tail;
    e->prev = l->tail.prev;
    l->tail.prev = e;
    e->prev->next = e;

    e->val = val;

    return e;
}

void
npr_dlist_remove(struct npr_dlist *l,
                 struct npr_dlist_elem *e)
{
    e->prev->next = e->next;
    e->next->prev = e->prev;

    npr_chunk_allocator_free(&l->a, e);
}