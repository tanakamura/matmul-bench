#include "npr/int-map.h"
#include "npr/chunk-alloc.h"
#include <stdlib.h>

#define MINSIZE 8
static struct npr_chunk_allocator npr_symtab_elem_allocator;

/* from st.c */
static const long primes[] = {
	8 + 3,
	16 + 3,
	32 + 5,
	64 + 3,
	128 + 3,
	256 + 27,
	512 + 9,
	1024 + 9,
	2048 + 5,
	4096 + 3,
	8192 + 27,
	16384 + 43,
	32768 + 3,
	65536 + 45,
	131072 + 29,
	262144 + 3,
	524288 + 21,
	1048576 + 7,
	2097152 + 17,
	4194304 + 15,
	8388608 + 9,
	16777216 + 43,
	33554432 + 35,
	67108864 + 15,
	134217728 + 29,
	268435456 + 3,
	536870912 + 11,
	1073741824 + 85,
	0
};

static int
new_size(int size)
{
    int i;

    int newsize;

    for (i = 0, newsize = MINSIZE;
	 i < (int )(sizeof(primes)/sizeof(primes[0]));
	 i++, newsize <<= 1)
    {
	if (newsize > size) return primes[i];
    }
    /* Ran out of polynomials */
    abort();
    return -1;			/* should raise exception */
}


static struct npr_symtab_entry *
alloc_entry(struct npr_symtab *tab, int *do_rehash)
{
    struct npr_symtab_entry *e = npr_chunk_allocator_alloc(&npr_symtab_elem_allocator);
    int ratio;
    tab->num_entry++;

    ratio = tab->num_entry/tab->num_bin;

    if (ratio > 3) {
        int old_size = tab->num_bin;
        unsigned int size = new_size(tab->num_bin + 1), i;
        struct npr_symtab_entry **new_entries = malloc(sizeof(struct npr_symtab_entry*) * size),
            **old_entries = tab->entries;
        *do_rehash = 1;

        for (i=0; i<size; i++) {
            new_entries[i] = NULL;
        }

        for (i=0; i<old_size; i++) {
            struct npr_symtab_entry *e = old_entries[i], *n;
            while (e) {
                struct npr_symbol *sym = e->sym;
                unsigned int new_idx = sym->hashcode % size;
                n = e->chain;

                e->chain = new_entries[new_idx];
                new_entries[new_idx] = e;

                e = n;
            }
        }

        free(old_entries);

        tab->num_bin = size;
        tab->entries = new_entries;
    } else {
        *do_rehash = 0;
    }

    return e;
}

struct npr_symtab_entry *
npr_symtab_lookup_entry(struct npr_symtab *tab,
                        struct npr_symbol *sym,
                        enum npr_lookup_command com)
{
    unsigned int hash = sym->hashcode;
    unsigned int idx = hash % tab->num_bin;
    int do_rehash = 0;
    struct npr_symtab_entry **ep = &tab->entries[idx], *e;

    while (*ep) {
        e = *ep;
        if (e->sym == sym) {
            return e;
        }
        ep = &e->chain;
    }

    if (com == NPR_LOOKUP_APPEND) {
        e = alloc_entry(tab, &do_rehash);
        e->sym = sym;
        e->data = NULL;

        if (do_rehash) {
            idx = hash % tab->num_bin;
            e->chain = tab->entries[idx];
            tab->entries[idx] = e;
        } else {
            e->chain = NULL;
            *ep = e;
        }

        return e;
    } else {
        return NULL;
    }
}

void
npr_symtab_init(struct npr_symtab *m,
                int size_hint)
{
    int size = new_size(size_hint);
    int i;

    m->entries = malloc(sizeof(struct npr_symtab_entry*) * size);

    m->num_bin = size;
    m->num_entry = 0;
    for (i=0; i<size; i++) {
        m->entries[i] = NULL;
    }
}

void
npr_symtab_fini(struct npr_symtab *m)
{
    int i, n = m->num_bin;
    for (i=0; i<n; i++) {
        struct npr_symtab_entry *e = m->entries[i], *n;
        while (e) {
            n = e->chain;
            npr_chunk_allocator_free(&npr_symtab_elem_allocator, e);
            e = n;
        }
    }

    free(m->entries);
}

void
npr_symtab_global_init()
{
    npr_chunk_allocator_init(&npr_symtab_elem_allocator,
                             sizeof(struct npr_symtab_entry),
                             512);
}

void
npr_symtab_global_fini()
{
    npr_chunk_allocator_fini(&npr_symtab_elem_allocator);
}

void
npr_symtab_stat(FILE *out,
                struct npr_symtab *tab,
                int verbose)
{
    int i;
    fprintf(out,
            "table(%p): num_bin=%d, num_entry=%d\n",
            tab,
            tab->num_bin,
            tab->num_entry);
    for (i=0; i<tab->num_bin; i++) {
        struct npr_symtab_entry *e = tab->entries[i];
        int len = 0;

        while (e) {
            len++;
            e = e->chain;
        }

        fprintf(out, "  entry%d: len=%d\n", i, len);
        if (verbose) {
            e = tab->entries[i];
            while (e) {
                fprintf(out, "    %s = %p\n",
                        e->sym->symstr,
                        e->data);
                e = e->chain;
            }
        }
    }
}

