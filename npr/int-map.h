#ifndef NPR_INT_MAP_H
#define NPR_INT_MAP_H

#include <stdio.h>
#include "npr/symbol.h"

#ifdef __cplusplus
extern "C" {
#endif

struct npr_symtab_entry {
    struct npr_symbol *sym;
    void *data;
    struct npr_symtab_entry *chain;
};

struct npr_symtab {
    int num_bin;
    int num_entry;
    struct npr_symtab_entry **entries;
};

void npr_symtab_init(struct npr_symtab *m,
                     int size_hint);
void npr_symtab_fini(struct npr_symtab *m);

enum npr_lookup_command{
    NPR_LOOKUP_APPEND,
    NPR_LOOKUP_FAIL
};

struct npr_symtab_entry *npr_symtab_lookup_entry(struct npr_symtab *tab,
                                                 struct npr_symbol *sym,
                                                 enum npr_lookup_command com);

void npr_symtab_stat(FILE *out,
                     struct npr_symtab *tab,
                     int verbose);

void npr_symtab_global_init(void);
void npr_symtab_global_fini(void);

#ifdef __cplusplus
}
#endif

#endif
