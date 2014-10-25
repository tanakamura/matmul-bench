#include <stdlib.h>
#include <string.h>
#include "npr/symbol.h"

struct npr_symbol *npr_plus_symbol, *npr_minus_symbol;

struct backet {
    struct backet *chain;
    struct npr_symbol *value;
};

struct table_t {
    unsigned int num_backets;
    struct backet **backets;
};

static struct table_t table;

static int
hash(const char *string, int len)
{
    unsigned int hval = NPR_SYMBOL_FNV1_32A_INIT;
    int i;

    /*
     * FNV-1a hash each octet in the buffer
     */
    for (i=0; i<len; i++) {
        npr_symbol_hash_push(&hval, string[i]);
    }
    return hval;
}

struct npr_symbol *
npr_intern_with_hash( const char * symstr, size_t str_len, unsigned int hval )
{
    unsigned int h = hval % table.num_backets;
    struct backet *chain, **begin = &(table.backets[h]);
    struct backet *b;
    struct npr_symbol *sym;

    chain = *begin;

    while ( chain ) {
        if ((chain->value->symstr_len == str_len) &&
            (memcmp(chain->value->symstr,symstr,str_len) == 0))
            return chain->value;

        chain = chain->chain;
    }

    sym = malloc( sizeof(struct npr_symbol) );
    sym->symstr = malloc( str_len+1 );
    memcpy( sym->symstr, symstr, str_len );
    sym->symstr[ str_len ] = '\0';
    sym->symstr_len = str_len;
    sym->hashcode = hval;
    sym->keyword = 0;

    sym->var_value = NULL;
    sym->tag_value = NULL;

    b = malloc( sizeof(struct backet) );
    b->value = sym;
    b->chain = table.backets[h];
    table.backets[h] = b;

    return sym;
}

struct npr_symbol *
npr_intern_with_length( const char *symstr, size_t str_len )
{
    unsigned int hval = hash(symstr,str_len);
    return npr_intern_with_hash(symstr, str_len, hval);
}

struct npr_symbol *
npr_intern( const char *symstr )
{
    return npr_intern_with_length( symstr, strlen(symstr) );
}

const char *
npr_intern_str(const char *symstr)
{
    struct npr_symbol *sym = npr_intern(symstr);
    return sym->symstr;
}

void
npr_symbol_init( void )
{
    static int init = 0;
    int i;
    if (init) {
        return;
    }

    init = 1;
    table.num_backets = 173;
    table.backets = (struct backet**)malloc( sizeof(struct backet*) * 173 );
    for ( i=0; i<173; i++ ) {
        table.backets[i] = NULL;
    }

#define KW(s,t) npr_intern( s )->keyword = t
}


void
npr_symbol_finish( void )
{
    int i;
    for ( i=0; i<173; i++ ) {
        struct backet *p = table.backets[i], *next;
        while ( p ) {
            next = p->chain;

            free( p->value->symstr );
            free( p->value );
            free( p );
            p = next;
        }
    }

    free( table.backets );
}
