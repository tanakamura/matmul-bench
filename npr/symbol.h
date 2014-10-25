#ifndef NPR_SYMBOL_H
#define NPR_SYMBOL_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct cmdkscript_var;

#define NPR_SYMBOL_FNV_32_PRIME 0x01000193
#define NPR_SYMBOL_FNV1_32A_INIT 0x811c9dc5

static __inline void
npr_symbol_hash_push(unsigned int *hashval, char c)
{
    unsigned int hval = *hashval;
    hval ^= (unsigned int)c;
    hval *= NPR_SYMBOL_FNV_32_PRIME;
    *hashval = hval;
}

struct npr_symbol {
    char *symstr;
    unsigned int symstr_len;
    unsigned int hashcode;
    unsigned int keyword;

    struct cmdkscript_var *var_value;
    struct cmdkscript_var *tag_value;

    void *native_code;
    void *value;
};

struct npr_symbol *npr_intern( const char * symstr );
struct npr_symbol *npr_intern_with_hash( const char * symstr, size_t strlen, unsigned int hash );
const char *npr_intern_str( const char * symstr );
struct npr_symbol *npr_intern_with_length( const char * symstr, size_t strlen );

void npr_symbol_init(void);
void npr_symbol_finish( void );


#ifdef __cplusplus
}
#endif



#endif
