#ifndef NPR_COMMON_BITS_HPP
#define NPR_COMMON_BITS_HPP

#include "xstdint.h"
#include "npr/compiler.h"

extern const unsigned char npr_popcnt_tbl[256];

static __inline unsigned int
roundup2(unsigned int x) {
 x = x - 1;
 x = x | (x >> 1);
 x = x | (x >> 2);
 x = x | (x >> 4);
 x = x | (x >> 8);
 x = x | (x >>16);
 return x + 1;
}

static __inline unsigned int
npr_popcnt32(uint32_t v)
{
    return npr_popcnt_tbl[v>>24] +
        npr_popcnt_tbl[(v>>16)&0xff] +
        npr_popcnt_tbl[(v>>8)&0xff] +
        npr_popcnt_tbl[(v>>0)&0xff];
}

static __inline unsigned int
npr_popcnt64(uint64_t v)
{
    return npr_popcnt32(v>>32) + npr_popcnt32(v&0xffffffff);
}

ALWAYS_INLINE
static unsigned int
npr_bsf32(unsigned int v)
{
#ifdef __GNUC__
    unsigned int ret;
    __asm__ ("bsf %1, %0"
             :"=r"(ret)
             :"rm"(v));
    return ret;
#else
    unsigned long m = v;
    unsigned long idx;
    _BitScanForward(&idx, m);
    return idx;
#endif
}

ALWAYS_INLINE
static unsigned int
npr_bsr32(unsigned int v)
{
#ifdef __GNUC__
    unsigned int ret;
    __asm__ ("bsr %1, %0"
             :"=r"(ret)
             :"g"(v));
    return ret;
#else
    unsigned long m = v;
    unsigned long idx;
    _BitScanReverse(&idx, m);
    return idx;
#endif
}


static __inline unsigned int
npr_bsf64(uint64_t v)
{
    if (v&0xffffffff) {
        return npr_bsf32(v&0xffffffff);
    } else {
        return npr_bsf32(v>>32) + 32;
    }
}

static __inline int
npr_find_bit_idx32_keymask(uint32_t bits,
                           uint32_t key)
{
    return npr_popcnt32((key-1)&bits);
}

static __inline int
npr_find_bit_idx32(uint32_t bits,
                   uint32_t pos)
{
    return npr_find_bit_idx32_keymask(bits, (1<<pos));
}


static __inline int
npr_is_x2(uint32_t bits)
{
    return (((bits-1) & bits) == 0);
}


#endif
