#ifndef NPR_ALIGN_H
#define NPR_ALIGN_H

#define NPR_ALIGN_UP(x,to) (((x)+((to)-1U))&~((to)-1U))
#define NPR_ALIGN_DOWN(x,to) ((x)&~((to)-1U))

#define NPR_CEIL_DIV(num,denom) (((num)+(denom-1))/(denom))

#ifdef _WIN32
#define aligned_malloc _aligned_malloc
#define aligned_free _aligned_free
#else
#include <malloc.h>
#define aligned_malloc(s,a) memalign(a,s)
#define aligned_free free
#endif

#endif
