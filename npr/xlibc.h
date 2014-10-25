#ifndef NPR_XLIBC_H
#define NPR_XLIBC_H

#ifdef __cplusplus
extern "C" {
#endif

#include <string.h>

#ifdef _MSC_VER
#define xstrdup _strdup
#else
#define xstrdup strdup
#endif

#ifdef __cplusplus
}
#endif

#endif
