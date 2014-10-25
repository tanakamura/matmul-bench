#ifndef NPR_STAT_H
#define NPR_STAT_H

#include "npr/error.h"

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

int npr_get_filesize(const char *path);
npr_errno_t npr_read_file(
    unsigned int *ret_size,
    char **ret_buf,
    const char *path);

#ifdef _WIN32
#define access _access
#endif

#ifdef __cplusplus
}
#endif

#endif
