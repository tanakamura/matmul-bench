#ifndef NPR_ERROR_H
#define NPR_ERROR_H

#include "npr/strbuf.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <errno.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#include <windows.h>

typedef DWORD npr_errno_t;

static __inline npr_errno_t
npr_errno() {
    return GetLastError();
}

static __inline int
npr_errno_ok(npr_errno_t e) {
    return (e == ERROR_SUCCESS);
}



#else
typedef int npr_errno_t;

static __inline npr_errno_t
npr_errno() {
    return errno;
}

static __inline int
npr_errno_ok(npr_errno_t e) {
    return (e == 0);
}

#endif

void npr_errno_message(struct npr_strbuf *sb, npr_errno_t err);

#ifdef __cplusplus
}
#endif

#endif
