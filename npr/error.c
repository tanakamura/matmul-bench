#include "npr/error.h"

#ifdef _WIN32
void
npr_errno_message(struct npr_strbuf *sb,
                  npr_errno_t err)
{
    LPVOID lpMsgBuf;
    FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER |
        FORMAT_MESSAGE_FROM_SYSTEM |
        FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        err,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPTSTR) &lpMsgBuf,
        0,
        NULL);

    npr_strbuf_puts(sb, (char*)lpMsgBuf);

    LocalFree(lpMsgBuf);
}

#else

#include <string.h>
void
npr_errno_message(struct npr_strbuf *sb,
                  npr_errno_t err)
{
    const char *p = strerror(err);
    npr_strbuf_puts(sb, p);
}

#endif
