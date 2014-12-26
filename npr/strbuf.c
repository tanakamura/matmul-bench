#include "npr/strbuf.h"
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include <stdio.h>

static void
reserve(struct npr_strbuf *buf, int len)
{
    size_t d = buf->buflen - buf->cur;
    if (d < len) {
        size_t req = buf->buflen + len;
        size_t alloc = req*2;
        buf->buflen = alloc;
        buf->buf = realloc(buf->buf, alloc);
    }
}

static void
append(struct npr_strbuf *buf, int c)
{
    buf->buf[buf->cur++] = c;
}

void
npr_strbuf_init(struct npr_strbuf *buf)
{
    buf->buf = NULL;
    buf->buflen = 0;
    buf->cur = 0;
}
void
npr_strbuf_fini(struct npr_strbuf *buf)
{
    free(buf->buf);
}

char *
npr_strbuf_strdup(struct npr_strbuf *buf)
{
    size_t len = buf->cur;
    char *ret = malloc(len + 1);
    memcpy(ret, buf->buf, len);
    ret[len] = '\0';
    return ret;
}

char *
npr_strbuf_strdup_pool(struct npr_strbuf *buf, struct npr_mempool *pool)
{
    size_t len = buf->cur;
    char *ret = npr_mempool_alloc(pool, len + 1);
    memcpy(ret, buf->buf, len);
    ret[len] = '\0';
    return ret;
}


void
npr_strbuf_putc(struct npr_strbuf *buf,
                int c)
{
    reserve(buf, 1);
    append(buf, c);
}

void
npr_strbuf_puts(struct npr_strbuf *buf,
                const char *s)
{
    size_t len = strlen(s);
    reserve(buf, len);
    memcpy(&buf->buf[buf->cur], s, len);
    buf->cur += len;
}
void
npr_strbuf_putsn(struct npr_strbuf *buf,
                 const char *s, int len)
{
    reserve(buf, len);
    memcpy(&buf->buf[buf->cur], s, len);
    buf->cur += len;
}

void
npr_strbuf_vprintf(struct npr_strbuf *buf,
                   const char *format,
                   va_list l0)
{
    int rem;
    int len;
    va_list tmp;

    /* should use va_copy on C99 */
#ifdef _MSC_VER
    tmp = l0;
#else
    __va_copy(tmp, l0);
#endif

#ifdef _WIN32
    (void)rem;
    len = _vscprintf(format, tmp);
    va_end(tmp);

#else
    rem = buf->buflen - buf->cur;
    len = vsnprintf(buf->buf + buf->cur, rem, format, tmp);
    va_end(tmp);

    if (len < rem) {
        buf->cur += len;
        return;
    }
#endif

    reserve(buf, len+1);
    vsnprintf(buf->buf + buf->cur, len+1, format, l0);
    buf->cur += len;
}

void
npr_strbuf_printf(struct npr_strbuf *buf,
                  const char *fmt,
                  ...)
{
    va_list l0;
    va_start(l0, fmt);
    npr_strbuf_vprintf(buf, fmt, l0);
    va_end(l0);
}


char *
npr_strbuf_c_str(struct npr_strbuf *buf)
{
    reserve(buf, 1);
    buf->buf[buf->cur] = '\0';

    return buf->buf;
}
