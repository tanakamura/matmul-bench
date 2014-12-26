#ifndef NPR_STRBUF_H
#define NPR_STRBUF_H

#include <stdarg.h>
#include "npr/mempool.h"

#ifdef __cplusplus
extern "C" {
#endif

struct npr_strbuf {
    int buflen;
    int cur;
    char *buf;
};

void npr_strbuf_init(struct npr_strbuf *sb);
void npr_strbuf_fini(struct npr_strbuf *sb);

char *npr_strbuf_c_str(struct npr_strbuf *sb);

char *npr_strbuf_strdup(struct npr_strbuf *sb);
char *npr_strbuf_strdup_pool(struct npr_strbuf *sb, struct npr_mempool *pool);

void npr_strbuf_putc(struct npr_strbuf *sb, int c);
void npr_strbuf_puts(struct npr_strbuf *sb, const char *s);
void npr_strbuf_putsn(struct npr_strbuf *sb, const char *s, int text_len);
void npr_strbuf_vprintf(struct npr_strbuf *sb, const char *fmt, va_list l);
void npr_strbuf_printf(struct npr_strbuf *buf, const char *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif
