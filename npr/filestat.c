#include "npr/filestat.h"

#include <stdlib.h>

#ifdef _WIN32
#include <windows.h>

int
npr_get_filesize(const char *path)
{
    struct _stat32 st;
    int r = _stat32(path, &st);
    if (r == -1) {
        return -1;
    }

    return st.st_size;
}

npr_errno_t
npr_read_file(unsigned int *ret_size,
              char **ret_buf,
              const char *path)
{
    HANDLE h = CreateFile(path, GENERIC_READ, FILE_SHARE_READ,
                          NULL, OPEN_EXISTING, 0, NULL);
    DWORD sz;
    char *ret;

    if (h == INVALID_HANDLE_VALUE) {
        return GetLastError();
    }

    sz = GetFileSize(h);
    ret = malloc(sz + 1);

    *ret_size = sz;
    *ret_buf = ret;

    ReadFile(h, ret, sz, NULL, NULL);
    ret[sz] = '\0';

    return ERROR_SUCCESS;
}

#else
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

int
npr_get_filesize(const char *path)
{
    struct stat st;
    int r = stat(path, &st);
    if (r == -1) {
        return -1;
    }

    return st.st_size;
}

npr_errno_t
npr_read_file(unsigned int *ret_size,
              char **ret_buf,
              const char *path)
{
    int fd = open(path, O_RDONLY);
    struct stat st;
    char *ret;

    if (fd == -1) {
        return errno;
    }

    fstat(fd, &st);
    ret = malloc(st.st_size + 1);

    *ret_size = st.st_size;
    *ret_buf = ret;

    read(fd, ret, st.st_size);
    ret[st.st_size] = '\0';

    return 0;
}


#endif


