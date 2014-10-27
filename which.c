#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <windows.h>
#include <Shlwapi.h>

static char *(ext_list[512]);
static char *(path_list[4096]);


int
main(int argc, char **argv)
{
    char *path = strdup(getenv("PATH"));
    char *pathext = strdup(getenv("PATHEXT"));
    char *path1, *pathext1;
    int num_path=0, num_ext=0,pi,ei;

    if (argc != 2) {
        puts("usage : which aa.exe");
        return 1;
    }

    path_list[0] = "";
    num_path = 1;

    path1 = strtok(path, ";");

    while (path1 != NULL) {
        size_t len = strlen(path1);
        int last_slash = 0;
        char *p0;

        if (path1[len-1] == '\\' || path1[len-1] == '/') {
            p0 = malloc(len+1);
            last_slash = 1;
        } else {
            p0 = malloc(len+2);
        }

        memcpy(p0, path1, len+1);

        if (! last_slash) {
            p0[len] = '\\';
            p0[len+1] = '\0';
        }

        path_list[num_path++] = p0;

        path1 = strtok(NULL,";");
    }


    ext_list[0] = "";
    num_ext=1;

    pathext1 = strtok(pathext, ";");
    while (pathext1 != NULL) {
        ext_list[num_ext++]=strdup(pathext1);
        pathext1 = strtok(NULL,";");
    }

    for (pi=0; pi<num_path; pi++) {
        char *buffer = malloc(1024*1024);
        char *p;
        char *path = path_list[pi];

        for (ei=0; ei<num_ext; ei++) {
            char *ext = ext_list[ei];

            p = strcpy(buffer,path);
            p = strcat(p, argv[1]);
            p = strcat(p, ext);

            if (PathFileExists(p)) {
                DWORD len = GetFullPathName(p, 0, NULL, NULL);
                char *ret = malloc(len);
                GetFullPathName(p, len, ret, NULL);
                puts(ret);
                return 0;
            }
        }
    }


    return 1;
}
