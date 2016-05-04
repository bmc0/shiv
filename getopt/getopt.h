#ifndef _GETOPT_H
#define _GETOPT_H

#ifdef __cplusplus
extern "C" {
#endif

extern char *optarg;
extern int optind, opterr, optopt;

int getopt(int argc, char **argv, const char *options);

#ifdef __cplusplus
}
#endif

#endif
