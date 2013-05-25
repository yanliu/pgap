#ifndef LOG_H
#define LOG_H

#define LOG_BUF_SIZE 1024*64
#define LOG_ENTRY_BUF_SIZE 1024*16
#define LOG_INTERVAL 20

#include "addr.h"
void glog_init(char *d, char *p, IDTYPE i);
void glog(char *format, ...);


#endif
