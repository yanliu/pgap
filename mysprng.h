#ifndef MYSPRNG_H
#define MYSPRNG_H

#include "sprng.h"

#ifdef PGAMODE
#define USE_MPI
#endif

double mysprng();
void mysprng_init();
void mysprng_free();

#endif
