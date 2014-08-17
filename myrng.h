#ifndef MYRNG_H
#define MYRNG_H

// macro: random number generator
#ifdef GSL_SPRNG
#include <gsl/gsl_rng.h>
#include "gsl-sprng.h"
#include <gsl/gsl_randist.h>
#define MYRANDI(kkk) ((int)(gsl_rng_uniform_int(gsl_sprng_r, (unsigned long int)(kkk))))
#define MYRANDF() (gsl_rng_uniform(gsl_sprng_r))
#else
#ifdef SPRNG
#include "mysprng.h"
#define MYRANDI(kkk) ((int)(mysprng() * (kkk)))
#define MYRANDF() (mysprng())
#else
#define MYRANDI(kkk) ((int)(drand48() * (kkk)))
#define MYRANDF() (drand48())
#endif
#endif

#endif
