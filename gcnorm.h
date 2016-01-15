#ifndef __GCNORM
#define __GCNORM

#include <stdio.h>
#include <string.h>

#include "globals.h"
#include "sam.h"
#include <math.h>
#include <float.h>

#define MAXCOPY 2.5
#define MAXCOPY_X 1.5
#define MINCOPY 1.5
#define MINCOPY_X 0.5


int norm_until_converges(enum WINDOWTYPE, float *, float *);
void normalize(enum WINDOWTYPE, float *, float *, float *, float *, float *, float *);
void calc_stat(enum WINDOWTYPE, float *, float *, char, float *, float *, float *, float *);
float clean_outlier(enum WINDOWTYPE, float *, float *);
void norm_wins(int, enum WINDOWTYPE, float *, float *);

#endif
