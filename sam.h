#ifndef __SAM
#define __SAM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <dirent.h>

#include "globals.h"
#include "utils.h"
#include "gcnorm.h"
#include "callcnv.h"

void read_mapfiles(char *, char *, char);
void readSAM(char *, char *, char, long *, long *, long *);

int insert_read_lw(char *, int, char *, int, long *);
int insert_read_sw(char *, int, char *, int, long *);
int insert_read_cw(char *, int, char *, int, long *);

int findChrom(char *, char *, int);

int binSearch(char *);

int windowSearch(struct window *, int, int);

void saveDepth(char *, int);

#endif
