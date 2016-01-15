#ifndef __PREP
#define __PREP

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "globals.h"
#include "utils.h"

typedef struct gapcell{
  char *chrom;
  int start;
  int end;
}_gapcell;


struct gapcell *gaptable;
struct gapcell *pseudotable;

void prep_genome(void);
int count_chrom(FILE *, int *);
void read_ref(FILE *, int, int);
void insert_chrom(int, char *, int, char *, int);
void windowmaker(int, int, char *, char *, int, int);

#endif
