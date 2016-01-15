#ifndef __UTILS
#define __UTILS


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <ctype.h>

#include "globals.h"
#include "prep.h"



void print_error(char *);

void set_runmode(enum MODETYPE);
void set_gender(enum GENDERTYPE SET);

void set_str(char **, char *);

void init_globals(void);

FILE *my_fopen(char *, char *, char);

void saveRefConfig(char *);
void loadRefConfig(char *);

int endswith(char *, char *);
void trimspace(char *);
int isPseudo(struct window, int);
int isAutosome(struct window *, int, int);

int compChr(const void *, const void *);
void * getMem(size_t);
double getMemUsage(void);
void freeMem(void *, size_t);

#endif
