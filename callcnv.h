#ifndef __CALLCNV
#define __CALLCNV

#define FILTER 0.5

#include <stdio.h>
#include "globals.h"
#include "utils.h"
#include "gcnorm.h"

enum CALLTYPE {DUPS, DELS};

void call_cnv(char *, char *, char *);
void readDepth(char *);
void print_copy_numbers(char *);
void dump_text_windows(char *, enum WINDOWTYPE);
void conc_depth(char **, int, char *);
void calculate_dups(char *);
void calculate_dels(char *);
void calculate_intervals(char *, enum CALLTYPE);
void trim_sw(int ,int ,int, enum CALLTYPE);
void write_dw(char* ,char *);
void filter_dw(int ,int ,int, enum CALLTYPE);
void read_gene(char*, char *);
//int q_compare(const void *, const void *);

#endif
