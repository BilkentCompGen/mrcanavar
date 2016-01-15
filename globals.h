
#ifndef __GLOBALS
#define __GLOBALS


#define MAX_STR 256
#define GC_BIN  1000

#define VERSION "1.0-beta"
#define LAST_UPDATE "December 29, 2014\nHer yer Taksim, her yer direniş"


enum MODETYPE {NONE, PREP, READSAM, CALL, CONC};

enum WINDOWTYPE {LW, SW, CW};

enum GENDERTYPE {AUTODETECT, MALE, FEMALE};

enum GENDERTYPE GENDER;

enum MODETYPE RUNMODE;

int ISNORMALIZED;

int num_chrom;
int num_pseudo;
int n_sam_b; // Begining sam file's index
int n_sam_e; // Terminating sam file's index

double memUsage;

int MULTGC;
float MAX_GC_CORR;
float MIN_GC_CORR;
int VERBOSE;
int CHECKSAM;

extern FILE *pseudo;
extern struct gapcell *pseudotable;
static const int magicNum = 74088;
static const int magicNumDepth = 1764;

char *GENOME_FASTA;
char *GENOME_GAPS;
char *GENOME_CONF;
char *GENOME_PSEUDO;

int LW_SIZE;
int SW_SIZE;
int CW_SIZE;
int LW_SLIDE;
int SW_SLIDE;

float LW_MEAN;
float LW_STD;
float LW_MEAN_X;
float LW_STD_X;

float SW_MEAN;
float SW_STD;
float SW_MEAN_X;
float SW_STD_X;

float CW_MEAN;
float CW_STD;
float CW_MEAN_X;
float CW_STD_X;

int CONT_WINDOW;
int CUT_WINDOW;
int MIN_DUP;

typedef struct window{
  int start;
  int end;
  float gc;
  float depth;
  char isControl;
}_window;

typedef struct chrom{
  char *name;
  int  length;
  int lw_cnt;
  int sw_cnt;
  int cw_cnt;
  int dw_cnt;
  struct window *sw;
  struct window *lw;
  struct window *cw;
  struct window *dw;
}_chrom;


struct chrom **chromosomes;

#endif
