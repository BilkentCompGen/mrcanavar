#include "utils.h"

void print_error(char *msg){
  fprintf(stderr, "\nmrCaNaVaR version %s.\nLast update: %s.\n", VERSION, LAST_UPDATE);

  fprintf(stderr, "\n%s\n", msg);
  if (!strstr(msg, "No run mode"))
    fprintf(stderr, "Kahrolsun bağzı şeyler.\n");
  fprintf(stderr, "Invoke parameter -h for help.\n");
  exit (0);
}

void set_runmode(enum MODETYPE NEWMODE){
  if (RUNMODE != NONE && RUNMODE != NEWMODE){
    print_error("Select one run mode only.\n");
  }
  RUNMODE = NEWMODE;
}

void set_gender(enum GENDERTYPE SET){
  if (GENDER != AUTODETECT && GENDER != SET){
    print_error("Select only one of --xx or --xy.\n");
  }
  GENDER = SET;
}

void set_str(char **target, char *source){

  if (*target != NULL) freeMem((*target), strlen(*target)+1);

  (*target) = (char *) getMem(sizeof(char) * (strlen(source)+1));

  strncpy((*target), source, (strlen(source)+1));

}


void init_globals(void){
  GENOME_FASTA = NULL;
  GENOME_CONF  = NULL;
  GENOME_GAPS  = NULL;
  GENOME_PSEUDO = NULL;
  RUNMODE      = NONE;
  GENDER       = AUTODETECT;
  MULTGC       = 0;
  MAX_GC_CORR  = 20.0;
  MIN_GC_CORR  = 0.05;
  VERBOSE      = 0;
  ISNORMALIZED = 0;

  num_chrom    = 0;
  memUsage     = 0;

  CHECKSAM     = 1;

  /* set window size defaults */

  LW_SIZE  = 5000;
  SW_SIZE  = 1000;
  CW_SIZE  = 1000;
  LW_SLIDE = 1000;
  SW_SLIDE = 1000;

  LW_MEAN  = 0.0;
  SW_MEAN  = 0.0;
  CW_MEAN  = 0.0;

  LW_STD   = 0.0;
  SW_STD   = 0.0;
  CW_STD   = 0.0;

  CONT_WINDOW = 7;
  CUT_WINDOW  = 6;

  MIN_DUP = 10000;

  PLOIDY = 2;
}


FILE *my_fopen(char *fname, char *mode, char GZ){
  FILE *fp;


  if (!GZ)
    fp = fopen(fname, mode);
  else
    fp = (FILE *) gzopen(fname, mode);

  if (fp == NULL){
    fprintf(stderr, "Unable to open file %s in %s mode.", fname, mode[0]=='w' ? "write" : "read");
    exit (0);
  }

  return fp;
}


void saveRefConfig(char *configFile){

  int i;
  int wcnt;
  char chrom_name_len;
  int retVal;
  FILE *config;

  config = my_fopen(configFile, "w", 0);

  /* sort the chromosomes pointer array */

  qsort(chromosomes, num_chrom, sizeof (struct chrom *), compChr);


  /* start with the magicNum, I will use this as a format check when loading */
  retVal = fwrite(&magicNum, sizeof(magicNum), 1, config);

  /* window sizes / slides */

  retVal = fwrite(&LW_SIZE,  sizeof(LW_SIZE),  1, config);
  retVal = fwrite(&SW_SIZE,  sizeof(SW_SIZE),  1, config);
  retVal = fwrite(&CW_SIZE,  sizeof(CW_SIZE),  1, config);
  retVal = fwrite(&LW_SLIDE, sizeof(LW_SLIDE), 1, config);
  retVal = fwrite(&SW_SLIDE, sizeof(SW_SLIDE), 1, config);


  /* reference genome numbers */


  /* number of pseudoautosomal regions */
  retVal = fwrite(&num_pseudo, sizeof(num_pseudo), 1, config);

  if (num_pseudo != 0){
    for (i=0; i<num_pseudo; i++){
      chrom_name_len = (char) strlen(pseudotable[i].chrom);
      retVal = fwrite(&chrom_name_len, sizeof(chrom_name_len), 1, config);
      retVal = fwrite(pseudotable[i].chrom, chrom_name_len * sizeof(char), 1, config);
      retVal = fwrite(&(pseudotable[i].start), sizeof(pseudotable[i].start), 1, config);
      retVal = fwrite(&(pseudotable[i].end), sizeof(pseudotable[i].end), 1, config);
    }
  }

  /* number of chromosomes */
  retVal = fwrite(&num_chrom, sizeof(num_chrom), 1, config);

  /* iterate through chromosomes, write their names, and window counts */

  for (i=0; i<num_chrom; i++){
    chrom_name_len = (char) strlen(chromosomes[i]->name);
    retVal = fwrite(&chrom_name_len, sizeof(chrom_name_len), 1, config);
    retVal = fwrite(chromosomes[i]->name, chrom_name_len * sizeof(char), 1, config);

    retVal = fwrite(&(chromosomes[i]->length), sizeof(chromosomes[i]->length), 1, config);
    retVal = fwrite(&(chromosomes[i]->lw_cnt), sizeof(chromosomes[i]->lw_cnt), 1, config);
    retVal = fwrite(&(chromosomes[i]->sw_cnt), sizeof(chromosomes[i]->sw_cnt), 1, config);
    retVal = fwrite(&(chromosomes[i]->cw_cnt), sizeof(chromosomes[i]->cw_cnt), 1, config);

    /* long windows */
    for (wcnt=0; wcnt<chromosomes[i]->lw_cnt; wcnt++){
      retVal = fwrite(&(chromosomes[i]->lw[wcnt].start), sizeof(chromosomes[i]->lw[wcnt].start), 1, config);
      retVal = fwrite(&(chromosomes[i]->lw[wcnt].end),   sizeof(chromosomes[i]->lw[wcnt].end),   1, config);
      retVal = fwrite(&(chromosomes[i]->lw[wcnt].gc),    sizeof(chromosomes[i]->lw[wcnt].gc),    1, config);
    }

    /* short windows */
    for (wcnt=0; wcnt<chromosomes[i]->sw_cnt; wcnt++){
      retVal = fwrite(&(chromosomes[i]->sw[wcnt].start), sizeof(chromosomes[i]->sw[wcnt].start), 1, config);
      retVal = fwrite(&(chromosomes[i]->sw[wcnt].end),   sizeof(chromosomes[i]->sw[wcnt].end),   1, config);
      retVal = fwrite(&(chromosomes[i]->sw[wcnt].gc),    sizeof(chromosomes[i]->sw[wcnt].gc),    1, config);
    }

    /* copy windows */
    for (wcnt=0; wcnt<chromosomes[i]->cw_cnt; wcnt++){
      retVal = fwrite(&(chromosomes[i]->cw[wcnt].start), sizeof(chromosomes[i]->cw[wcnt].start), 1, config);
      retVal = fwrite(&(chromosomes[i]->cw[wcnt].end),   sizeof(chromosomes[i]->cw[wcnt].end),   1, config);
      retVal = fwrite(&(chromosomes[i]->cw[wcnt].gc),    sizeof(chromosomes[i]->cw[wcnt].gc),    1, config);
    }

  }

  fclose(config);

  if (VERBOSE && retVal == 0)
    fprintf(stderr, "There was potentially an error in fwrite during saveConfig.\n");

}




void loadRefConfig(char *configFile){

  int i;
  int wcnt;
  char chrom_name_len;
  char readString[MAX_STR];
  int retVal;

  int isMagicNum;

  FILE *config;

  config = my_fopen(configFile, "r", 0);

  /* start with the magicNum,  use this as a format check when loading */
  retVal = fread(&isMagicNum, sizeof(isMagicNum), 1, config);

  if (isMagicNum != magicNum)
    print_error("Reference configuration file seems to be invalid or corrupt.\nYou might have created the configuration file with an earlier version of mrCaNaVaR.\n");


  fprintf(stdout, "Loading reference configuration, hold on ... ");
  fflush(stdout);

  /* window sizes / slides */

  retVal = fread(&LW_SIZE,  sizeof(LW_SIZE),  1, config);
  retVal = fread(&SW_SIZE,  sizeof(SW_SIZE),  1, config);
  retVal = fread(&CW_SIZE,  sizeof(CW_SIZE),  1, config);
  retVal = fread(&LW_SLIDE, sizeof(LW_SLIDE), 1, config);
  retVal = fread(&SW_SLIDE, sizeof(SW_SLIDE), 1, config);


  /* reference genome numbers */

  /* number of pseudoautosomal regions */
  retVal = fread(&num_pseudo, sizeof(num_pseudo), 1, config);

  if (num_pseudo != 0){

    pseudotable = (struct gapcell *) getMem(sizeof(struct gapcell) * num_pseudo);

    for (i=0; i<num_pseudo; i++){
      retVal = fread(&chrom_name_len, sizeof(chrom_name_len), 1, config);

      retVal = fread(readString, chrom_name_len * sizeof(char), 1, config);
      readString[(int)chrom_name_len] = 0;
      trimspace(readString);

      set_str(&(pseudotable[i].chrom), readString);

      retVal = fread(&(pseudotable[i].start), sizeof(pseudotable[i].start), 1, config);
      retVal = fread(&(pseudotable[i].end), sizeof(pseudotable[i].end), 1, config);
    }
  }

  /* number of chromosomes */
  retVal = fread(&num_chrom, sizeof(num_chrom), 1, config);

  /* create chromosomes data structure */

  chromosomes = (struct chrom **) getMem(sizeof(struct chrom *) * num_chrom);
  for (i=0; i<num_chrom; i++)
    chromosomes[i] = (struct chrom *) getMem(sizeof(struct chrom) * num_chrom);

  /* iterate through chromosomes, write their names, and window counts */

  for (i=0; i<num_chrom; i++){

    retVal = fread(&chrom_name_len, sizeof(chrom_name_len), 1, config);

    retVal = fread(readString, chrom_name_len * sizeof(char), 1, config);
    readString[(int)chrom_name_len] = 0;
    trimspace(readString);

    set_str(&(chromosomes[i]->name), readString);

    retVal = fread(&(chromosomes[i]->length), sizeof(chromosomes[i]->length), 1, config);
    retVal = fread(&(chromosomes[i]->lw_cnt), sizeof(chromosomes[i]->lw_cnt), 1, config);
    retVal = fread(&(chromosomes[i]->sw_cnt), sizeof(chromosomes[i]->sw_cnt), 1, config);
    retVal = fread(&(chromosomes[i]->cw_cnt), sizeof(chromosomes[i]->cw_cnt), 1, config);

    /* create windows structures */


    chromosomes[i]->lw = (struct window *) getMem (sizeof(struct window) * chromosomes[i]->lw_cnt);
    chromosomes[i]->sw = (struct window *) getMem (sizeof(struct window) * chromosomes[i]->sw_cnt);
    chromosomes[i]->cw = (struct window *) getMem (sizeof(struct window) * chromosomes[i]->cw_cnt);
    chromosomes[i]->dw = (struct window *) getMem (sizeof(struct window) * chromosomes[i]->lw_cnt); //Duplicate Windows



    /* long windows */
    for (wcnt=0; wcnt<chromosomes[i]->lw_cnt; wcnt++){
      retVal = fread(&(chromosomes[i]->lw[wcnt].start), sizeof(chromosomes[i]->lw[wcnt].start), 1, config);
      retVal = fread(&(chromosomes[i]->lw[wcnt].end),   sizeof(chromosomes[i]->lw[wcnt].end),   1, config);
      retVal = fread(&(chromosomes[i]->lw[wcnt].gc),    sizeof(chromosomes[i]->lw[wcnt].gc),    1, config);
      // this is the default, auto-control-picker will fix this
      chromosomes[i]->lw[wcnt].isControl = 1;
      chromosomes[i]->lw[wcnt].depth = 0.0;
    }

    /* short windows */
    for (wcnt=0; wcnt<chromosomes[i]->sw_cnt; wcnt++){
      retVal = fread(&(chromosomes[i]->sw[wcnt].start), sizeof(chromosomes[i]->sw[wcnt].start), 1, config);
      retVal = fread(&(chromosomes[i]->sw[wcnt].end),   sizeof(chromosomes[i]->sw[wcnt].end),   1, config);
      retVal = fread(&(chromosomes[i]->sw[wcnt].gc),    sizeof(chromosomes[i]->sw[wcnt].gc),    1, config);
      chromosomes[i]->sw[wcnt].isControl = 1;
      chromosomes[i]->sw[wcnt].depth = 0.0;
    }

    /* copy windows */
    for (wcnt=0; wcnt<chromosomes[i]->cw_cnt; wcnt++){
      retVal = fread(&(chromosomes[i]->cw[wcnt].start), sizeof(chromosomes[i]->cw[wcnt].start), 1, config);
      retVal = fread(&(chromosomes[i]->cw[wcnt].end),   sizeof(chromosomes[i]->cw[wcnt].end),   1, config);
      retVal = fread(&(chromosomes[i]->cw[wcnt].gc),    sizeof(chromosomes[i]->cw[wcnt].gc),    1, config);
      chromosomes[i]->cw[wcnt].isControl = 1;
      chromosomes[i]->cw[wcnt].depth = 0.0;
    }

    /* fill in the rest of the chromosomes structure */

  }

  fprintf(stdout, "[OK]. %d chromosomes loaded.\n", num_chrom);

  fclose(config);

  if (VERBOSE && retVal == 0)
    fprintf(stderr, "There was potentially an error in fread during loadConfig.\n");

}



int compChr(const void *p1, const void *p2){
  /* compare function to sort the chromosome pointer array */
  struct chrom *a, *b;

  a = *((struct chrom **)p1);
  b = *((struct chrom **)p2);


  return (strcmp ( a->name, b->name) );

}


int endswith(char *src, char *end){


  int endlen = strlen(end);
  int srclen = strlen(src);
  int copyindex;


  if (endlen > srclen)
    return 0;

  copyindex = srclen - endlen;

  if (memcmp(src+copyindex, end, endlen))
    return 0;
  else
    return 1;


}

void trimspace(char *str){
  int len;
  char *cutstr;

  len = strlen(str) - 1;
  while(1){
    if (isspace(str[len]))
      str[len--]=0;
    else
      break;
  }
  
  cutstr = strchr(str, ' ');
  if (cutstr != NULL) 
    *cutstr = 0;
}

int isPseudo(struct window win, int chrom_id){
  int i;
  for (i=0 ; i < num_pseudo ; i++){
    if(!strcmp(chromosomes[chrom_id]->name, pseudotable[i].chrom))
      if(win.start >= pseudotable[i].start && win.end <= pseudotable[i].end)
	return 1;
  }
  return 0;
}

int isAutosome(struct window *win, int win_id, int chrom_id){
  int ret = 1;

  if (GENDER == FEMALE)
    return 1;

  if (PLOIDY == 1)
    return 1;
  
  if (win != NULL)
    if (isPseudo(win[win_id], chrom_id))
      return 1;


  if (strstr(chromosomes[chrom_id]->name, "chrX") || strstr(chromosomes[chrom_id]->name, "chrY"))
    ret=0;
  if (strstr(chromosomes[chrom_id]->name, "X_random") || strstr(chromosomes[chrom_id]->name, "Y_random"))
    ret=0;
  if (!strcmp(chromosomes[chrom_id]->name, "X") || !strcmp(chromosomes[chrom_id]->name, "Y"))
    ret=0;

  return ret;
}

void * getMem(size_t size){
  void *ret;
  ret = malloc(size);
  if (ret == NULL){
    fprintf(stderr, "Cannot allocate memory. Currently addressed memory = %0.2f MB, requested memory = %0.2f MB.\nCheck the available main memory, and if you have user limits (ulimit -v).\n", getMemUsage(), (float)(size/1048576.0));
    exit(0);
  }
  memUsage+=size;
  return ret;
}

double getMemUsage(){
  return memUsage/1048576.0;
}

void freeMem(void *ptr, size_t size){
  memUsage-=size;
  free(ptr);
}
