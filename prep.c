#include "prep.h"

void prep_genome(void){
  FILE *fasta;
  FILE *gaps;
  FILE *pseudo;

  int num_gaps;

  char chrom[MAX_STR];
  int start, end;
  int max_len;

  int i;

  if (GENOME_FASTA == NULL)
    print_error("Select genome fasta file (input) through -fasta parameter.\n");
  if (GENOME_GAPS == NULL)
    print_error("Select genome gaps file (input) through -gaps parameter.\n");
  if (GENOME_CONF == NULL)
    print_error("Select genome fasta file (output) through -conf parameter.\n");


  if (LW_SIZE <= 0)
    print_error ("Large window size (-lw_size) must be a poitive integer.\n");
  if (SW_SIZE <= 0)
    print_error ("Small window size (-sw_size) must be a poitive integer.\n");
  if (CW_SIZE <= 0)
    print_error ("Copy number window size (-cw_size) must be a poitive integer.\n");
  if (LW_SLIDE <= 0)
    print_error ("Large window slide (-lw_slide) must be a poitive integer.\n");
  if (SW_SLIDE <= 0)
    print_error ("Small window slide (-sw_slide) must be a poitive integer.\n");

  if (SW_SIZE >= LW_SIZE)
    print_error ("Small window size (-sw_size) should be smaller than large window size (-lw_size).\n");
  if (SW_SLIDE > LW_SIZE)
    print_error ("Small window slide (-sw_slide) should be smaller than or equal to large window slide (-lw_slide).\n");


  fasta = my_fopen(GENOME_FASTA, "r", 0);
  gaps = my_fopen(GENOME_GAPS, "r", 0);
  num_gaps = 0;


  /* initialize gap table */
  while (fscanf (gaps, "%s\t%d\t%d\n", chrom, &start, &end) > 0)
    num_gaps ++;

  gaptable = (struct gapcell *) getMem(sizeof(struct gapcell) * num_gaps);

  rewind(gaps);


  /* read gap table */

  i = 0;
  while (fscanf (gaps, "%s\t%d\t%d\n", chrom, &start, &end) > 0){
    trimspace(chrom);
    gaptable[i].chrom = NULL;
    set_str(&(gaptable[i].chrom), chrom);
    gaptable[i].start = start;
    gaptable[i].end = end;
    i++;
  }
  fclose(gaps);

  if (GENOME_PSEUDO != NULL){
    pseudo = my_fopen(GENOME_PSEUDO, "r",0);
    num_pseudo = 0;
    /* initialize pseudo table */
    while (fscanf (pseudo, "%s\t%d\t%d\n", chrom, &start, &end) > 0)
      num_pseudo ++;
    
    pseudotable = (struct gapcell *) getMem(sizeof(struct gapcell) * num_pseudo);
    
    rewind(pseudo);
    
    /* read pseudo table */
    
    i = 0;
    while (fscanf (pseudo, "%s\t%d\t%d\n", chrom, &start, &end) > 0){
      trimspace(chrom);
      pseudotable[i].chrom = NULL;
      set_str(&(pseudotable[i].chrom), chrom);
      pseudotable[i].start = start;
      pseudotable[i].end = end;
      i++;
    }
    fclose(pseudo);
    
  }
  
  fprintf (stdout, "Scanning the reference genome.\n");
  /* count basic numbers from the reference genome */
  max_len = 0;
  num_chrom = count_chrom(fasta, &max_len);
  rewind(fasta);

  chromosomes = (struct chrom **) getMem (sizeof (struct chrom *) * num_chrom);

  for (i=0; i<num_chrom; i++)
    chromosomes[i] = (struct chrom *) getMem (sizeof (struct chrom) * num_chrom);

  fprintf(stdout, "Total of %d chromosomes in the reference genome, with %d gaps.\nLongest chromosome is %d bp\n", num_chrom, num_gaps, max_len);

  fprintf (stdout, "Reading the reference genome.\n");

  read_ref(fasta, max_len, num_gaps);

  fprintf (stdout, "Saving the reference configuration.\n");
  saveRefConfig(GENOME_CONF);

}



int count_chrom(FILE *fasta, int *max_len){
  int length=0; int maxlength=0;
  char ch;
  int cnt = 0;
  char skip[MAX_STR];
  char *retStr = NULL;

  //while (fscanf(fasta, "%c", &ch) > 0){
  while (1){

    ch = fgetc(fasta);
    if (feof(fasta)) break;

    if (ch == '>') {
      cnt ++;
      retStr = fgets(skip, MAX_STR, fasta);
      fprintf(stdout, ".");
      fflush(stdout);
      if (length > maxlength) maxlength = length;
      length = 0;
    }
    else if (!isspace(ch)) length++;
  }
  fprintf(stdout, "\n");

  if (length > maxlength)
    maxlength = length;

  *max_len = maxlength;

  if (VERBOSE && retStr == NULL)
    fprintf(stderr, "There was potentially an error in fgets during count_chrom.\n");

  return cnt;
}



void read_ref(FILE *fasta, int max_len, int num_gaps){

  char ch;
  int cnt = 0;
  char chrom_name[MAX_STR];
  int length;
  char *chrom_seq;
  char *retStr = NULL;

  chrom_seq = (char *) getMem (sizeof(char) * (max_len+1));
  cnt = -1;
  length = 0;

  while (1){

    ch = fgetc(fasta);
    if (feof(fasta)) break;

    if (ch == '>') {

      if (cnt != -1){
	/* insert the previous chromosome */
	chrom_seq[length]=0;
	insert_chrom(cnt, chrom_name, length, chrom_seq, num_gaps);
      }

      cnt ++;
      retStr = fgets(chrom_name, MAX_STR, fasta);
      chrom_name[strlen(chrom_name)-1] = 0;
      trimspace(chrom_name);
      length = 0;
    }
    else if (!isspace(ch)){
      chrom_seq[length++] = ch;
    }
  }

  /* insert the last chromosome */
  chrom_seq[length]=0;
  trimspace(chrom_name);
  insert_chrom(cnt, chrom_name, length, chrom_seq, num_gaps);

  freeMem(chrom_seq, sizeof(char) * (max_len+1));

  if (VERBOSE && retStr == NULL)
    fprintf(stderr, "There was potentially an error in fgets during read_ref.\n");}


void insert_chrom(int cnt, char *chrom_name, int length, char *chrom_seq, int num_gaps){

  fprintf(stdout, "Chromosome %s (%d bp). Counting windows ... ", chrom_name, length);
  fflush(stdout);

  chromosomes[cnt]->name = NULL;
  set_str(&(chromosomes[cnt]->name), chrom_name);
  chromosomes[cnt]->length = length;

  windowmaker(num_gaps, cnt, chrom_name, chrom_seq, length, 0);
  fprintf(stdout, "\t\tRecalculating windows ... ");
  fflush(stdout);
  windowmaker(num_gaps, cnt, chrom_name, chrom_seq, length, 1);

  fprintf(stdout, "\n");
}


void windowmaker(int num_gaps, int chrom_id, char *chrom_name, char *chrom_seq, int length, int flag){
  int i;
  int j;
  int nchar;
  int s, e;

  int gc;
  int lw_cnt = 0; // large window (5Kb)
  int sw_cnt = 0; // small window (1Kb)
  int cw_cnt = 0; // copy number window (1kb non-ovp);

  /* if flag = 0; don't record it */

  if (!flag){
    for (i=0;i<num_gaps;i++){
      if (!strcmp(gaptable[i].chrom, chrom_name)){
	if (gaptable[i].start == gaptable[i].end)
	  chrom_seq[gaptable[i].start] = 'X';
	else{
	  for (j=gaptable[i].start; (j<gaptable[i].end && j<length); j++)
	    chrom_seq[j] = 'X';
	}
      }
    }
  }


  /* count cw_cnt */

  nchar = 0;
  s = 0 ; e = 0;
  gc = 0;

  for (i=0; i<length; i++){

    if (chrom_seq[i] != 'N' && chrom_seq[i] != 'X')
      nchar++;
    if (chrom_seq[i] == 'G' || chrom_seq[i] == 'C')
      gc++;

    if ((nchar == CW_SIZE || chrom_seq[i] == 'X') && nchar != 0){
      e = i;
      if (chrom_seq[i] == 'X'){
	nchar = 0;
	gc = 0;
      }
      else if (flag == 1 && nchar == CW_SIZE){
	/* save */
	chromosomes[chrom_id]->cw[cw_cnt].start = s;
	chromosomes[chrom_id]->cw[cw_cnt].end = e+1;
	chromosomes[chrom_id]->cw[cw_cnt].gc = (float)gc / (float)nchar;
	chromosomes[chrom_id]->cw[cw_cnt].depth = 0; // for readability/completeness
	cw_cnt++;
	gc = 0;
      }
      else if (flag == 0 && nchar == CW_SIZE){
	cw_cnt ++;
      }

      s = i + 1;
      nchar = 0;
    }

    if (chrom_seq[i] == 'X'){
      s = i + 1;
      nchar = 0;
      gc = 0;
    }

  }


  /* count sw_cnt */

  i = 0;
  nchar = 0;
  s = 0 ; e = 0;
  gc = 0;

  while (i < length){

    if (chrom_seq[i] != 'N' && chrom_seq[i] != 'X')
      nchar++;
    if (chrom_seq[i] == 'G' || chrom_seq[i] == 'C')
      gc++;

    if ((nchar == SW_SIZE || chrom_seq[i] == 'X') && nchar != 0){
      e = i;
      if (chrom_seq[i] == 'X'){
        nchar = 0;
	gc = 0;
      }
      else  if (flag == 1 && nchar == SW_SIZE){
	/* save */
	chromosomes[chrom_id]->sw[sw_cnt].start = s;
	chromosomes[chrom_id]->sw[sw_cnt].end = e+1;
	chromosomes[chrom_id]->sw[sw_cnt].gc = (float)gc / (float)nchar;
	chromosomes[chrom_id]->sw[sw_cnt].depth = 0; // for readability/completeness
	sw_cnt++;
	gc = 0;
      }
      else if (flag == 0 && nchar == SW_SIZE)
	sw_cnt++;

      s = s + SW_SLIDE;
      i = s - 1;
      nchar = 0;
    }
    if (chrom_seq[i]=='X'){
      s = i + 1;
      i = s - 1;
      gc = 0;
    }

    i++;
  }


  nchar = 0;
  i = 0;
  s = 0 ; e = 0;
  gc = 0;
  /* count lw_cnt */

  while (i<length){

    if (chrom_seq[i] != 'N' && chrom_seq[i] != 'X')
      nchar++;
    if (chrom_seq[i] == 'G' || chrom_seq[i] == 'C')
      gc++;

    if ((nchar == LW_SIZE || chrom_seq[i] == 'X') && nchar != 0){
      e = i;
      if (chrom_seq[i] == 'X'){
        nchar = 0;
	gc = 0;
      }
      else if (flag == 1 && nchar == LW_SIZE){
	/* save */
	chromosomes[chrom_id]->lw[lw_cnt].start = s;
	chromosomes[chrom_id]->lw[lw_cnt].end = e+1;
	chromosomes[chrom_id]->lw[lw_cnt].gc = (float)gc / (float)nchar;
	chromosomes[chrom_id]->lw[lw_cnt].depth = 0; // for readability/completeness
	lw_cnt++;
	gc = 0;
      }
      else if (flag == 0 && nchar == LW_SIZE)
	lw_cnt++;

      s = s + LW_SLIDE;
      i = s - 1;
      nchar = 0;
    }

    if (chrom_seq[i] == 'X') {
      s = i + 1;
      i = s - 1;
      gc = 0;
    }

    i++;
  }



  if (flag == 0){

    printf("\n\t\tLW: %d\tSW: %d\tCW: %d\n", lw_cnt, sw_cnt, cw_cnt);

    chromosomes[chrom_id]->lw_cnt = lw_cnt;
    chromosomes[chrom_id]->sw_cnt = sw_cnt;
    chromosomes[chrom_id]->cw_cnt = cw_cnt;

    chromosomes[chrom_id]->cw = (struct window *) getMem (sizeof(struct window) * cw_cnt);
    chromosomes[chrom_id]->lw = (struct window *) getMem (sizeof(struct window) * lw_cnt);
    chromosomes[chrom_id]->sw = (struct window *) getMem (sizeof(struct window) * sw_cnt);


  }


}

