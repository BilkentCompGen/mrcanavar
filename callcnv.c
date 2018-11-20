#include "callcnv.h"

//Comparison function for quick sort
static int q_compare(const void *p1, const void *p2){
  float a, b;

  a = *((float *)p1);
  b = *((float *)p2);

  if (a>b) 
    return 1;
  else
    return -1;
}


void call_cnv(char *depthFile, char *out_prefix,char* glist){

  int i, j;// k;

  float *gclookup;
  float *gclookup_x;

  char *fname;
  char logname[2 * MAX_STR];
  FILE *log;
  //float maxcutlw, maxcutsw;
  //float max;
  //float max_x;
  //float min;
  //float min_x;

  if (GENOME_CONF == NULL)
    print_error("Select genome configuration file (input) through the -conf parameter.\n");
  if (depthFile == NULL)
    print_error("Select read depth output file through the -depth parameter.\n");
  if (out_prefix == NULL)
    print_error("Select output file prefix through the -o parameter.\n");


  loadRefConfig(GENOME_CONF);

  readDepth(depthFile);


  sprintf(logname, "%s.log", out_prefix);
  log = my_fopen(logname, "w", 0);

  gclookup   = (float *) getMem(sizeof(float) * GC_BIN);
  gclookup_x = (float *) getMem(sizeof(float) * GC_BIN);


  /* trivial control removal */

  fprintf(stdout, "Control region cleanup...");
  fflush(stdout);


  /* add stdev calculation here */


  for (i=0; i<num_chrom; i++){
    for (j=0; j<chromosomes[i]->lw_cnt; j++){
      if (chromosomes[i]->lw[j].depth > (float) LW_MEAN * 2.0 || chromosomes[i]->lw[j].depth < (float) LW_MEAN / 10.0)
        chromosomes[i]->lw[j].isControl = 0;
      if (strstr(chromosomes[i]->name, "random") || strstr(chromosomes[i]->name, "Y") || strstr(chromosomes[i]->name, "hap") || strstr(chromosomes[i]->name, "Un"))
        chromosomes[i]->lw[j].isControl = 0;
    }
    for (j=0; j<chromosomes[i]->sw_cnt; j++){
      if (chromosomes[i]->sw[j].depth > (float) SW_MEAN * 2.0 || chromosomes[i]->sw[j].depth < (float) SW_MEAN / 10.0)
        chromosomes[i]->sw[j].isControl = 0;
      if (strstr(chromosomes[i]->name, "random") || strstr(chromosomes[i]->name, "Y") || strstr(chromosomes[i]->name, "hap") || strstr(chromosomes[i]->name, "Un"))
        chromosomes[i]->sw[j].isControl = 0;
    }

    for (j=0; j<chromosomes[i]->cw_cnt; j++){
      if (chromosomes[i]->cw[j].depth > (float) CW_MEAN * 2.0 || chromosomes[i]->cw[j].depth < (float) CW_MEAN / 10.0)
        chromosomes[i]->cw[j].isControl = 0;
      if (strstr(chromosomes[i]->name, "random") || strstr(chromosomes[i]->name, "Y") || strstr(chromosomes[i]->name, "hap") || strstr(chromosomes[i]->name, "Un"))
        chromosomes[i]->cw[j].isControl = 0;
    }
  }


  fprintf(stdout, "\n");

  norm_until_converges(CW, gclookup, gclookup_x);

  norm_until_converges(LW, gclookup, gclookup_x);

  norm_until_converges(SW, gclookup, gclookup_x);

  fprintf (stdout, "Writing normalized CW depth to: %s.cw_norm.bed.\n", out_prefix);

  fname = (char *) getMem(sizeof (char) * (strlen(out_prefix) + strlen(".copynumber.bed") + 2));

  sprintf (fname, "%s.cw_norm.bed", out_prefix);
  dump_text_windows(fname, CW);

  fprintf (stdout, "Writing normalized LW depth to: %s.lw_norm.bed.\n", out_prefix);

  sprintf (fname, "%s.lw_norm.bed", out_prefix);
  dump_text_windows(fname, LW);

  fprintf (stdout, "Writing normalized SW depth to: %s.sw_norm.bed.\n", out_prefix);

  sprintf (fname, "%s.sw_norm.bed", out_prefix);
  dump_text_windows(fname, SW);

  sprintf (fname, "%s.copynumber.bed", out_prefix);
  fprintf (stdout, "Writing copy numbers to: %s.copynumber.bed. \n", out_prefix);
  print_copy_numbers(fname);

  freeMem(fname, (strlen(out_prefix) + strlen(".copynumber.bed") + 2));



  fprintf (log, "\nmrCaNaVaR version %s\nLast update: %s\n\n", VERSION, LAST_UPDATE);
  fprintf (log, "Calculating library %s\n", out_prefix);
  fprintf (log, "GC correction mode: %s\n", MULTGC == 1 ? "MULTIPLICATIVE" : "ADDITIVE");

  fprintf (log, "\nAfter GC Correction:\n--------------------\n");
  fprintf (log, "Sample Gender: %s.\n", GENDER == MALE ? "Male" : "Female");
  fprintf (log, "CW Average Read Depth: %f, Standard Deviation: %f\n", CW_MEAN, CW_STD);


  if (GENDER == MALE)
    fprintf (log, "CW Average chrX Read Depth: %f, Standard Deviation: %f\n", CW_MEAN_X, CW_STD_X);

  fprintf (log, "LW Average Read Depth: %f, Standard Deviation: %f\n", LW_MEAN, LW_STD);
  if (GENDER == MALE)
    fprintf (log, "LW Average chrX Read Depth: %f, Standard Deviation: %f\n", LW_MEAN_X, LW_STD_X);

  fprintf (log, "SW Average Read Depth: %f, Standard Deviation: %f\n", SW_MEAN, SW_STD);
  if (GENDER == MALE)
    fprintf (log, "SW Average chrX Read Depth: %f, Standard Deviation: %f\n", SW_MEAN_X, SW_STD_X);

  fclose(log);

  calculate_dups(out_prefix);
  
  /* recalculate stats for deletion calls */
  
  fprintf(stderr, "Recalculating stats for deletion calls. LW: %f/%f SW: %f/%f\n", LW_MEAN, LW_STD, SW_MEAN, SW_STD);
  
  /*

  for (i=0; i<num_chrom; i++){
    if (isAutosome(NULL, -1, i)){
      maxcutlw = LW_MEAN * 2;
      maxcutsw = SW_MEAN * 2;
    }
    else{
      maxcutlw = LW_MEAN_X * 1.5;
      maxcutsw = SW_MEAN_X * 1.5;
    }
    for (j=0; j<chromosomes[i]->lw_cnt; j++){
      if (chromosomes[i]->lw[j].depth > maxcutlw){
	chromosomes[i]->lw[j].isControl = 0;
	k = j;
	while (k > 0 && chromosomes[i]->lw[k-1].end >= chromosomes[i]->lw[j].start){
	  chromosomes[i]->lw[--k].isControl = 0;
	}

	k = j;
	
	while (k < chromosomes[i]->lw_cnt-1 && chromosomes[i]->lw[k+1].start <= chromosomes[i]->lw[j].end){
	  chromosomes[i]->lw[++k].isControl = 0;
	}
      }
    }

    for (j=0; j<chromosomes[i]->sw_cnt; j++){
      if (chromosomes[i]->sw[j].depth > maxcutsw){
	chromosomes[i]->sw[j].isControl = 0;
	k = j;
	while (k > 0 && chromosomes[i]->sw[k-1].end >= chromosomes[i]->sw[j].start){
	  chromosomes[i]->sw[--k].isControl = 0;
	}

	k = j;
	
	while (k < chromosomes[i]->sw_cnt-1 && chromosomes[i]->sw[k+1].start <= chromosomes[i]->sw[j].end){
	  chromosomes[i]->sw[++k].isControl = 0;
	}
      }
    }
  }
  
  calc_stat(LW, gclookup, gclookup_x, 0, &max, &max_x, &min, &min_x);
  fprintf(stderr, "max %f maxx %f min %f minx %f\n", max, max_x, min, min_x);
  calc_stat(SW, gclookup, gclookup_x, 0, &max, &max_x, &min, &min_x);

  fprintf(stderr,"Now. LW: %f/%f SW: %f/%f\n", LW_MEAN, LW_STD, SW_MEAN, SW_STD);
  */
  calculate_dels(out_prefix);

  freeMem(gclookup, sizeof(float) * GC_BIN);
  freeMem(gclookup_x, sizeof(float) * GC_BIN);


  if(glist!=NULL)
    read_gene(glist, out_prefix);
}

void calculate_dups(char *out_prefix){

  printf("Calculating duplications....\n");
  calculate_intervals(out_prefix, DUPS);

}

void calculate_dels(char *out_prefix){

  printf("Calculating deletions....\n");
  // init dw here...
  calculate_intervals(out_prefix, DELS);
 
}

void calculate_intervals(char *out_prefix, enum CALLTYPE type){


  char* fname;

  int i, j,k,window;
  int cur_dup=0;

  int num_dups=0;

  float cutoff;
  float cutoff_x;

  if (type == DUPS){
    fname = (char *) getMem(sizeof (char) * (strlen(out_prefix) + strlen(".dups.bed") + 2));
    sprintf(fname, "%s.dups.bed", out_prefix);
    cutoff = LW_MEAN + 4*LW_STD;
    cutoff_x = LW_MEAN_X + 4*LW_STD_X;
  }
  else {
    fname = (char *) getMem(sizeof (char) * (strlen(out_prefix) + strlen(".dels.bed") + 2));
    sprintf(fname, "%s.dels.bed", out_prefix);
    cutoff = LW_MEAN - 4*LW_STD;
    cutoff_x = LW_MEAN_X - 4*LW_STD_X;
  }

  for (i = 0; i < num_chrom; i++){

    j=0;
    while (j < chromosomes[i]->lw_cnt){
      window=j+CONT_WINDOW;
      if ( window > chromosomes[i]->lw_cnt)
        window=chromosomes[i]->lw_cnt-1;

      //If the gene is Autosome go into this.
      if (isAutosome(chromosomes[i]->lw, j, i)){

        for( k=j; k< window; k++){	
	  if (type == DUPS){
	    if(chromosomes[i]->lw[k].depth > cutoff ){
	      num_dups++;
	    }
	  }
	  else {
	    if(chromosomes[i]->lw[k].depth < cutoff ){
	      num_dups++;
	    }	      
	  }
	}

	if(num_dups >= CUT_WINDOW){
	  chromosomes[i]->dw[cur_dup].start=chromosomes[i]->lw[j].start;
	  chromosomes[i]->dw[cur_dup].end=chromosomes[i]->lw[--k].end;
	  //Hit the first windows and now extend.
	  
	  if (type == DUPS){
	    while( chromosomes[i]->lw[++k].depth > cutoff && chromosomes[i]->dw[cur_dup].end > chromosomes[i]->lw[k].start ) {
	      chromosomes[i]->dw[cur_dup].end=chromosomes[i]->lw[k].end;
	    }
	  }
	  else{
	    while( chromosomes[i]->lw[++k].depth < cutoff && chromosomes[i]->dw[cur_dup].end > chromosomes[i]->lw[k].start ) {
	      chromosomes[i]->dw[cur_dup].end=chromosomes[i]->lw[k].end;
	    }
	  }
	  //Trim from start and end with sw's.
	  trim_sw(i,cur_dup,0, type);
	  
	  //Filter the ones that do not overlap a certain amount of CW's	
	  filter_dw(i,cur_dup,0, type);
	  cur_dup++;
	  
	}
	
	num_dups=0;
	j=k+1;
	
      }

      else{
	if(GENDER == FEMALE){
	  if (strstr(chromosomes[i]->name, "chrY"))
	    break;
	  if (strstr(chromosomes[i]->name, "Y_random"))
	    break;
	  if ( !strcmp(chromosomes[i]->name, "Y"))
            break;
        } 
	
        for( k=j; k< window; k++){	
          if(chromosomes[i]->lw[k].depth > (LW_MEAN_X + 4*LW_STD_X) ){
	    
            num_dups++;

          }
        }
	// We have a hit and we have to extend the window now.
        if(num_dups >= CUT_WINDOW){
          chromosomes[i]->dw[cur_dup].start=chromosomes[i]->lw[j].start;
          chromosomes[i]->dw[cur_dup].end=chromosomes[i]->lw[k-1].end;

	  if (type == DUPS){
	    while(  chromosomes[i]->lw[++k].depth > cutoff_x && chromosomes[i]->dw[cur_dup].end > chromosomes[i]->lw[k].start){
	      chromosomes[i]->dw[cur_dup].end=chromosomes[i]->lw[k].end;
	    }
	  }
	  else{
	    while(  chromosomes[i]->lw[++k].depth < cutoff_x && chromosomes[i]->dw[cur_dup].end > chromosomes[i]->lw[k].start){
	      chromosomes[i]->dw[cur_dup].end=chromosomes[i]->lw[k].end;
	    }
	  }
	  
          trim_sw(i, cur_dup, 1, type);

          filter_dw(i, cur_dup, 1, type);				
          cur_dup++;

        }

        num_dups=0;
        j=k;
      }
    }
    
    chromosomes[i]->dw_cnt=cur_dup;
    num_dups=0;
    cur_dup=0;
  }
  
  write_dw(fname,out_prefix);
  free(fname);
}

int is_cn_ok(int copy_num, int isAut, enum CALLTYPE type){

  int cnok;

  if (isAut){
    if (type==DUPS)
      cnok = 2.5;
    else
      cnok = 1.5;
  }
  else{
    if (type==DUPS)
      cnok = 1.5;
    else
      cnok = 0.5;
  }

  if (type == DUPS){
    if(copy_num < cnok){
      return 0;
    }
    else{
      return 1;
    }
  }
  else{
    if(copy_num > cnok){
      return 0;
    }
    else{
      return 1;
    }
  }
  
}

void filter_dw(int i, int cur_dup, int isAut, enum CALLTYPE type){
  int j=0;
  float copy_num;
  float isOK=0.0;	// Number of CW's that pass the filter.
  float isNOK=0.0;      // Number of CW's that fail to pass the filter.

  //Scan all CW's for overlapping intervals and count them.
  for(j=0; j< chromosomes[i]->cw_cnt; j++){
    if(isAut==0){
      //CN calculation
      copy_num = (chromosomes[i]->cw[j].depth / CW_MEAN) * 2;
      if(chromosomes[i]->dw[cur_dup].start<=chromosomes[i]->cw[j].start && chromosomes[i]->dw[cur_dup].end >= chromosomes[i]->cw[j].end){

	if (type == DUPS){
	  if(copy_num < 2.5){
	    isNOK++;
	  }
	  else{
	    isOK++;
	  }
	}
	else{
	  if(copy_num > 1.5){
	    isNOK++;
	  }
	  else{
	    isOK++;
	  }
	}
	continue;
      }

      if(chromosomes[i]->dw[cur_dup].end > chromosomes[i]->cw[j].start && chromosomes[i]->dw[cur_dup].end < chromosomes[i]->cw[j].end ){
	if((chromosomes[i]->dw[cur_dup].end - chromosomes[i]->cw[j].start) / (chromosomes[i]->cw[j].end -chromosomes[i]->cw[j].start) >= FILTER && is_cn_ok(copy_num, 1, type)){
	  isOK++;
	}else{
	  isNOK++;
	}
	continue;
      }

      if(chromosomes[i]->dw[cur_dup].start > chromosomes[i]->cw[j].start && chromosomes[i]->dw[cur_dup].start < chromosomes[i]->cw[j].end){
	if((chromosomes[i]->cw[j].end - chromosomes[i]->dw[cur_dup].start) / (chromosomes[i]->cw[j].end -chromosomes[i]->cw[j].start) >= FILTER && is_cn_ok(copy_num, 1, type)){
	  isOK++;
	}else{
	  isNOK++;
	}
	continue;
      }
    }
    else{       //If not autosomal
      //CN calculation
      copy_num = chromosomes[i]->cw[j].depth / CW_MEAN_X;
      if(chromosomes[i]->dw[cur_dup].start<=chromosomes[i]->cw[j].start && chromosomes[i]->dw[cur_dup].end >= chromosomes[i]->cw[j].end){
	//	if(copy_num < 1.5){
	if (!is_cn_ok(copy_num, 0, type)){ 
	  isNOK++;
	}else{
	  isOK++;
	}
	continue;
      }

      if(chromosomes[i]->dw[cur_dup].end > chromosomes[i]->cw[j].start && chromosomes[i]->dw[cur_dup].end < chromosomes[i]->cw[j].end ){
	if((chromosomes[i]->dw[cur_dup].end - chromosomes[i]->cw[j].start) / (chromosomes[i]->cw[j].end -chromosomes[i]->cw[j].start) >= FILTER && is_cn_ok(copy_num, 0, type)){
	  isOK++;
	}else{
	  isNOK++;
	}
	continue;
      }

      if(chromosomes[i]->dw[cur_dup].start > chromosomes[i]->cw[j].start && chromosomes[i]->dw[cur_dup].start < chromosomes[i]->cw[j].end){
	if((chromosomes[i]->cw[j].end - chromosomes[i]->dw[cur_dup].start) / (chromosomes[i]->cw[j].end -chromosomes[i]->cw[j].start) >= FILTER && is_cn_ok(copy_num, 0, type)){
	  isOK++;
	}else{
	  isNOK++;
	}
	continue;
      }
    }
  }
  //Filter is applied in here.
  if((isOK/isNOK) < (1.0/3.0)){
    chromosomes[i]->dw[cur_dup].start=-1;
    chromosomes[i]->dw[cur_dup].end=-1;
  }
 
}

void trim_sw(int i, int cur_dup, int isAut, enum CALLTYPE type){
  int j=0;

  float cutoff, cutoff_x;

  if (type == DUPS){
    cutoff = SW_MEAN + 4*SW_STD;
    cutoff_x = SW_MEAN_X + 4*SW_STD_X;
  }
  else {
    cutoff = SW_MEAN - 4*SW_STD;
    cutoff_x = SW_MEAN_X - 4*SW_STD_X;
  }


  for(j=0; j< chromosomes[i]->sw_cnt; j++){
    if(isAut==0){
      if(chromosomes[i]->dw[cur_dup].start==chromosomes[i]->sw[j].start && chromosomes[i]->dw[cur_dup].end > chromosomes[i]->sw[j].end){
	if (type == DUPS){
	  if(chromosomes[i]->sw[j].depth < cutoff){
	    chromosomes[i]->dw[cur_dup].start=chromosomes[i]->sw[j-1].end;
	    continue;			
	  }
	}
	else{
	  if(chromosomes[i]->sw[j].depth > cutoff){
	    chromosomes[i]->dw[cur_dup].start=chromosomes[i]->sw[j-1].end;
	    continue;			
	  }
	}
      }

      if(chromosomes[i]->dw[cur_dup].end == chromosomes[i]->sw[j].end ){
	if (type == DUPS){
	  if(chromosomes[i]->sw[j].depth < cutoff){
	    chromosomes[i]->dw[cur_dup].end=chromosomes[i]->sw[j+1].end;
	    continue;			
	  }
	}
	else{
	  if(chromosomes[i]->sw[j].depth > cutoff){
	    chromosomes[i]->dw[cur_dup].end=chromosomes[i]->sw[j+1].end;
	    continue;			
	  }
	}
      }
    }
    else{
      if(chromosomes[i]->dw[cur_dup].start==chromosomes[i]->sw[j].start ){
	if (type == DUPS){
	  if(chromosomes[i]->sw[j].depth < cutoff_x){
	    chromosomes[i]->dw[cur_dup].start=chromosomes[i]->sw[j-1].end;
	    continue;			
	  }
	}
	else{
	  if(chromosomes[i]->sw[j].depth > cutoff_x){
	    chromosomes[i]->dw[cur_dup].start=chromosomes[i]->sw[j-1].end;
	    continue;			
	  }
	}
      }

      if(chromosomes[i]->dw[cur_dup].end == chromosomes[i]->sw[j].end ){
	if (type == DUPS){
	  if(chromosomes[i]->sw[j].depth < cutoff_x){
	    chromosomes[i]->dw[cur_dup].end=chromosomes[i]->sw[j+1].start;
	    continue;			
	  }
	}
	else{
	  if(chromosomes[i]->sw[j].depth > cutoff_x){
	    chromosomes[i]->dw[cur_dup].end=chromosomes[i]->sw[j+1].start;
	    continue;			
	  }
	}
      }
    }
  }
}

void write_dw(char* fname,char *out_prefix){

  FILE *txtDepth;
  int i,j,k;
  int collapse;

  txtDepth = my_fopen(fname, "w", 0);

  fprintf(txtDepth, "#%s\t%s\t%s\t\n\n", "CHROM", "START", "END");

  printf("Writing to : %s\n",fname);

  for (i = 0; i < num_chrom; i++){
    if (GENDER == FEMALE && (strstr(chromosomes[i]->name, "Y")))
      continue;

    j = 0;		
    while(j < chromosomes[i]->dw_cnt){
      collapse = j;	
      while ((chromosomes[i]->dw[collapse].end > chromosomes[i]->dw[++collapse].start) && (collapse < chromosomes[i]->dw_cnt))
        ;
      collapse--;
      chromosomes[i]->dw[j].end=chromosomes[i]->dw[collapse].end;
      for(k=j+1; k<= collapse; k++){
        chromosomes[i]->dw[k].start=-1;
        chromosomes[i]->dw[k].end=-1;
      }
      if (chromosomes[i]->dw[j].start != -1 && chromosomes[i]->dw[j].end != -1 && chromosomes[i]->dw[j].end - chromosomes[i]->dw[j].start >= MIN_DUP){

        fprintf(txtDepth, "%s\t%d\t%d\n", chromosomes[i]->name, chromosomes[i]->dw[j].start, chromosomes[i]->dw[j].end);
        fflush(txtDepth);
      }
      j = collapse+1;				
    }

  }
  fclose(txtDepth);

}

void readDepth(char *depthFile){

  FILE *binDepth;
  int i, j;
  int retVal;

  double lw_total;
  double sw_total;
  double cw_total;

  int lw_cnt;
  int sw_cnt;
  int cw_cnt;

  int isMagicNum;


  binDepth = my_fopen(depthFile, "r", 0);

  retVal = fread(&isMagicNum, sizeof(isMagicNum), 1, binDepth);


  if (isMagicNum == magicNumDepth)
    ISNORMALIZED = 1;
  else if (isMagicNum == magicNum)
    fprintf(stdout, "Read depth file is not yet normalized.\n");
  else
    print_error("Read depth file seems to be invalid or corrupt.\n");



  lw_total = 0.0;
  sw_total = 0.0;
  cw_total = 0.0;

  lw_cnt   = 0;
  sw_cnt   = 0;
  cw_cnt   = 0;

  float lwdepth;
  float swdepth;
  float cwdepth;

  /* read LW */

  for (i = 0; i < num_chrom; i++){
    for (j = 0; j < chromosomes[i]->lw_cnt; j++){
      retVal = fread(&(lwdepth), sizeof(chromosomes[i]->lw[j].depth), 1, binDepth);
      chromosomes[i]->lw[j].depth += lwdepth;
      lw_total += lwdepth;
      //      lw_total += chromosomes[i]->lw[j].depth;
      chromosomes[i]->lw[j].isControl = 1;
      lw_cnt++;
    }
  }

  /* read SW */

  for (i = 0; i < num_chrom; i++){
    for (j = 0; j < chromosomes[i]->sw_cnt; j++){
      retVal = fread(&(swdepth), sizeof(chromosomes[i]->sw[j].depth), 1, binDepth);
      chromosomes[i]->sw[j].depth += swdepth;
      sw_total += swdepth;
      //sw_total += chromosomes[i]->sw[j].depth;
      chromosomes[i]->sw[j].isControl = 1;
      sw_cnt++;
    }
  }
  /* read CW */

  for (i = 0; i < num_chrom; i++){
    for (j = 0; j < chromosomes[i]->cw_cnt; j++){
      retVal = fread(&(cwdepth), sizeof(chromosomes[i]->cw[j].depth), 1, binDepth);
      chromosomes[i]->cw[j].depth += cwdepth;
      //     cw_total += chromosomes[i]->cw[j].depth;
      chromosomes[i]->cw[j].isControl = 1;
      cw_total +=cwdepth; 
      cw_cnt++;
    }
  }

  LW_MEAN = lw_total / lw_cnt;
  SW_MEAN = sw_total / sw_cnt;
  CW_MEAN = cw_total / cw_cnt;

  fprintf(stdout, "[OK] depth file %s is loaded.\n", depthFile);
  if (VERBOSE){
    fprintf(stdout, "LW_MEAN: %f\tSW_MEAN: %f\tCW_MEAN:%f\n",  LW_MEAN, SW_MEAN, CW_MEAN);
    if (retVal == 0)
      fprintf(stderr, "There was potentially an error in fread.\n");
  }
  fclose(binDepth);
}


void dump_text_windows(char *fname, enum WINDOWTYPE wt){

  FILE *txtDepth;
  int i, j;

  txtDepth = my_fopen(fname, "w", 0);

  fprintf(txtDepth, "#%s\t%s\t%s\t%s\t%s\t%s\n\n", "CHROM", "START", "END", "GC\%", "READ_DEPTH", "IS_CONTROL");

  switch (wt){
  case LW:
    for (i = 0; i < num_chrom; i++){
      for (j = 0; j < chromosomes[i]->lw_cnt; j++)
        if (chromosomes[i]->lw[j].isControl == 1)
          fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\tY\n", chromosomes[i]->name, chromosomes[i]->lw[j].start, chromosomes[i]->lw[j].end, chromosomes[i]->lw[j].gc, chromosomes[i]->lw[j].depth);
        else
          fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\tN\n", chromosomes[i]->name, chromosomes[i]->lw[j].start, chromosomes[i]->lw[j].end, chromosomes[i]->lw[j].gc, chromosomes[i]->lw[j].depth);
    }
    break;

  case SW:
    for (i = 0; i < num_chrom; i++){
      for (j = 0; j < chromosomes[i]->sw_cnt; j++)
	if (chromosomes[i]->sw[j].isControl == 1)
	  fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\tY\n", chromosomes[i]->name, chromosomes[i]->sw[j].start, chromosomes[i]->sw[j].end, chromosomes[i]->sw[j].gc, chromosomes[i]->sw[j].depth);
	else
	  fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\tN\n", chromosomes[i]->name, chromosomes[i]->sw[j].start, chromosomes[i]->sw[j].end, chromosomes[i]->sw[j].gc, chromosomes[i]->sw[j].depth);
    }
    break;

  case CW:
    for (i = 0; i < num_chrom; i++){
      for (j = 0; j < chromosomes[i]->cw_cnt; j++)
	if (chromosomes[i]->cw[j].isControl == 1)
	  fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\tY\n", chromosomes[i]->name, chromosomes[i]->cw[j].start, chromosomes[i]->cw[j].end, chromosomes[i]->cw[j].gc, chromosomes[i]->cw[j].depth);
	else
	  fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\tN\n", chromosomes[i]->name, chromosomes[i]->cw[j].start, chromosomes[i]->cw[j].end, chromosomes[i]->cw[j].gc, chromosomes[i]->cw[j].depth);
    }
    break;
  }

  fclose(txtDepth);

}

void  print_copy_numbers(char *fname){

  FILE *txtDepth;
  int i, j;
  float copy_num;

  txtDepth = my_fopen(fname, "w", 0);
	
  fprintf(txtDepth, "#%s\t%s\t%s\t%s\t%s\n\n", "CHROM", "START", "END", "GC\%", "COPYNUMBER");

  for (i = 0; i < num_chrom; i++){
    for (j = 0; j < chromosomes[i]->cw_cnt; j++){
      if (isAutosome(chromosomes[i]->cw, j, i))
	copy_num = (chromosomes[i]->cw[j].depth / CW_MEAN) * 2;
      else
	copy_num = chromosomes[i]->cw[j].depth / CW_MEAN_X;

      if (GENDER == FEMALE && (strstr(chromosomes[i]->name, "Y")))
	continue;

      fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\n", chromosomes[i]->name, chromosomes[i]->cw[j].start, chromosomes[i]->cw[j].end, chromosomes[i]->cw[j].gc, copy_num);
    }
  }


  fclose(txtDepth);


}

void conc_depth(char **c_depths, int fileCount, char *depthFile){
  int i;

  if (GENOME_CONF == NULL)
    print_error("Select genome configuration file (input) through the -conf parameter.\n");
  if (c_depths == NULL) 
    print_error("Select read depth input files through the -concdepth parameter.\n");  
  if (depthFile == NULL)
    print_error("Select read depth output file through the -depth parameter.\n");


  loadRefConfig(GENOME_CONF);

  for(i = 0; i < fileCount; i++)
    readDepth(c_depths[i]);

  saveDepth(depthFile, 1);
}
void read_gene(char *glist, char *out_prefix){
     
  FILE *g_out;
  FILE* geneList;
  char* out_fname;
  int i,j,k,start,end,count,cur;
  float copy_num,total;
  char rest_of_list[1000];
  float copy_n[10000];
  char  chr[100];

  char *ret = 0;
     
  out_fname = (char *) getMem(sizeof (char) * (strlen(glist) + strlen(".out.bed") + 2));

  geneList=fopen(glist, "r");
     
  sprintf (out_fname, "%s.genes.bed", out_prefix);

  g_out= my_fopen(out_fname, "w", 0);
	
  fprintf(stdout, "Reading gene list: %s ....\n",glist);
  fprintf(stdout, "Writing to: %s (This may take a while)\n",out_fname);
     
  fprintf(g_out, "#%s\t%s\t%s\t%s\t%s\t%s\n\n", "CHROM", "START", "END", "BP_LENGTH" ,"MEDIAN", "AVERAGE");

  for (i = 0; i < num_chrom; i++){
        
    while (fscanf(geneList, "%s\t%d\t%d", chr, &start, &end) > 0){
      //Gets the portion of the genelist file after 'END'
      ret = fgets(rest_of_list, 1000, geneList);
      if (ret == NULL){
	fprintf(stderr, "Gene file read error.\n");
	exit(-1);
      }
      rest_of_list[strlen(rest_of_list)-1] = 0;
	 
      if (strcmp(chr, chromosomes[i]->name)) 
	continue;

      if (start > end){
	fprintf(stderr, "Start coordinate (%d) is bigger then end (%d). Exiting.\n", start, end); 
	exit(0); 
      }
      count = end-start+1;
      //Find the interleaving CW's of gene list
      for (j = 0; j < chromosomes[i]->cw_cnt; j++){
	if (isAutosome(chromosomes[i]->cw, j, i)){
	  copy_num = (chromosomes[i]->cw[j].depth / CW_MEAN) * 2;
	  if(chromosomes[i]->cw[j].start >= start && chromosomes[i]->cw[j].end <= end){
	    copy_n[cur]=copy_num;
	    cur++;
	  }
	  if(chromosomes[i]->cw[j].start < start && chromosomes[i]->cw[j].end > start){
	    copy_n[cur]=copy_num;
	    cur++;
	  }
	  if(chromosomes[i]->cw[j].start < end && chromosomes[i]->cw[j].end > end){
	    copy_n[cur]=copy_num;
	    cur++;
	  }
	}
	else{
	  copy_num = chromosomes[i]->cw[j].depth / CW_MEAN_X;
	  if(chromosomes[i]->cw[j].start >= start && chromosomes[i]->cw[j].end <= end){
	    copy_n[cur]=copy_num;
	    cur++;
	  }
	  if(chromosomes[i]->cw[j].start < start && chromosomes[i]->cw[j].end > start){
	    copy_n[cur]=copy_num;
	    cur++;
	  }
	  if(chromosomes[i]->cw[j].start < end && chromosomes[i]->cw[j].end > end){
	    copy_n[cur]=copy_num;
	    cur++;
	  }
	}
      }
      total=0.0;
      for(k=0;k < cur;k++){
	total+=copy_n[k];
      }
	 
      qsort(copy_n, cur, sizeof (float), q_compare);
      //Print the output of gene list
      if (total != 0)
	fprintf(g_out, "%s\t%d\t%d\t%d\t%s\t%f\t%f\n", chr, start, end, count, rest_of_list, copy_n[cur/2], (total/cur) );
      else
	fprintf(g_out, "%s\t%d\t%d\t%d\t%s\t%f\t%f\n", chr, start, end, count, rest_of_list, 0.0, 0.0);
	  
      cur=0;
      count=0;
    }
	
    rewind(geneList);
  }

}

