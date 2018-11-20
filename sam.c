#include "sam.h"


void read_mapfiles(char *sams, char *depthFile, char gzSAM, int dir_or_list){

  DIR *dp;

  int fileCnt;
  int totalFile;
  int n_iter;
  int cntr;
  int i,j;

  struct dirent *ep;
  float *gclookup;
  float *gclookup_x;

  long lw_total;
  long sw_total;
  long cw_total;

  int lw_cnt;
  int sw_cnt;
  int cw_cnt;
  char *token;
  char *safe_sams = NULL;
  char **filenames;
  
  lw_total = 0.0;
  sw_total = 0.0;
  cw_total = 0.0;

  lw_cnt   = 0;
  sw_cnt   = 0;
  cw_cnt   = 0;

  gclookup = NULL;
  gclookup_x = NULL;


  if (GENOME_CONF == NULL)
    print_error("Select genome configuration file (input) through the -conf parameter.\n");
  if (sams == NULL)
    print_error("Select input directory that contains the SAM files through the -samdir parameter.\n");
  if (depthFile == NULL)
    print_error("Select read depth output file through the -depth parameter.\n");


  loadRefConfig(GENOME_CONF);
  set_str(&safe_sams, sams);
  
  fileCnt = 0;
  totalFile = 0;
  n_iter = 0;
  cntr = 0;

  if (dir_or_list == SAMDIR){
    fprintf(stdout, "Scanning the SAM directory: %s\n", sams);
    
    dp = opendir(sams);
    
    
    if (dp == NULL)
      print_error("SAM directory not found.\n");
   
    
    while((ep=readdir(dp))){
      if (ep->d_name[0] == '.')
	continue;
      if (ep->d_type == DT_DIR)
	continue;
      
      if (CHECKSAM && !endswith(ep->d_name, ".sam") && !endswith(ep->d_name, ".sam.gz"))
	continue;
      
      /*
	if (!strstr(ep->d_name, "sam"))
	continue;
      */
      
      i = strlen(ep->d_name)-1;
      
      if ((ep->d_name[i] == 'z' && ep->d_name[i-1] == 'g' && ep->d_name[i-2] == '.') && gzSAM == 0){
	print_error("File name ends with .gz yet --gz option is not selected. Are you sure?\nIf the files are indeed uncompressed, rename the files.\n");
      }
      totalFile++;
    }
    
    rewinddir(dp);
  }

  else {
    token = strtok(safe_sams, ",");
    while (token != NULL){
      if (endswith(token, ".gz") && gzSAM == 0)
	print_error("File name ends with .gz yet --gz option is not selected. Are you sure?\nIf the files are indeed uncompressed, rename the files.\n");
      totalFile++;
      token = strtok (NULL, ",");
    }
    filenames = (char **) getMem(sizeof(char *) * totalFile);
    memset(filenames, 0, sizeof(char *) * totalFile);
    i = 0;
    free(safe_sams);
    safe_sams = NULL;
    set_str(&safe_sams, sams);

    token = strtok(safe_sams, ",");
    while (token != NULL){
      if (endswith(token, ".gz") && gzSAM == 0)
	print_error("File name ends with .gz yet --gz option is not selected. Are you sure?\nIf the files are indeed uncompressed, rename the files.\n");
      set_str(&(filenames[i++]), token);
      token = strtok (NULL, ",");
    }
  }

  
  
  n_iter = totalFile; // Defult case for SAM's in samdir

  if(n_sam_e != 0 && n_sam_e <= totalFile )
    n_iter = (n_sam_e-n_sam_b) + 1;

  if (dir_or_list == SAMDIR){
    while((ep=readdir(dp)) && n_iter > 0){
      
      if (ep->d_name[0] == '.')
	continue;
      if (ep->d_type == DT_DIR)
	continue;
      
      if (CHECKSAM && !endswith(ep->d_name, ".sam") && !endswith(ep->d_name, ".sam.gz"))
	continue;
      
      if ((n_sam_b-1) > cntr){
	++cntr;
	continue;
      }
      
      
      fprintf(stdout, "\r                                                      \rLoading file %d of total %d: %s...", (fileCnt+1), totalFile, ep->d_name);
      fflush(stdout);
      
      readSAM(sams, dir_or_list, ep->d_name, gzSAM, &cw_total, &sw_total, &lw_total);
      
      if (VERBOSE) fprintf(stderr, "\nREADSAM: %ld\t%ld\t%ld\n", lw_total, sw_total, cw_total);
      
      fileCnt++;
      n_iter--;
    }
    
    closedir(dp);
  }
  
  else {        
    // process list
    /*
    token = strtok(safe_sams, ",");
    while (token != NULL){
      fprintf(stdout, "\r                                                      \rLoading file %d of total %d: %s...", (fileCnt+1), totalFile, token);
      fflush(stdout);
      if (CHECKSAM && !endswith(token, ".sam") && !endswith(token, ".sam.gz"))
	continue;
      readSAM(sams, dir_or_list, token, gzSAM, &cw_total, &sw_total, &lw_total);
      token = strtok (NULL, ",");
      fileCnt++;
      }*/
    for (i=0; i<totalFile; i++){
      fprintf(stdout, "\r                                                      \rLoading file %d of total %d: %s...", (fileCnt+1), totalFile, filenames[i]);
      fflush(stdout);
      if (CHECKSAM && !endswith(filenames[i], ".sam") && !endswith(filenames[i], ".sam.gz"))
	continue;
      readSAM(sams, dir_or_list, filenames[i], gzSAM, &cw_total, &sw_total, &lw_total);
      fileCnt++;
    }

    for (i=0; i<totalFile; i++)
      free(filenames[i]);
    free(filenames);
    
  }

  free(safe_sams);
  
  if (fileCnt == 0)
    print_error("SAM directory or list does not contain any files with extensions \".sam\" or \".sam.gz\".\nIf you do have the correct files, please rename them to have either \".sam\" or \".sam.gz\" extensions.\n");
  else
    fprintf(stdout, "\n\n%d file%s loaded.\n", fileCnt, fileCnt>1 ? "s" : "");
  

  //  if(RUNMODE != CONC){
		 
  
    /* trivial control removal */
    fprintf(stdout, "Control region cleanup...");
    fflush(stdout);


    /* add stdev calculation here */

    for (i=0; i<num_chrom; i++){
      lw_cnt+=chromosomes[i]->lw_cnt;
      sw_cnt+=chromosomes[i]->sw_cnt;
      cw_cnt+=chromosomes[i]->cw_cnt;
    }
    
    LW_MEAN = (double) lw_total / lw_cnt;
    SW_MEAN = (double) sw_total / sw_cnt;
    CW_MEAN = (double) cw_total / cw_cnt;

    if (VERBOSE){
      fprintf(stdout, "\nLW_TOTAL: %ld\tSW_TOTAL: %ld\tCW_TOTAL:%ld\n",  lw_total, sw_total, cw_total);
      fprintf(stdout, "\nLW_CNT: %d\tSW_CNT: %d\tCW_CNT:%d\n",  lw_cnt, sw_cnt, cw_cnt);
      fprintf(stdout, "\nLW_MEAN: %f\tSW_MEAN: %f\tCW_MEAN:%f\n",  LW_MEAN, SW_MEAN, CW_MEAN);
    }



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
	else if (strstr(chromosomes[i]->name, "random") || strstr(chromosomes[i]->name, "Y") || strstr(chromosomes[i]->name, "hap") || strstr(chromosomes[i]->name, "Un"))
	  chromosomes[i]->cw[j].isControl = 0;
      }
    }


    fprintf(stdout, "[OK].\n");

    gclookup   = (float *) getMem(sizeof(float) * GC_BIN);
    gclookup_x = (float *) getMem(sizeof(float) * GC_BIN);

    fprintf(stdout, "Writing non-normalized files.\n");

    saveDepth(depthFile, 0);

    norm_until_converges(CW, gclookup, gclookup_x);

    norm_until_converges(LW, gclookup, gclookup_x);

    norm_until_converges(SW, gclookup, gclookup_x);

    //  }

  saveDepth(depthFile, 1);

  if (gclookup != NULL)
    freeMem(gclookup, sizeof(float) * GC_BIN);
  if (gclookup_x != NULL)
    freeMem(gclookup_x, sizeof(float) * GC_BIN);
  
}



void readSAM(char *indirSAM, int dir_or_list, char *fname, char gzSAM, long *cw_total, long *sw_total, long *lw_total){
  FILE *sam;

  char *samfile;

  char samLine[4 * MAX_STR];
  char prevChrom[MAX_STR];
  char chrom[MAX_STR];
  char *token;

  int pos;
  int chrom_id;
  int prev_chrom_id = -1;

  samfile = (char *) getMem (sizeof(char) * (strlen(indirSAM) + strlen(fname) + 2));

  if (dir_or_list == SAMDIR){
    sprintf(samfile, "%s/%s", indirSAM, fname);
    sam = my_fopen(samfile, "r", gzSAM);
  }
  else
    sam = my_fopen(fname, "r", gzSAM);

  prevChrom[0] = 0;

  while (1){

    /* read entire line from the SAM file */
    if (gzSAM){
      gzgets((gzFile) sam, samLine, (4 * MAX_STR));
      if (gzeof((gzFile) sam))
	break;
    }
    else{
      token = fgets(samLine, (4 * MAX_STR), sam);
      if (feof(sam))
	break;
    }

    /* ignore if SAM header  */

    if (samLine[0]=='@')
      continue;

    /* parse chrom */
    token = NULL;

    token = strtok(samLine, "\t"); // read name
    if (token == NULL){
      fprintf(stderr, "\n\n** WARNING ** The SAM file %s seems to be corrupt. Skipping it.\n\n", fname);
      if (gzSAM)
	gzclose((gzFile)sam);
      else
	fclose(sam);
      freeMem(samfile, sizeof(char) * (strlen(indirSAM) + strlen(fname) + 2));
      return;
    }
    token = strtok(NULL,    "\t"); // flag
    if (token == NULL){
      fprintf(stderr, "\n\n** WARNING ** The SAM file %s seems to be corrupt. Skipping it.\n\n", fname);
      if (gzSAM)
	gzclose((gzFile)sam);
      else
	fclose(sam);
      freeMem(samfile, sizeof(char) * (strlen(indirSAM) + strlen(fname) + 2));
      return;
    }
    token = strtok(NULL,    "\t"); //chrom
    if (token == NULL){
      fprintf(stderr, "\n\n** WARNING ** The SAM file %s seems to be corrupt. Skipping it.\n\n", fname);
      if (gzSAM)
	gzclose((gzFile)sam);
      else
	fclose(sam);
      freeMem(samfile, sizeof(char) * (strlen(indirSAM) + strlen(fname) + 2));
      return;
    }

    strcpy(chrom, token);
    trimspace(chrom);

    token = strtok(NULL,    "\t"); // pos
    if (token == NULL){
      fprintf(stderr, "\n\n** WARNING ** The SAM file %s seems to be corrupt. Skipping it.\n\n", fname);
      if (gzSAM)
	gzclose((gzFile)sam);
      else
	fclose(sam);
      freeMem(samfile, sizeof(char) * (strlen(indirSAM) + strlen(fname) + 2));
      return;
    }

    pos = atoi(token) - 1;  // SAM file is 1-based; our format is 0-based

    if (pos < 0)
      continue;

    /*
      debug if needed
      fprintf(stdout, "%s\t%d\n", chrom, pos);
    */



    chrom_id = insert_read_lw(chrom, pos, prevChrom, prev_chrom_id, lw_total);

    /*
      if (VERBOSE)
      fprintf(stdout, "[READSAM]\t%s\t%d\tchrom_id: %d\n", chrom, pos, chrom_id);
    */

    if (chrom_id == -1) // this chromosome is not in the config file
      continue;

    chrom_id = insert_read_sw(chrom, pos, prevChrom, chrom_id, sw_total);

    if (chrom_id == -1) // this chromosome is not in the config file; it shouldn't come to this
      continue;

    chrom_id = insert_read_cw(chrom, pos, prevChrom, chrom_id, cw_total);



    if (chrom_id == -1) // this chromosome is not in the config file; it shouldn't come to this
      continue;

    strcpy(prevChrom, chrom);
    prev_chrom_id = chrom_id;


  }

  if (gzSAM)
    gzclose((gzFile)sam);
  else
    fclose(sam);

  freeMem(samfile, sizeof(char) * (strlen(indirSAM) + strlen(fname) + 2));
}





int insert_read_lw(char *chrom, int pos, char *prevChrom, int prev_chrom_id, long *lw_total){
  int chrom_id;
  int window_id;

  int flank;

  chrom_id = findChrom(chrom, prevChrom, prev_chrom_id);

  if (chrom_id != -1){

    if (pos >= chromosomes[chrom_id]->length){
      return chrom_id;
    }

    window_id = windowSearch(chromosomes[chrom_id]->lw, chromosomes[chrom_id]->lw_cnt, pos);


    if (window_id != -1){

      chromosomes[chrom_id]->lw[window_id].depth += 1;
      (*lw_total)++;

      /* iterate left */
      flank = window_id - 1;

      while (flank >= 0 && (pos >= chromosomes[chrom_id]->lw[flank].start && pos <= chromosomes[chrom_id]->lw[flank].end)){
	chromosomes[chrom_id]->lw[flank--].depth += 1;
	(*lw_total)++;
      }

      /* iterate right */
      flank = window_id + 1;

      while (flank < chromosomes[chrom_id]->lw_cnt && (pos >= chromosomes[chrom_id]->lw[flank].start && pos <= chromosomes[chrom_id]->lw[flank].end)){
	chromosomes[chrom_id]->lw[flank++].depth += 1;
	(*lw_total)++;
      }

    }

  }

  return chrom_id;
}



int insert_read_sw(char *chrom, int pos, char *prevChrom, int prev_chrom_id, long *sw_total){
  int chrom_id;
  int window_id;

  int flank;

  chrom_id = findChrom(chrom, prevChrom, prev_chrom_id);

  if (chrom_id != -1){

    if (pos >= chromosomes[chrom_id]->length)
      return chrom_id;

    window_id = windowSearch(chromosomes[chrom_id]->sw, chromosomes[chrom_id]->sw_cnt, pos);

    if (window_id != -1){

      chromosomes[chrom_id]->sw[window_id].depth += 1;
      (*sw_total)++;
      /* iterate left */
      flank = window_id - 1;

      while (flank >= 0 && (pos >= chromosomes[chrom_id]->sw[flank].start && pos <= chromosomes[chrom_id]->sw[flank].end)){
	chromosomes[chrom_id]->sw[flank--].depth += 1;
	(*sw_total)++;
      }

      /* iterate right */
      flank = window_id + 1;

      while (flank < chromosomes[chrom_id]->sw_cnt && (pos >= chromosomes[chrom_id]->sw[flank].start && pos <= chromosomes[chrom_id]->sw[flank].end)){
	chromosomes[chrom_id]->sw[flank++].depth += 1;
	(*sw_total)++;
      }


    }

  }

  return chrom_id;
}


int insert_read_cw(char *chrom, int pos, char *prevChrom, int prev_chrom_id, long *cw_total){
  int chrom_id;
  int window_id;

  chrom_id = findChrom(chrom, prevChrom, prev_chrom_id);

  if (chrom_id != -1){

    if (pos >= chromosomes[chrom_id]->length)
      return chrom_id;

    window_id = windowSearch(chromosomes[chrom_id]->cw, chromosomes[chrom_id]->cw_cnt, pos);

    if (window_id != -1){

      chromosomes[chrom_id]->cw[window_id].depth += 1;
      (*cw_total)++;
    }

  }
  return chrom_id;
}


int findChrom(char *chrom, char *prevChrom, int prev_chrom_id){
  int chrom_id;


  if (!strcmp(chrom, prevChrom)){
    chrom_id = prev_chrom_id;
  }

  else if (prev_chrom_id == -1 && !strcmp(chrom, chromosomes[0]->name)){ // the first entry in the file
    chrom_id = 0;
  }
  else if (prev_chrom_id != -1 && prev_chrom_id < num_chrom -1 && !strcmp(chrom, chromosomes[prev_chrom_id+1]->name)){
    chrom_id = prev_chrom_id + 1;
  }
  else{
    chrom_id = binSearch(chrom);
  }

  return chrom_id;

}

int binSearch(char *chrom){
  int start;
  int end;
  int med;
  int strtest;

  start = 0;

  end = num_chrom - 1;

  med = (start + end) / 2;


  while (1){

    if (start > end)
      return -1;

    if (start == med || end == med){
      if (!strcmp(chrom, chromosomes[start]->name))
	return start;
      else if (!strcmp(chrom, chromosomes[end]->name))
	return end;
      else{
	return -1;
      }
    }

    strtest = strcmp(chrom, chromosomes[med]->name);

    if(strtest == 0)
      return med;

    else if (strtest < 0){
      end = med;
      med = (start + end) / 2;
    }

    else {
      start = med;
      med = (start + end) / 2;
    }

  }
}



int windowSearch(struct window *searchWin, int winCnt, int pos){
  int start;
  int end;
  int med;

  start = 0;
  end = winCnt - 1;

  med = (start + end) / 2;

  while (1){

    if (start > end)
      return -1;

    if (start == med || end == med){

      if (pos >= searchWin[start].start && pos < searchWin[start].end)
	return start;

      else if (pos >= searchWin[end].start && pos < searchWin[end].end)
	return end;

      else return -1;

    }

    if (pos >= searchWin[med].start && pos < searchWin[med].end)
      return med;

    else if (pos < searchWin[med].start){
      end = med;
      med = (start + end) / 2;
    }

    else {
      start = med;
      med = (start + end) / 2;
    }
  }
}



void saveDepth(char *depthFile, int is_norm){
  int i;
  int j;
  int retVal = 1;

  char *fname;
  FILE *txtDepth;
  FILE *binDepth;

	
  if (is_norm){
    fname = (char *) getMem (sizeof(char) * (strlen(depthFile) + 20));
    binDepth = my_fopen(depthFile, "w", 0);
    /* start with the magicNum, I will use this as a format check when loading */
    retVal = fwrite(&magicNumDepth, sizeof(magicNumDepth), 1, binDepth);
    for (i = 0; i < num_chrom; i++){
      for (j = 0; j < chromosomes[i]->lw_cnt; j++){
	retVal = fwrite(&(chromosomes[i]->lw[j].depth), sizeof(chromosomes[i]->lw[j].depth), 1, binDepth);
      }
    }
    
    for (i = 0; i < num_chrom; i++){
      for (j = 0; j < chromosomes[i]->sw_cnt; j++){
	retVal = fwrite(&(chromosomes[i]->sw[j].depth), sizeof(chromosomes[i]->cw[j].depth), 1, binDepth);
      }
    }
    
    for (i = 0; i < num_chrom; i++){
      for (j = 0; j < chromosomes[i]->cw_cnt; j++){
	retVal = fwrite(&(chromosomes[i]->cw[j].depth), sizeof(chromosomes[i]->cw[j].depth), 1, binDepth);
      }
    }
    
    fclose(binDepth);
    
    fprintf (stdout, "Writing normalized CW depth to: %s.cw_norm.bed.\n", depthFile);

    sprintf (fname, "%s.cw_norm.bed", depthFile);
    dump_text_windows(fname, CW);
    
    fprintf (stdout, "Writing normalized LW depth to: %s.lw_norm.bed.\n", depthFile);
    
    sprintf (fname, "%s.lw_norm.bed", depthFile);
    dump_text_windows(fname, LW);
    
    fprintf (stdout, "Writing normalized SW depth to: %s.sw_norm.bed.\n", depthFile);
    
    sprintf (fname, "%s.sw_norm.bed", depthFile);
    dump_text_windows(fname, SW);
    freeMem(fname, sizeof(char) * (strlen(depthFile) + 20));
    return;
  }

  fname = (char *) getMem (sizeof(char) * (strlen(depthFile) + 10));
  sprintf(fname, "%s.lw.txt", depthFile);
  txtDepth = my_fopen(fname, "w", 0);

  fprintf(stdout, "Saving depth file %s\n", fname);

  for (i = 0; i < num_chrom; i++){
    for (j = 0; j < chromosomes[i]->lw_cnt; j++){
      fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\n", chromosomes[i]->name, chromosomes[i]->lw[j].start, chromosomes[i]->lw[j].end, chromosomes[i]->lw[j].gc, (chromosomes[i]->lw[j].depth));
    }
  }

  fclose(txtDepth);

  sprintf(fname, "%s.sw.txt", depthFile);
  txtDepth = my_fopen(fname, "w", 0);
  fprintf(stdout, "Saving depth file %s\n", fname);

  for (i = 0; i < num_chrom; i++){
    for (j = 0; j < chromosomes[i]->sw_cnt; j++){
      fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\n", chromosomes[i]->name, chromosomes[i]->sw[j].start, chromosomes[i]->sw[j].end,  chromosomes[i]->sw[j].gc, (chromosomes[i]->sw[j].depth));
    }
  }

  fclose(txtDepth);

  sprintf(fname, "%s.cw.txt", depthFile);
  txtDepth = my_fopen(fname, "w", 0);
  fprintf(stdout, "Saving depth file %s\n", fname);

  for (i = 0; i < num_chrom; i++){
    for (j = 0; j < chromosomes[i]->cw_cnt; j++){
      fprintf(txtDepth, "%s\t%d\t%d\t%f\t%f\n", chromosomes[i]->name, chromosomes[i]->cw[j].start, chromosomes[i]->cw[j].end,  chromosomes[i]->cw[j].gc, (chromosomes[i]->cw[j].depth));
    }
  }

  fclose(txtDepth);

  freeMem(fname, sizeof(char) * (strlen(depthFile) + 10));

  if (VERBOSE && retVal == 0)
    fprintf(stderr, "There was potentially an error in fwrite during saveDepth.\n");
  
}
