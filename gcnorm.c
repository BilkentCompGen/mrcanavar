#include "gcnorm.h"
#define STDMULT 3

int norm_until_converges (enum WINDOWTYPE wt, float *gclookup, float *gclookup_x){
  float max;
  float max_x;
  float min;
  float min_x;

  float p_max;
  float p_max_x;
  float p_min;
  float p_min_x;

  int i, j;
  float mean    = 0.0;
  float mean_x  = 0.0;
  float stdev   = 0.0;
  float stdev_x = 0.0;
  int iter;

  float maxcut;
  float maxcut_x;
  float mincut;
  float mincut_x;

  iter = 1;
  calc_stat(wt, gclookup, gclookup_x, 0, &max, &max_x, &min, &min_x);

  switch(wt){
  case LW:
    mean    = LW_MEAN;
    mean_x  = LW_MEAN_X;
    break;
  case SW:
    mean    = SW_MEAN;
    mean_x  = SW_MEAN_X;
    break;
  case CW:
    mean    = CW_MEAN;
    mean_x  = CW_MEAN_X;
    break;
  }


  if (GENDER == AUTODETECT){
    if (VERBOSE)
      fprintf(stdout, "MEAN: %f\tMEAN_X: %f\n", mean, mean_x);


    if (mean_x / mean < 0.75){ // magical ratio
      GENDER = MALE;

      if (VERBOSE)
	fprintf(stdout, "Autodetect Gender: Male\n");

      i = findChrom("chrX", "", -1);
      
      if (i==-1)
        i = findChrom("X", "", -1);

      if (i != -1){
	switch(wt){
	  
	case CW:
	  for (j=0; j<chromosomes[i]->cw_cnt; j++)
	    if (chromosomes[i]->cw[j].depth > (float) CW_MEAN_X * 2 || chromosomes[i]->cw[j].depth < (float) CW_MEAN_X / 10.0)
	      chromosomes[i]->cw[j].isControl = 0;
	  break;
	  
	case LW:
	  for (j=0; j<chromosomes[i]->lw_cnt; j++)
	    if (chromosomes[i]->lw[j].depth > (float) LW_MEAN_X * 2 || chromosomes[i]->lw[j].depth < (float) LW_MEAN_X / 10.0)
	      chromosomes[i]->lw[j].isControl = 0;
	  break;
	  
	case SW:
	  for (j=0; j<chromosomes[i]->sw_cnt; j++)
	    if (chromosomes[i]->sw[j].depth > (float) SW_MEAN_X * 2 || chromosomes[i]->sw[j].depth < (float) SW_MEAN_X / 10.0)
	      chromosomes[i]->sw[j].isControl = 0;
	  break;
	  
	}
      }
    }

    else{
      GENDER = FEMALE;
      if (VERBOSE)
	fprintf(stdout, "Autodetect Gender: Female\n");
    }
  }

  if (!ISNORMALIZED)
    normalize(wt, gclookup, gclookup_x, &max, &max_x, &min, &min_x);


  do{

    if (VERBOSE){
      fprintf(stdout, "Control regions %s iteration %d ", (wt == LW ? "LW" : (wt == SW ? "SW" : "CW")), iter);
      fflush(stdout);
    }

    p_min   = min;
    p_max   = max;
    p_min_x = min_x;
    p_max_x = max_x;

    calc_stat(wt, gclookup, gclookup_x, 1, &max, &max_x, &min, &min_x);

    switch(wt){
    case LW:
      mean    = LW_MEAN;
      stdev   = LW_STD;
      mean_x  = LW_MEAN_X;
      stdev_x = LW_STD_X;
      break;
    case SW:
      mean    = SW_MEAN;
      stdev   = SW_STD;
      mean_x  = SW_MEAN_X;
      stdev_x = SW_STD_X;
      break;
    case CW:
      mean    = CW_MEAN;
      stdev   = CW_STD;
      mean_x  = CW_MEAN_X;
      stdev_x = CW_STD_X;
      break;
    }

    if (GENDER == FEMALE){
      max_x   = max;
      mean_x  = mean;
      stdev_x = stdev;
      min_x   = min;
    }

    iter++;

    maxcut   = mean + STDMULT * stdev;
    mincut   = mean - STDMULT * stdev;
    maxcut_x = mean_x + STDMULT * stdev_x;
    mincut_x = mean_x - STDMULT * stdev_x;

    /*
    maxcut   = mean * 2.5;// + STDMULT * stdev;
    mincut   = mean / 2.5; //- STDMULT * stdev;
    maxcut_x = mean_x * 1.5;//+ STDMULT * stdev_x;
    mincut_x = mean_x / 1.5; //- STDMULT * stdev_x;
    */

    /*
    if (GENDER != FEMALE){
      maxcut_x = mean_x + STDMULT * stdev_x;
      mincut_x = mean_x - STDMULT * stdev_x;
      }*/

    if (mincut < 0.0) mincut = mean / 10.0;
    if (mincut_x < 0.0) mincut_x = mean_x / 10.0;

    /*
    if (maxcut - mean > mean - mincut)
      maxcut = mean + (mean - mincut);

    if (GENDER != FEMALE){
      if (maxcut_x - mean_x > mean_x - mincut_x)
	maxcut_x = mean_x + (mean_x - mincut_x);
	}*/


    if (VERBOSE){
      fprintf(stdout, "mean: %f\tstdev: %f\tmax: %f (cut: %f)\tmin: %f (cut: %f) \n", mean, stdev, max, maxcut, min, mincut);
      if (GENDER != FEMALE)
	fprintf(stdout, "mean_x: %f\tstdev_x: %f\tmax_x: %f (cut: %f)\tmin_x: %f (cut: %f)\n", mean_x, stdev_x, max_x, maxcut_x, min_x, mincut_x);
    }

    if (p_min == min && p_max == max && p_min_x == min_x && p_max_x == max_x)
      break;

    //  } while (max >= maxcut || max_x >= maxcut_x || min <= mincut || min_x <= mincut_x);  
  } while (stdev > 0.25 * mean * (0.5 * PLOIDY));


  fprintf(stdout, "%s Normalization completed.\n",  (wt == LW ? "LW" : (wt == SW ? "SW" : "CW")));

  return 1;
}


void normalize (enum WINDOWTYPE wt, float *gclookup, float *gclookup_x, float *max, float *max_x, float *min, float *min_x){

  int i;

  for (i=0; i<num_chrom; i++)
    norm_wins(i, wt, gclookup, gclookup_x);

  calc_stat(wt, gclookup, gclookup_x, 0, max, max_x, min, min_x);

}


void calc_stat(enum WINDOWTYPE wt, float *gclookup, float *gclookup_x, char doClean, float *_max, float *_max_x, float *_min, float *_min_x){
  int i, j, k;

  float lw_var;
  float lw_var_x;

  int norm_win_cnt;
  float max;

  int norm_win_cnt_x;
  float max_x;

  float min;
  float min_x;

  int gc_index;

  float gc_total[GC_BIN]; // total depth by GC
  int gc_wincount[GC_BIN]; // count of windows with the given GC content

  float gc_total_x[GC_BIN]; // total depth by GC
  int gc_wincount_x[GC_BIN]; // count of windows with the given GC content

  float MEAN   = 0.0;
  float MEAN_X = 0.0;


  float maxcut, mincut;
  float maxcut_x, mincut_x;

  struct window *win = NULL;
  float STD          = 0.0;
  float STD_X        = 0.0;
  float this_total   = 0.0;
  float this_total_x = 0.0;
  float win_cnt;

  lw_var = 0.0;
  norm_win_cnt = 0;

  lw_var_x = 0.0;
  norm_win_cnt_x = 0;
  
  win_cnt = 0;

  for (i=0; i<GC_BIN; i++){
    gc_total[i]    = 0.0;
    gclookup[i]    = 0.0;
    gc_wincount[i] = 0;

    gc_total_x[i]    = 0.0;
    gclookup_x[i]    = 0.0;
    gc_wincount_x[i] = 0;
  }

  max = 0.0;
  max_x = 0.0;

  min = FLT_MAX;
  min_x = FLT_MAX;

  if (doClean){
    
    for (i=0; i<num_chrom; i++){
      
      switch (wt){
      case LW:
	win     = chromosomes[i]->lw;
	win_cnt = chromosomes[i]->lw_cnt;
	MEAN    = LW_MEAN; 
	MEAN_X  = LW_MEAN_X;
	STD     = LW_STD;
	STD_X   = LW_STD_X;
	break;
      case SW:
	win     = chromosomes[i]->sw;
	win_cnt = chromosomes[i]->sw_cnt;
	MEAN    = SW_MEAN; 
	MEAN_X  = SW_MEAN_X;
	STD     = SW_STD;
	STD_X   = SW_STD_X;
	break;
      case CW:
	win     = chromosomes[i]->cw;
	win_cnt = chromosomes[i]->cw_cnt;
	MEAN    = CW_MEAN; 
	MEAN_X  = CW_MEAN_X;
	STD     = CW_STD;
	STD_X   = CW_STD_X;
	break;
      }


      for (j=0; j<win_cnt; j++){
	
	if (win[j].isControl == 1){

	  /* AUTOSOMES -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */
	  if (isAutosome(NULL, -1, i)){
	    maxcut = MEAN + STDMULT * STD;
	    mincut = MEAN - STDMULT * STD;
	    if (mincut < MEAN / 10.0) mincut = MEAN / 10.0;

	    /*    if (maxcut - MEAN > MEAN - mincut)
		  maxcut = MEAN + (MEAN - mincut); */

	    if (win[j].depth > maxcut || win[j].depth < mincut){

	      /* Remove this window and its neighbors from controls */
	      win[j].isControl = 0;
		
	      if (wt != CW){
		k = j;

		while (k > 0 && win[k-1].end >= win[j].start){
		  win[--k].isControl = 0;
		}
		  
		k = j;
		  
		while (k < win_cnt-1 && win[k+1].start <= win[j].end){
		  win[++k].isControl = 0;
		}
	      }
	      else{
		if (j > 0)
		  win[j-1].isControl = 0;
		if (j < win_cnt-1)
		  win[j+1].isControl = 0;
	      }
	    }

	    else{
	      this_total += win[j].depth;
	      norm_win_cnt++;
	    }
	  }  // if AUTOSOME


	  /* chrX -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */
	  else{
	    maxcut_x = MEAN_X + 2 * STD_X;
	    mincut_x = MEAN_X - 2 * STD_X;

	    if (mincut_x < MEAN_X / 10.0) mincut_x = MEAN_X / 10.0;

	    if (win[j].depth > maxcut_x || win[j].depth < mincut_x){

	      /* Remove this window and its neighbors from controls */
	      win[j].isControl = 0;

	      if (wt != CW){
		k = j;

		while (k > 0 && win[k-1].end >= win[j].start){
		  win[--k].isControl = 0;
		}
		  
		k = j;
		  
		while (k < win_cnt-1 && win[k+1].start <= win[j].end){
		  win[++k].isControl = 0;
		}
	      }
	      else{
		if (j > 0)
		  win[j-1].isControl = 0;
		if (j < win_cnt-1)
		  win[j+1].isControl = 0;
	      }
	    }

	    else{
	      this_total_x += win[j].depth;
	      norm_win_cnt_x++;
	    }
	  }  // if chrX

	}

      }
    }

    MEAN_X = this_total_x / norm_win_cnt_x;
    MEAN = this_total / norm_win_cnt;

    switch (wt){
    case LW:
      LW_MEAN    = MEAN; 
      LW_MEAN_X  = MEAN_X;
      break;
    case SW:
      SW_MEAN    = MEAN; 
      SW_MEAN_X  = MEAN_X;
      break;
    case CW:
      CW_MEAN    = MEAN; 
      CW_MEAN_X  = MEAN_X;
      break;
    }
   
  }  // do clean

  this_total      = 0.0;
  norm_win_cnt    = 0;
  this_total_x    = 0.0;
  norm_win_cnt_x  = 0;

  for (i=0; i<num_chrom; i++){

    switch (wt){
    case LW:
      win     = chromosomes[i]->lw;
      win_cnt = chromosomes[i]->lw_cnt;
      break;
    case SW:
      win     = chromosomes[i]->sw;
      win_cnt = chromosomes[i]->sw_cnt;
      break;
    case CW:
      win     = chromosomes[i]->cw;
      win_cnt = chromosomes[i]->cw_cnt;
      break;
    }

    for (j=0; j<win_cnt; j++){

      if (win[j].isControl == 1){

	/* AUTOSOMES -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */
	if (isAutosome(NULL, -1, i)){
	  lw_var += (MEAN - win[j].depth) * (MEAN - win[j].depth);
	  norm_win_cnt++;

	  this_total += win[j].depth;
	  gc_index = win[j].gc * GC_BIN;
	  gc_total[gc_index] += win[j].depth;
	  gc_wincount[gc_index]++;

	  if (win[j].depth > max)
	    max = win[j].depth;
	  if (win[j].depth < min)
	    min = win[j].depth;

	} // AUTOSOMES

	  /* chrX -- NOTE: chrY is NOT a control; it is prefiltered, no need to check for it here */
	else{ // if (strstr(chromosomes[i]->name, "chrX") && GENDER != FEMALE){

	  lw_var_x += (MEAN_X - win[j].depth) * (MEAN_X - win[j].depth);
	  norm_win_cnt_x++;

	  this_total_x += win[j].depth;
	  gc_index = win[j].gc * GC_BIN;
	  gc_total_x[gc_index] += win[j].depth;
	  gc_wincount_x[gc_index]++;

	  if (win[j].depth > max_x)
	    max_x = win[j].depth;
	  if (win[j].depth < min_x)
	    min_x = win[j].depth;

	} // chrX


      } // if control
    } // inner for j
  }     //    outer for i

  MEAN_X = this_total_x / norm_win_cnt_x;
  STD_X = sqrt(lw_var_x / norm_win_cnt_x);

  MEAN = this_total / norm_win_cnt;
  STD = sqrt(lw_var / norm_win_cnt);

  switch (wt){
  case LW:
    LW_MEAN    = MEAN; 
    LW_MEAN_X  = MEAN_X;
    LW_STD     = STD;
    LW_STD_X   = STD_X;
    break;
  case SW:
    SW_MEAN    = MEAN; 
    SW_MEAN_X  = MEAN_X;
    SW_STD     = STD;
    SW_STD_X   = STD_X;
    break;
  case CW:
    CW_MEAN    = MEAN; 
    CW_MEAN_X  = MEAN_X;
    CW_STD     = STD;
    CW_STD_X   = STD_X;
    break;
  }
    
    
  /* calculate the gclookup table */

  for (i=0; i<GC_BIN; i++){

    /* AUTOSOMES */

    if (gc_wincount[i] == 0 || gc_total[i] == 0){
      j=i-1;

      while (j >= 0 && gc_total[j] == 0)
	j--;

      if (j >= 0 && (gc_total[j] == 0 || gc_wincount[j] == 0)){
	j=i+1;
	while (j < GC_BIN && gc_total[j] == 0) j++;
      }

      gc_total[i]    = gc_total[j];
      gc_wincount[i] = gc_wincount[j];
    }

    if (gc_total[i] != 0.0){
      gclookup[i] = gc_total[i] / (float) gc_wincount[i];

      if (MULTGC){

	gclookup[i] = MEAN / gclookup[i];

	if (gclookup[i] > MAX_GC_CORR)
	  gclookup[i] = MAX_GC_CORR;
	else if (gclookup[i] < MIN_GC_CORR)
	  gclookup[i] = MIN_GC_CORR;
      }

    }

    /* chrX */

    if (GENDER != FEMALE){

      if (gc_wincount_x[i] == 0 || gc_total_x[i] == 0){
	j=i-1;

	while (j >= 0 && gc_total_x[j] == 0)
	  j--;

    if (j < 0) j = 0;

	if (gc_total_x[j] == 0 || gc_wincount_x[j] == 0){
	  j=i+1;
	  while (j < GC_BIN && gc_total_x[j] == 0) j++;
	}

	gc_total_x[i]    = gc_total_x[j];
	gc_wincount_x[i] = gc_wincount_x[j];
      }

      if (gc_total_x[i] != 0.0){
	gclookup_x[i] = gc_total_x[i] / (float) gc_wincount_x[i];

	if (MULTGC){
	  gclookup_x[i] = MEAN_X / gclookup_x[i];
	  if (gclookup_x[i] > MAX_GC_CORR)
	    gclookup_x[i] = MAX_GC_CORR;
	  else if (gclookup_x[i] < MIN_GC_CORR)
	    gclookup_x[i] = MIN_GC_CORR;
	}

      }

    }


  }


  *_max   = max;
  *_max_x = max_x;

  *_min   = min;
  *_min_x = min_x;

}


void norm_wins(int chrom_id, enum WINDOWTYPE wt, float *gclookup, float *gclookup_x){
  int j;
  int gc_index;
  float new_depth;
  struct window *win = NULL;
  int win_cnt        = 0;
  float mean         = 0.0;
  float mean_x       = 0.0;

  switch (wt){
  case LW:
    win     = chromosomes[chrom_id]->lw;
    win_cnt = chromosomes[chrom_id]->lw_cnt;
    mean    = LW_MEAN; 
    mean_x  = LW_MEAN_X;
    break;
  case SW:
    win     = chromosomes[chrom_id]->sw;
    win_cnt = chromosomes[chrom_id]->sw_cnt;
    mean    = SW_MEAN; 
    mean_x  = SW_MEAN_X;
    break;
  case CW:
    win     = chromosomes[chrom_id]->cw;
    win_cnt = chromosomes[chrom_id]->cw_cnt;
    mean    = CW_MEAN; 
    mean_x  = CW_MEAN_X;
    break;
  }

  for (j=0; j<win_cnt; j++){
    gc_index  = win[j].gc * (float) GC_BIN;
    if (isAutosome(win, j, chrom_id)){
      if (MULTGC)
	new_depth = gclookup[gc_index] * win[j].depth;
      else
	new_depth = win[j].depth - (gclookup[gc_index] - mean);
    }
    else{
      if (MULTGC)
	new_depth = gclookup_x[gc_index] * win[j].depth;
      else
	new_depth = win[j].depth - (gclookup_x[gc_index] - mean_x);
    }
    
    if (new_depth < 0) new_depth = 0;
    
    win[j].depth = new_depth;
  }
  
}

