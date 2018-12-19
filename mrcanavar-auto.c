/* this is an autopilot for mrsfast mapping & mrcanavar calling. Intended for automated use in cloud. It works locally, but assumes a lot of things. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "utils.h"

char *ref_genome;
char *input_files;
char *cnvr_conf;
char *gene_list;
int num_threads;
int skip_map;
char is_cloud[64];
int dry_run;

void print_help(void){
  fprintf(stderr, "Usage:\n\n");
  fprintf(stderr, "mrcanavar-auto\n\t--input [list]: comma-separated FASTQ(.gz) files.\n\t--ref [ref.fa]: repeat masked reference genome. mrsFAST index (ref.fa.index) should be in the same directory.\n");
  fprintf(stderr, "\t--conf [ref.cnvr]: mrCaNaVaR config file.\n\t--gene [genes.bed]: List of genes.\n\t--threads [int]: number of threads for mrsFAST.\n");
  fprintf(stderr, "\t--skip-mapping: Skip mrsFAST mapping. Use this only if you have the mapping output and rerunning mrCaNaVaR for some reason.\n");
  fprintf(stderr, "\t--cloud: Cloud mode, the directory info from the input file names will be stripped for output file generation.\n");
  fprintf(stderr, "\t--dry-run: Do not execute commands, just print them.\n\n");
}

int parse_command_line( int argc, char** argv)
{
  
  int o;
  int index;
  static int do_cloud = 0;

  
  ref_genome = NULL;  input_files = NULL;  cnvr_conf = NULL;  gene_list = NULL;
  
  static struct option long_options[] =
    {
      {"input" , required_argument, 0, 'i'},
      {"ref" , required_argument, 0, 'f'},
      {"conf" , required_argument, 0, 'c'},
      {"gene" , required_argument, 0, 'g'},    
      {"threads" , required_argument, 0, 't'},    
      {"skip-mapping" , no_argument, &skip_map, 1},    
      {"cloud" , no_argument, &do_cloud, 1},    
      {"dry-run" , no_argument, &dry_run, 1},    
      {0 , 0, 0, 0 }
    };
  
  if( argc == 1)
    {
      print_help();
      return -1;
    }
  
  while( ( o = getopt_long( argc, argv, "i:f:c:g:t", long_options, &index)) != -1)
    {
      switch( o)
	{
	case 'i':
	  set_str(&input_files, optarg);
	  break;
	case 'f':
	  set_str(&ref_genome, optarg);
	  break;
	case 'g':
	  set_str(&gene_list, optarg);
	  break;
	case 'c':
	  set_str(&cnvr_conf, optarg);
	  break;
	case 't':
          num_threads = atoi(optarg);
	  break;
	}
      
    }
  
  if (input_files == NULL){
    fprintf (stderr, "Input files (FASTQ) are missing.\n");
    return -1;
  }
  if (ref_genome == NULL){
    fprintf (stderr, "Reference genome (repeatmasked FASTA) is missing.\n");
    return -1;
  }
  if (gene_list == NULL){
    fprintf (stderr, "Gene list (BED) is missing.\n");
    return -1;
  }
  if (cnvr_conf == NULL){
    fprintf (stderr, "Config file (CNVR) is missing.\n");
    return -1;
  }

  if (do_cloud)
    strcpy(is_cloud, "--cloud");
  
  return 0;
}

char * pacman_directories(char *fname){
  char *last_slash;
  //fprintf(stderr, "Eating %s  ->  ", fname);
  last_slash = strrchr(fname, '/');
  if (last_slash != NULL){
    //fprintf(stderr, "%s\n", last_slash+1);
    return last_slash+1;
  }
  else
    ;
  //fprintf(stderr, "%s\n", fname);

  return fname;
}

int main (int argc, char **argv){
  int ret;
  int i = 0;
  
  int number_of_fastq = 0;
  char **fastq_files;
  char *safe_fastq_list = NULL;
  char *token;
  
  char cmd_line[1048576];
  char tmp_line[1048576];
  is_cloud[0] = 0;

  dry_run = 0;
  skip_map = 0;
  num_threads = 4;
  
  ret = parse_command_line(argc, argv);

  if (ret == -1)
    return ret;
  
  set_str(&safe_fastq_list, input_files);
  token = strtok(safe_fastq_list, ",");
  while (token != NULL){
    number_of_fastq++;
    token = strtok(NULL, ",");
  }
  
  fastq_files = (char **) getMem(sizeof(char *) * number_of_fastq);
  memset(fastq_files, 0, sizeof(char *) * number_of_fastq);
  
  set_str(&safe_fastq_list, input_files);
  token = strtok(safe_fastq_list, ",");
  while(token != NULL){
    set_str(&(fastq_files[i++]), token);
    token = strtok (NULL, ",");
  }

  /* map */
  for (i=0; i<number_of_fastq; i++){
    
    if (endswith (fastq_files[i], ".gz")) 
      sprintf(cmd_line, "mrsfast %s --search %s --seq %s --seqcomp --outcomp -t %d --mem 6 --crop 36 -e 2", is_cloud, ref_genome, fastq_files[i], num_threads);
    else
      sprintf(cmd_line, "mrsfast %s --search %s --seq %s --outcomp -t %d --mem 6 --crop 36 -e 2", is_cloud,  ref_genome, fastq_files[i], num_threads);

    fprintf(stderr, "[MAPPING] %s\n", cmd_line);

    if (!skip_map && !dry_run)
      ret = system(cmd_line);
  }

  fprintf(stderr, "[PROGRESS] Read mapping done. Now reading SAM files.\n");


  if (is_cloud[0] == '-'){
    for (i=0; i<number_of_fastq; i++)
      fastq_files[i] = pacman_directories(fastq_files[i]);
  }
  
  /* read */

  if (number_of_fastq == 1){
    sprintf(cmd_line, "mrcanavar --read --gz -conf %s -samlist %s.sam.gz -depth %s.depth", cnvr_conf, fastq_files[0], fastq_files[0]);
  }
  else{
    tmp_line[0] = 0;
    sprintf(tmp_line, "%s", fastq_files[0]);
    for (i=1; i<number_of_fastq; i++){
      strcat(tmp_line, ".sam.gz,");
      strcat(tmp_line, fastq_files[i]);
    }
    strcat(tmp_line, ".sam.gz");

  sprintf(cmd_line, "mrcanavar --read --gz -conf %s -samlist %s -depth %s.depth", cnvr_conf, tmp_line, fastq_files[0]);

  }

  fprintf(stderr, "[PROGRESS] SAM files loaded. Now calling copy numbers and CNVs.\n");
  
  fprintf(stderr, "[READ SAM FILES] %s\n", cmd_line);
  if (!dry_run)
    ret = system(cmd_line);
  /* call */

  sprintf(cmd_line, "mrcanavar --call -conf %s -depth %s.depth -gene %s -o %s-out", cnvr_conf, fastq_files[0], gene_list, fastq_files[0]);
  fprintf(stderr, "[CALL COPY NUMBERS] %s\n", cmd_line);
  if (!dry_run)
    ret = system(cmd_line);

  fprintf(stderr, "[PROGRESS] Done.\n");
}
