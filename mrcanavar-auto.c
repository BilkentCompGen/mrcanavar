
/* this is an autopilot for mrsfast mapping & mrcanavar calling. Intended for automated use in cloud. It works locally, but assumes a lot of things. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "utils.h"

char *ref_genome;
char *unmasked_ref_genome;
char *input_files;
char *cnvr_conf;
char *gene_list;
int num_threads;
int skip_map;
char is_cloud[64];
int dry_run;
char *bam_input;
int mem;
int uncompressed;
int nosam;
int kmerlen;

void print_help(void){
  fprintf(stderr, "Usage:\n\n");
  fprintf(stderr, "mrcanavar-auto\n\t--input [list]: comma-separated FASTQ(.gz) files.\n\t--aln-input [list]: comma-separated BAM/CRAM files (alternative to --input).\n\t--ref [ref.fa]: repeat masked reference genome. mrsFAST index (ref.fa.index) should be in the same directory.\n");
  fprintf(stderr, "\t--conf [ref.cnvr]: mrCaNaVaR config file.\n\t--gene [genes.bed]: List of genes.\n\t--threads [int]: number of threads for mrsFAST.\n");
  fprintf(stderr, "\t--kmer [int]: Cropping length for mrsFAST mapping. Default is 36.\n");
  fprintf(stderr, "\t--skip-mapping: Skip mrsFAST mapping. Use this only if you have the mapping output and rerunning mrCaNaVaR for some reason.\n");
  fprintf(stderr, "\t--cloud: Cloud mode, the directory info from the input file names will be stripped for output file generation.\n");
  fprintf(stderr, "\t--no-gz: Do not compress mrsFAST output. This option will generate larger files, but it will save some run time.\n");
  fprintf(stderr, "\t--no-sam: Do not generate SAM(.gz) output files. Pipe the mrsFAST output directly into mrCaNaVaR. This option does not work with multiple FASTQ files as input.\n");
  fprintf(stderr, "\t--dry-run: Do not execute commands, just print them.\n\n");
}

int parse_command_line( int argc, char** argv)
{
  
  int o;
  int index;
  static int do_cloud = 0;

  
  ref_genome = NULL;  input_files = NULL;  cnvr_conf = NULL;  gene_list = NULL; bam_input = NULL; unmasked_ref_genome = NULL;
  
  static struct option long_options[] =
    {
      {"input" , required_argument, 0, 'i'},
      {"aln-input" , required_argument, 0, 'a'},
      {"ref" , required_argument, 0, 'f'},
      {"unmasked-ref" , required_argument, 0, 'u'},
      {"conf" , required_argument, 0, 'c'},
      {"gene" , required_argument, 0, 'g'},    
      {"threads" , required_argument, 0, 't'},    
      {"mem" , required_argument, 0, 'm'},    
      {"kmer" , required_argument, 0, 'k'},    
      {"skip-mapping" , no_argument, &skip_map, 1},    
      {"no-gz" , no_argument, &uncompressed, 1},    
      {"no-sam" , no_argument, &nosam, 1},    
      {"dry-run" , no_argument, &dry_run, 1},    
      {"cloud" , no_argument, &do_cloud, 1},    
      {0 , 0, 0, 0 }
    };
  
  if( argc == 1)
    {
      print_help();
      return -1;
    }
  
  while( ( o = getopt_long( argc, argv, "i:a:f:u:c:g:m:k:t", long_options, &index)) != -1)
    {
      switch( o)
	{
	case 'i':
	  set_str(&input_files, optarg);
	  break;
	case 'a':
	  set_str(&bam_input, optarg);
	  break;
	case 'f':
	  set_str(&ref_genome, optarg);
	  break;
	case 'u':
	  set_str(&unmasked_ref_genome, optarg);
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
	case 'm':
          mem = atoi(optarg);
	  break;
	case 'k':
          kmerlen = atoi(optarg);
	  break;
	}
      
    }
  
  if (input_files == NULL && bam_input == NULL){
    fprintf (stderr, "Input files (FASTQ or BAM/CRAM) are missing.\n");
    return -1;
  }
  if (ref_genome == NULL){
    fprintf (stderr, "Reference genome (repeatmasked FASTA) is missing.\n");
    return -1;
  }
  if (bam_input != NULL && unmasked_ref_genome == NULL){
    fprintf (stderr, "BAM/CRAM input requires unmasked reference genome, but it is missing.\n");
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
  char out_file[1024];
  
  int bam_mode;
  char outcomp[32];

  kmerlen = 36;
  nosam = 0;
  mem = 6;
  
  is_cloud[0] = 0;

  dry_run = 0;
  skip_map = 0;
  uncompressed = 0;
  num_threads = 4;
  
  ret = parse_command_line(argc, argv);

  
  if (ret == -1)
    return ret;

  if (bam_input != NULL){
    bam_mode = 1;
    fprintf (stderr, "BAM/CRAM mode.\n");
  }
  else{
    fprintf (stderr, "FASTQ mode.\n");
    bam_mode = 0;
  }

  if (uncompressed){
    outcomp[0] = 0;
    fprintf(stderr, "Uncompressed SAM files.\n");
  }
  else{
    strcpy(outcomp, "--outcomp");
    fprintf(stderr, "Compressing SAM files.\n");
  }



  if ( input_files != NULL)
    set_str(&safe_fastq_list, input_files);
  else if ( bam_input != NULL)
    set_str(&safe_fastq_list, bam_input);
  else{
    fprintf (stderr, "Input files (FASTQ or BAM/CRAM) are missing.\n");
    return -1;
  }


  token = strtok(safe_fastq_list, ",");
  while (token != NULL){
    number_of_fastq++;
    token = strtok(NULL, ",");
  }
  
  fastq_files = (char **) getMem(sizeof(char *) * number_of_fastq);
  memset(fastq_files, 0, sizeof(char *) * number_of_fastq);
  
  if ( input_files != NULL)
    set_str(&safe_fastq_list, input_files);
  else if ( bam_input != NULL)
    set_str(&safe_fastq_list, bam_input);
  else{
    fprintf (stderr, "Input files (FASTQ or BAM/CRAM) are missing.\n");
    return -1;
  }
  
  token = strtok(safe_fastq_list, ",");
  while(token != NULL){
    set_str(&(fastq_files[i++]), token);
    token = strtok (NULL, ",");
  }


  if (nosam && number_of_fastq > 1){
    fprintf(stderr, "[WARNING]: --no-sam option works only one FASTQ or BAM/CRAM file is supplied. Removing the parameter from consideration.\n");
    nosam = 0;
  }

  /* map */
  if (!nosam){
    for (i=0; i<number_of_fastq; i++){
      
      if ( bam_mode == 0 ){
	if (endswith (fastq_files[i], ".gz")){
	  sprintf(cmd_line, "mrsfast %s --search %s --seq %s --seqcomp %s -t %d --mem %d --crop %d -e 2 --disable-nohits", is_cloud, ref_genome, fastq_files[i], outcomp, num_threads, mem, kmerlen);
	}
	else
	  sprintf(cmd_line, "mrsfast %s --search %s --seq %s %s -t %d --mem %d --crop %d -e 2 --disable-nohits", is_cloud,  ref_genome, fastq_files[i], outcomp, num_threads, mem, kmerlen);
      }
      else if (number_of_fastq == 1){
	/* samtools and pipe */
	sprintf(cmd_line, "samtools view -@ %d -T %s %s | awk '{print \">\"$1\"\\n\"$10}' | mrsfast %s --search %s --seq /dev/stdin %s -t %d --mem %d --crop %d -e 2 -o %s --disable-nohits", num_threads, unmasked_ref_genome, fastq_files[i], is_cloud,  ref_genome,  outcomp, num_threads, mem, kmerlen, fastq_files[i]);
      }
      
      fprintf(stderr, "[MAPPING] %s\n", cmd_line);
      
      if (!skip_map && !dry_run)
	ret = system(cmd_line);
    }
    
    fprintf(stderr, "[PROGRESS] Read mapping done. Now reading SAM files.\n");
    
    
    if (is_cloud[0] == '-'){
      for (i=0; i<number_of_fastq; i++)
	fastq_files[i] = pacman_directories(fastq_files[i]);
    }

    strcpy(out_file, fastq_files[0]);
    fprintf(stderr, "[OUT_FILE] %s\n", out_file);
    
    /* read */
    
    if (number_of_fastq == 1){
      if (uncompressed)
	sprintf(cmd_line, "mrcanavar --read -conf %s -samlist %s.sam -depth %s.depth", cnvr_conf, fastq_files[0], out_file);
      else
	sprintf(cmd_line, "mrcanavar --read --gz -conf %s -samlist %s.sam.gz -depth %s.depth", cnvr_conf, fastq_files[0], out_file);
    }
    else{
      tmp_line[0] = 0;
      sprintf(tmp_line, "%s", fastq_files[0]);
      for (i=1; i<number_of_fastq; i++){
	if (uncompressed)
	  strcat(tmp_line, ".sam,");
	else
	  strcat(tmp_line, ".sam.gz,");
	strcat(tmp_line, fastq_files[i]);
      }
      if (uncompressed)
	strcat(tmp_line, ".sam");
      else
	strcat(tmp_line, ".sam.gz");
      
      if (uncompressed)
	sprintf(cmd_line, "mrcanavar --read -conf %s -samlist %s -depth %s.depth", cnvr_conf, tmp_line, out_file);
      else
	sprintf(cmd_line, "mrcanavar --read --gz -conf %s -samlist %s -depth %s.depth", cnvr_conf, tmp_line, out_file);
      
    }

    fprintf(stderr, "[PROGRESS] SAM files generated. Now calling copy numbers and CNVs.\n");
    
    fprintf(stderr, "[READ SAM FILES] %s\n", cmd_line);
    if (!dry_run)
      ret = system(cmd_line);

  }

  else{  // pipe directly into mrcanavar
    if ( bam_mode == 0 ){
      if (endswith (fastq_files[i], ".gz")){
	sprintf(tmp_line, "mrsfast %s --search %s --seq %s --seqcomp -t %d --mem %d --crop %d -e 2 --disable-nohits -o /dev/stdout ", is_cloud, ref_genome, fastq_files[i], num_threads, mem, kmerlen);
      }
      else
	sprintf(tmp_line, "mrsfast %s --search %s --seq %s -t %d --mem %d --crop %d -e 2 --disable-nohits -o /dev/stdout", is_cloud,  ref_genome, fastq_files[i], num_threads, mem, kmerlen);
    }
    else { // bam/cram to mrcanavar
      sprintf(tmp_line, "samtools view -@ %d -T %s %s | awk '{print \">\"$1\"\\n\"$10}' | mrsfast %s --search %s --seq /dev/stdin -t %d --mem %d --crop %d -e 2 -o /dev/stdout --disable-nohits", num_threads, unmasked_ref_genome, fastq_files[0], is_cloud,  ref_genome, num_threads, mem, kmerlen);      
    }
    
    sprintf(cmd_line, "%s | mrcanavar --read -conf %s -samstdin -depth %s.depth", tmp_line, cnvr_conf, out_file);

    fprintf(stderr, "[MAP AND READ] %s\n", cmd_line);
    if (!dry_run)
      ret = system(cmd_line);
    
  }
  
  /* call */

  sprintf(cmd_line, "mrcanavar --call -conf %s -depth %s.depth -gene %s -o %s-out", cnvr_conf, out_file, gene_list, out_file);
  fprintf(stderr, "[CALL COPY NUMBERS] %s\n", cmd_line);
  if (!dry_run)
    ret = system(cmd_line);

  fprintf(stderr, "[PROGRESS] Done.\n");
}
