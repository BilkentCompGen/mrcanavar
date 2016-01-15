#include "mrcanavar.h"



int main(int argc, char **argv){

  int i;
  int ctr;
  int fileCount;
  int fileEnd;
  int fileStart;

  char gzSAM;

  char *indirSAM;
  char *depthFile;
  char *out_prefix;
  char *gene;
  char **c_depths;

  init_globals();

  gzSAM      = 0;
  indirSAM   = NULL;
  depthFile  = NULL;
  out_prefix = NULL;
  gene=NULL;
  n_sam_e = 0;
  n_sam_b =0;
  fileCount = 0;
  fileEnd = 0; 
  fileStart = 1;
  
  c_depths = NULL;

  for (i=1; i<argc; i++){


    /* modes */
    if (!strcmp(argv[i], "--prep"))
      set_runmode(PREP);
    else if (!strcmp(argv[i], "--read"))
      set_runmode(READSAM);
    else if (!strcmp(argv[i], "--call"))
      set_runmode(CALL);
    else if (!strcmp(argv[i], "--conc"))
      set_runmode(CONC);

    /* compulsory parameters for all modes */
    else if (!strcmp(argv[i], "-conf"))
      set_str(&GENOME_CONF, argv[i+1]);

    /* compulsory parameters for the PREP mode */
    else if (!strcmp(argv[i], "-fasta"))
      set_str(&GENOME_FASTA, argv[i+1]);
    else if (!strcmp(argv[i], "-gaps"))
      set_str(&GENOME_GAPS, argv[i+1]);

    /* optional parameters for the PREP mode */
    else if (!strcmp(argv[i], "-pseudoa"))
      set_str(&GENOME_PSEUDO , argv[i+1]);
    else if (!strcmp(argv[i], "-lw_size"))
      LW_SIZE = atoi(argv[i+1]);
    else if (!strcmp(argv[i], "-sw_size"))
      SW_SIZE = atoi(argv[i+1]);
    else if (!strcmp(argv[i], "-cw_size"))
      CW_SIZE = atoi(argv[i+1]);
    else if (!strcmp(argv[i], "-lw_slide"))
      LW_SLIDE = atoi(argv[i+1]);
    else if (!strcmp(argv[i], "-sw_slide"))
      SW_SLIDE = atoi(argv[i+1]);


    /* compulsory parameters for  READ, CALL and CONC modes */
    else if (!strcmp(argv[i], "-depth"))
      set_str(&depthFile, argv[i+1]);
    
    /* compulsory paramaters for CONC mode*/
    else if (!strcmp(argv[i], "-concdepth")){
      
      fileStart=i+1; fileEnd=i+1; i++;

      while (i<argc && argv[i][0]!='-')
	fileEnd=i++;
      i--;

      fileCount = fileEnd-fileStart+1;

      if (fileEnd == fileStart ){
	fprintf(stderr, "No files are given.\n");
	return 0;
      }

      c_depths = getMem( fileCount * sizeof * c_depths);

      for ( ctr=0; ctr<fileCount; ctr++){
	set_str(&c_depths[ctr], argv[ctr+fileStart]);
      }

    }	 
      
 

    /* compulsory parameters for the READ mode */

    else if (!strcmp(argv[i], "-samdir"))
      set_str(&indirSAM, argv[i+1]);
  
  
    /* optional parameters for the SAM mode */
    else if (!strcmp(argv[i], "--gz"))
      gzSAM = 1;
    /* optional parameters for the SAM mode */
    else if (!strcmp(argv[i], "-nsam")){
      n_sam_b = atoi(argv[i+1]);
      n_sam_e = atoi(argv[i+2]);
    }
    /* compulsory parameters for the CALL mode */
    else if (!strcmp(argv[i], "-o"))
      set_str(&out_prefix, argv[i+1]);

    /* optional parameters for the CALL mode */
    else if (!strcmp(argv[i], "-cont_win"))
      CONT_WINDOW = atoi(argv[i+1]);
    else if (!strcmp(argv[i], "-cut_win"))
      CUT_WINDOW = atoi(argv[i+1]);
    else if (!strcmp(argv[i], "--xx"))
      set_gender(FEMALE);
    else if (!strcmp(argv[i], "--xy"))
      set_gender(MALE);
    else if (!strcmp(argv[i], "--multgc"))
      MULTGC = 1;
    else if (!strcmp(argv[i], "-mindup"))
      MIN_DUP = atoi(argv[i+1]);
    else if (!strcmp(argv[i], "-gene"))
      set_str(&gene, argv[i+1]);


    /* generic stuff */
    else if (!strcmp(argv[i], "-v")){
      fprintf (stdout, "\nmrCaNaVaR version %s\nLast update: %s\n", VERSION, LAST_UPDATE);
      return 0;
    }
    else if (!strcmp(argv[i], "-h")){
      printHelp(argv[0]);
      return 0;
    }
    else if (!strcmp(argv[i], "--verbose"))
      VERBOSE = 1;
    else if (!strcmp(argv[i], "--disable-sam-check"))
      CHECKSAM = 0;
  }


  switch(RUNMODE){
  case NONE:
    print_error("No run mode selected.\nSelect mode: --prep, ---read,--call or --conc\n");
    break;
  case PREP:
    fprintf(stdout, "Mode: Prepare genome...\n");
    prep_genome();
    break;
  case READSAM:
    fprintf(stdout, "Mode: Read SAM ...\n");
    read_mapfiles(indirSAM, depthFile, gzSAM);
    break;
  case CALL:
    fprintf(stdout, "Mode: Call copy numbers ...\n");
    call_cnv(depthFile, out_prefix,gene);
    break;
  case CONC:
    fprintf(stdout, "Mode: Merge depth files...\n");
    conc_depth(c_depths, fileCount, depthFile);
    break;

  default:
    break;
  }

  return 0;

}


void printHelp(char *binfile){

  fprintf (stdout, "\nmrCaNaVaR version %s\nLast update: %s\n", VERSION, LAST_UPDATE);
  fprintf (stdout, "%s --<mode>  [options] \n\n", binfile);

  fprintf (stdout, "========            RUN MODES            ========\n\n");

  fprintf (stdout, "\t--prep                  : Prepare reference genome configuration file.\n");
  fprintf (stdout, "\t--read                  : Read mapping information from SAM files.\n");
  fprintf (stdout, "\t--conc                  : Merge .depth files that are given as parameters (optional mode).\n");
  fprintf (stdout, "\t--call                  : Call CNVs and predict copy numbers.\n\n");

  fprintf (stdout, "======== PREP MODE COMPULSORY PARAMETERS ========\n\n");
  fprintf (stdout, "\t-fasta <fasta_file>     : FASTA file for the reference genome.\n");
  fprintf (stdout, "\t-gaps  <gaps_file>      : Gap coordinates of the reference genome in BED format.\n");
  fprintf (stdout, "\t-conf  <config_file>    : Reference configuration file (output).\n\n");

  fprintf (stdout, "========  PREP MODE OPTIONAL PARAMETERS  ========\n\n");
  fprintf (stdout, "\t-lw_size   <lw_size>    : Long window span size. Default is 5000.\n");
  fprintf (stdout, "\t-lw_slide  <lw_slide>   : Long window slide size. Default is 1000.\n");
  fprintf (stdout, "\t-sw_size   <sw_size>    : Short window span size. Default is 1000.\n");
  fprintf (stdout, "\t-sw_slide  <sw_slide>   : Short window slide size. Default is 1000.\n");
  fprintf (stdout, "\t-cw_size   <cw_size>    : Copy number window size. Default is 1000.\n");
  fprintf (stdout, "\t-pseudoa   <bed_file>   : Coordinates for pseudoautosomal regions in the reference genome in BED format.\n\n");

  fprintf (stdout, "======== READ MODE COMPULSORY PARAMETERS ========\n\n");
  fprintf (stdout, "\t-conf   <config_file>   : Reference configuration file (input).\n");
  fprintf (stdout, "\t-samdir <sam_dir>       : Directory that contains SAM files for mapping information.\n");
  fprintf (stdout, "\t-depth  <depth_file>    : Read depth file (output).\n\n");
  fprintf (stdout, "========  READ MODE OPTIONAL PARAMETERS  ========\n\n");
  fprintf (stdout, "\t--gz                    : Indicates that the SAM files are compressed in gzip format.\n");
  fprintf (stdout, "\t-nsam <first> <second>  : Choose a sequence of SAM files between two indexes.\n\n");

  fprintf (stdout, "======== CONC MODE COMPULSORY PARAMETERS ========\n\n");
  fprintf (stdout, "\t-conf     <config_file> : Reference configuration file (input).\n");
  fprintf (stdout, "\t-concdepth <depth_file> : Read depth file(s).\n");
  fprintf (stdout, "\t-depth     <depth_file> : Read depth file (output).\n\n");

  fprintf (stdout, "======== CALL MODE COMPULSORY PARAMETERS ========\n\n");
  fprintf (stdout, "\t-conf   <config_file>   : Reference configuration file (input).\n");
  fprintf (stdout, "\t-depth  <depth_file>    : Read depth file (input).\n");
  fprintf (stdout, "\t-o      <out_prefix>    : Prefix for the output file names.\n\n");

  fprintf (stdout, "========  CALL MODE OPTIONAL PARAMETERS  ========\n\n");
  fprintf (stdout, "\t--xx                    : Set gender of the sequenced sample as female. Mammalian genomes only. Default is autodetect.\n");
  fprintf (stdout, "\t--xy                    : Set gender of the sequenced sample as male. Mammalian genomes only. Default is autodetect.\n");
  fprintf (stdout, "\t--multgc                : Perform multiplicative GC correction. Default is additive.\n");
  fprintf (stdout, "\t-mindup <min_dup_len>   : Minimum duplication length. Default is 10000.\n");
  fprintf (stdout, "\t-gene <genelist>        : Coordinates for genes for calculating gene-based copy numbers.\n");
  fprintf (stdout, "\t--verbose               : Verbose output.\n\n");


  /*

    fprintf (stdout, "========       NOT IMPLEMENTED YET       ========\n\n");
    fprintf (stdout, "\t-cont_win <cwin_number> : Contiguous window number to look for high/low depth windows. Default is 7.\n");
    fprintf (stdout, "\t-cut_win  <cwin_cutoff> : Window number cutoff. Default is 6.\n\n");


  */

}
