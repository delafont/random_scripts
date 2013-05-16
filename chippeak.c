/*
  chip_peak.c

  Signal Peak Tool.
  The program locates signal peaks within a SGA file.

  # Arguments:
  # feature type, strand, integration range (Win1)
  # minimal distance (Win2), counts threshold (Thres)

  Giovanna Ambrosini, EPFL/ISREC, Giovanna.Ambrosini@epfl.ch

#define DEBUG
*/
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <ctype.h>
#ifdef DEBUG
#include <mcheck.h>
#endif

/*#define BUF_SIZE 4096 */
#define BUF_SIZE 8192
#define LINE_SIZE 1024
#define FT_MAX  64
#define SEQ_ID  32
#define POS_MAX 16
#define CNT_MAX 16
#define EXT_MAX 256

typedef struct _options_t {
  int help;
  int debug;
  int refine;
} options_t;

static options_t options;

typedef struct _feature_t {
   char seq_id[SEQ_ID]; // chrom
   char *ft;
   char ft_str;
   char **name;
   int *pos;            // position
   int *cnt;            // read count
   int *npo;
   int *nct;
   int *ptr;
} feature_t, *feature_p_t;

feature_t ref_ft;
int strand_flag = 0;

char *Feature = NULL;
int Win1 = 0; // window
int Win2 = 0; // vicinity
int Thres = 100; // threshold
int Coff = 10; // cutoff

void
locate_peaks(int len)
{
  /* Compute sum in window Win1; store high values in shorter arrays (npo, nct) */
  unsigned long long sum = 0;
  int cnts = 0;
  int i, k;
  int j = 0;
  size_t mLen = BUF_SIZE;
  unsigned int size;
  int *lm;
  char str = '0';
  for (i = 1; i <= len; i++) {
    sum = ref_ft.cnt[i];
    /* Sum up all tag counts within the window range +-Win1/2 */
    for (k = i - 1; ref_ft.pos[k] >= ref_ft.pos[i] - Win1/2; k--) {
      sum += ref_ft.cnt[k];
    }
    for (k = i + 1; ref_ft.pos[k] <= ref_ft.pos[i] + Win1/2; k++) {
      sum += ref_ft.cnt[k];
    }
    if ((unsigned int)j >= mLen - 1) {
      mLen *= 2;
      if (( ref_ft.npo = (int *)realloc(ref_ft.npo, mLen * sizeof(int))) == NULL) {
        perror("locate_peaks: realloc");
        exit(1);
      }
      if (( ref_ft.nct = (int *)realloc(ref_ft.nct, mLen * sizeof(int))) == NULL) {
        perror("locate_peaks: realloc");
        exit(1);
      }
      if (( ref_ft.ptr = (int *)realloc(ref_ft.ptr, mLen * sizeof(int))) == NULL) {
        perror("locate_peaks: realloc");
        exit(1);
      }
    }
    if ( sum >= (unsigned long long)Thres) {
      j++;
      ref_ft.npo[j] = ref_ft.pos[i];
      ref_ft.nct[j] = sum;
      ref_ft.ptr[j] = i;
    }
  }
  /* Initialize Local Maxima Array lm */
  size = (unsigned int)j + 1;
  if (( lm = (int*)calloc(size, sizeof(int))) == NULL) {
    perror("locate_peaks: calloc");
    exit(1);
  }
  /* Keep only one maximum value (peak) within a vicinity range +-Win2/2 */
  /* Record local maxima in lm flag Array */
  /* We distinguish three different cases :
      - local maxima within Win2/2 distance in forward direction (lm = 1)
      - local maxima within Win2/2 distance in backward direction (lm = 2)
      - local maxima within +-Win2/2 in both forw/back directions (lm = 3)
  */
  /* Select maxima forward path (alike segmentation algorithm) */
  int max = 1;
  for (i = 2; i <= j; i++) {
    if (ref_ft.npo[i] > ref_ft.npo[max] + Win2/2) { /* id the distance between two local
                              maxima (i, max) is greater than Win2
                          keep pos max as a local maxima and
                          increment lm flag */
      lm[max]++;
      max = i;
    } else if (ref_ft.nct[i] > ref_ft.nct[max]) {  /* Else, max is not a local maxima */
      max = i;
    }
  }
  lm[max]++;
  /* Select maxima backward path */
  max = j;
  for (i = j - 1; i >= 0; i--) {
    if (ref_ft.npo[i] < ref_ft.npo[max] - Win2/2) { /* id the distance between two local
                              maxima (i, max) is greater than Win2
                          keep pos max as a local maxima and
                          increment by 2 lm flag */
      lm[max] += 2;
      max = i;
    } else if (ref_ft.nct[i] >= ref_ft.nct[max]) {  /* Else, max is not a local maxima */
      max = i;
    }
  }
  lm[max] += 2;

  /* Print out local maxima positions */
  /* Only positions with lm[i]=3 should be considered as signal peaks */
  for (i = 1; i <= j; i++) {
    if (options.debug) {
      if (ref_ft.ft != NULL) {
        fprintf(stderr,"dbg: %s\t%s\t%d\t%c\t%d\t%d\n", ref_ft.seq_id, ref_ft.ft, ref_ft.npo[i], str, ref_ft.nct[i], lm[i]);
      } else {
        fprintf(stderr,"dbg: %s\t%s\t%d\t%c\t%d\t%d\n", ref_ft.seq_id, ref_ft.name[i], ref_ft.npo[i], str, ref_ft.nct[i], lm[i]);
      }
    }
    if (lm[i] == 3) {
      if (!options.refine) {
        if (ref_ft.ft != NULL) {
          printf("%s\t%s\t%d\t%c\t%d\n", ref_ft.seq_id, ref_ft.ft, ref_ft.npo[i], str, ref_ft.nct[i]);
        } else {
          printf("%s\t%s\t%d\t%c\t%d\n", ref_ft.seq_id, ref_ft.name[i], ref_ft.npo[i], str, ref_ft.nct[i]);
        }
      } else {
        /* Refine peak positions */
    /*
    printf("\ndbg: Refine peak pos for POS %d i=%d : \tk=%d : ref_ft.cnt %d\tref_ft.pos  %d\n", ref_ft.npo[i], i, ref_ft.ptr[i], ref_ft.cnt[ref_ft.ptr[i]], ref_ft.pos[ref_ft.ptr[i]]);
    */
        sum = ref_ft.npo[i] * ref_ft.cnt[ref_ft.ptr[i]];
        cnts = ref_ft.cnt[ref_ft.ptr[i]];
    /*
    printf("\ndbg:  INIT: sum %llu   cnts %d \n\n", sum, cnts);
    */
        /* Recompute peak position within the window range +-(Win1/2) */
        for (k = ref_ft.ptr[i] - 1; ref_ft.pos[k] >= ref_ft.npo[i] - Win1/2; k--) {
      /*
      printf("i=%d ref_ft.npo %d\tk=%d : ref_ft.cnt %d\tref_ft.pos  %d\n", i, ref_ft.npo[i], k, ref_ft.cnt[k], ref_ft.pos[k]);
      */
          sum += ref_ft.cnt[k]*ref_ft.pos[k];
          cnts += ref_ft.cnt[k];
        }
        for (k = ref_ft.ptr[i] + 1; ref_ft.pos[k] <= ref_ft.npo[i] + Win1/2; k++) {
      /*
      printf("i=%d ref_ft.npo %d\tk=%d : ref_ft.cnt %d\tref_ft.pos  %d\n", i, ref_ft.npo[i], k, ref_ft.cnt[k], ref_ft.pos[k]);
      */
          sum += ref_ft.cnt[k]*ref_ft.pos[k];
          cnts += ref_ft.cnt[k];
        }
    /*
    printf("\ndbg: OLD POS %d : \t", ref_ft.npo[i]);
    printf("SUM = %llu CNTS = %d NEW POS : %d \n\n", sum, cnts, (int)sum/cnts);
    */
    ref_ft.npo[i] = sum/cnts;
        if (ref_ft.ft != NULL) {
          printf("%s\t%s\t%d\t%c\t%d\n", ref_ft.seq_id, ref_ft.ft, ref_ft.npo[i], str, ref_ft.nct[i]);
        } else {
          printf("%s\t%s\t%d\t%c\t%d\n", ref_ft.seq_id, ref_ft.name[i], ref_ft.npo[i], str, ref_ft.nct[i]);
        }
      } /* If refine peak pos */
    }
  }
}

int
process_sga(char *iFile)
{
  FILE *f = fopen(iFile, "r");
  char seq_id_prev[SEQ_ID] = "";
  int pos, cnt;
  size_t mLen = BUF_SIZE;
  char *s, *res, *buf;
  size_t bLen = LINE_SIZE;
  unsigned int k = 0;

  if (f == NULL) {
    fprintf(stderr, "Could not open file %s: %s(%d)\n",
            iFile, strerror(errno), errno);
    return 1;
  }
  // `ref_ft` is a `feature_t` struct
  // Allocate memory for all its elements
  if (( ref_ft.pos = (int*)calloc(mLen, sizeof(int))) == NULL) {
    perror("process_sga: malloc");
    exit(1);
  }
  if (( ref_ft.cnt = (int*)calloc(mLen, sizeof(int))) == NULL) {
    perror("process_sga: malloc");
    exit(1);
  }
  if (( ref_ft.npo = (int*)calloc(mLen, sizeof(int))) == NULL) {
    perror("process_sga: malloc");
    exit(1);
  }
  if (( ref_ft.nct = (int*)calloc(mLen, sizeof(int))) == NULL) {
    perror("process_sga: malloc");
    exit(1);
  }
  if (( ref_ft.ptr = (int*)calloc(mLen, sizeof(int))) == NULL) {
    perror("process_sga: malloc");
    exit(1);
  }
  if (ref_ft.ft == NULL) {
    if (( ref_ft.name = (char**)calloc(mLen, sizeof(*(ref_ft.name)))) == NULL) {
      perror("process_sga: malloc");
      exit(1);
    }
  }
  if ((s = malloc(bLen * sizeof(char))) == NULL) {
    perror("process_sga: malloc");
    exit(1);
  }
  if (options.debug)
    fprintf(stderr, "Processing file %s\n", iFile);
#ifdef DEBUG
  int c = 1;
#endif
/*
  // Equivalent to
  while (fscanf(f,"%s %s %d %c %d", seq_id, ft, &pos, &strand, &cnt) != EOF) {
  */
  while ((res = fgets(s, (int) bLen, f)) != NULL) {
    // `s` is either the whole next line, of a string of length `bLen`
    // `res` is useless (?)
    size_t cLen = strlen(s);     // length of the current line
    char seq_id[SEQ_ID] = "";    // chromosome
    char ft[FT_MAX] = "";        // feature name
    char position[POS_MAX] = ""; // position
    char count[CNT_MAX] = "";    // count
    char ext[EXT_MAX] = "";      // ?
    char strand = '\0';          // strand
    unsigned int i = 0;

    // While bLen is not enough to capture the whole line, allocate more memory
    while (cLen + 1 == bLen && s[cLen - 1] != '\n') {
      bLen *= 2;
      if ((s = realloc(s, bLen)) == NULL) {
        perror("process_file: realloc");
        exit(1);
      }
      res = fgets(s + cLen, (int) (bLen - cLen), f);
      cLen = strlen(s);
    }
    // Make the string valid (?)
    if (s[cLen - 1] == '\n')
      s[cLen - 1] = 0;

    buf = s; // same pointer, intended to `strcopy`?
    /* Get SGA fields */
    /* SEQ ID */
    while (*buf != 0 && !isspace(*buf)) {
      // Store in `seq_id` all chars of `s` as long as they are neither space-like or \0
      if (i >= SEQ_ID) {
        fprintf(stderr, "Seq ID is too long \"%s\" \n", buf);
        exit(1);
      }
      seq_id[i++] = *buf++;
    }
    while (isspace(*buf)) // escape empty spaces
      buf++;
    /* FEATURE */
    i = 0;
    while (*buf != 0 && !isspace(*buf)) {
      if (i >= FT_MAX) {
        fprintf(stderr, "Feature is too long \"%s\" \n", buf);
        exit(1);
      }
      ft[i++] = *buf++;
    }
    while (isspace(*buf))
      buf++;
    /* Position */
    i = 0;
    while (isdigit(*buf)) {
      if (i >= POS_MAX) {
        fprintf(stderr, "Position is too large \"%s\" \n", buf);
        exit(1);
      }
      position[i++] = *buf++;
    }
    position[i] = 0;
    pos = atoi(position);
    while (isspace(*buf))
      buf++;
    /* Strand */
    strand = *buf++;
    while (isspace(*buf))
      buf++;
    /* Counts */
    i = 0;
    while (isdigit(*buf)) {
      if (i >= CNT_MAX) {
        fprintf(stderr, "Count is too large \"%s\" \n", buf);
        exit(1);
      }
      count[i++] = *buf++;
    }
    count[i] = 0;
    cnt = atoi(count);
    while (isspace(*buf))
      buf++;
    /* SGA Extension */
    i = 0;
    while (*buf != 0) {
      if (i >= EXT_MAX) {
        fprintf(stderr, "Extension is too long \"%s\" \n", buf);
        exit(1);
      }
      ext[i++] = *buf++;
    }

#ifdef DEBUG
    printf(" [%d] seq ID: %s   Feat: %s (%c)  Pos: %d  Cnts: %d Ext: %s\n", c++, seq_id, ft, strand, pos, cnt, ext);
#endif
    // More f...ing memory allocation tests
    if (k >= mLen - 1) {
      mLen *= 2; // double buffer size if necessary, then reallocate `ref_ft` slots
#ifdef DEBUG
      fprintf(stderr, "reallocating memory for ref_ft.pos ref_ft.strand ref_ft.cnt (k=%d, size=%d)\n", j, mLen);
#endif
      if (( ref_ft.pos = (int *)realloc(ref_ft.pos, mLen * sizeof(int))) == NULL) {
        perror("process_sga: realloc");
        exit(1);
      }
      if (( ref_ft.cnt = (int *)realloc(ref_ft.cnt, mLen * sizeof(int))) == NULL) {
        perror("process_sga: realloc");
        exit(1);
      }
      if (ref_ft.ft == NULL) {
        if (( ref_ft.name = (char**)realloc(ref_ft.name, mLen * sizeof(*(ref_ft.name)))) == NULL) {
          perror("process_sga: malloc");
          exit(1);
        }
      }
    }

    /* Check Chromosome BEGINNING, process previous signal peaks and printout results*/
    // `ref_ft` is a `feature_t` struct
    // If chromosome changes
    if (strcmp(seq_id, seq_id_prev) != 0) {
      ref_ft.pos[0] = ref_ft.pos[1] - Win1/2 - 1;
      ref_ft.pos[k + 1] = ref_ft.pos[k] + Win1/2 + 1;
      //===== MAIN CALL ====//
      locate_peaks((int)k);
      //==== \MAIN CALL ====//
      strcpy(seq_id_prev, seq_id); // the current id becomes the previous
      k = 0;
    }
    // ?
    if (ref_ft.ft == NULL) {
        k++;
        strcpy(ref_ft.seq_id, seq_id);
        ref_ft.name[k] = malloc(strlen(ft) + 1);
        strcpy(ref_ft.name[k], ft);
        ref_ft.pos[k] = pos;
        if (cnt > Coff)
            ref_ft.cnt[k] = Coff;
        else
            ref_ft.cnt[k] = cnt;
    } else if (ref_ft.ft_str == '\0') {
      if (strcmp(ft, ref_ft.ft) == 0) {
        k++;
        strcpy(ref_ft.seq_id, seq_id);
        strcpy(ref_ft.ft, ft);
        ref_ft.pos[k] = pos;
        if (cnt > Coff)
            ref_ft.cnt[k] = Coff;
        else
            ref_ft.cnt[k] = cnt;
        }
    } else if (strand_flag == 1) {
      if (strand == ref_ft.ft_str) {
        k++;
        strcpy(ref_ft.seq_id, seq_id);
        strcpy(ref_ft.ft, ft);
        ref_ft.pos[k] = pos;
        if (cnt > Coff)
            ref_ft.cnt[k] = Coff;
        else
            ref_ft.cnt[k] = cnt;
        }
    } else {
      if (strcmp(ft, ref_ft.ft) == 0  && strand == ref_ft.ft_str) {
        k++;
        strcpy(ref_ft.seq_id, seq_id);
        strcpy(ref_ft.ft, ft);
        ref_ft.pos[k] = pos;
        if (cnt > Coff)
            ref_ft.cnt[k] = Coff;
        else
            ref_ft.cnt[k] = cnt;
      }
    }
  } /* End of While */
  /* Locate signal peaks for the last chromosome */
  ref_ft.pos[0] = ref_ft.pos[1] - Win1/2 - 1;
  ref_ft.pos[k + 1] = ref_ft.pos[k] + Win1/2 + 1;
  //===== MAIN CALL ====//
  locate_peaks((int)k);
  //==== \MAIN CALL ====//
  fclose(f);
  return 0;
}

int
main(int argc, char *argv[])
{
#ifdef DEBUG
  mcheck(NULL);
  mtrace();
#endif

  while (1) {
    int c = getopt(argc, argv, "f:drhw:v:t:c:");
    if (c == -1)
      break;
    switch (c) {
      case 'f':
        Feature = optarg;
        break;
      case 'd':
        options.debug = 1;
        break;
      case 'r':
        options.refine = 1;
        break;
      case 'h':
        options.help = 1;
        break;
    case 'w':
        Win1 = atoi(optarg);
    break;
    case 'v':
        Win2 = atoi(optarg);
    break;
    case 't':
        Thres = atoi(optarg);
    break;
    case 'c':
        Coff = atoi(optarg);
    break;
      default:
        printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }
  if (optind == argc || options.help == 1 || Win1 == 0 || Win2 == 0 ) {
    fprintf(stderr, "Usage: %s [options] -f <feature name> -w <window> -v <vicinity> <SGA File>\n"
             "      where options are:\n"
         "  \t\t -h     Show this message\n"
         "  \t\t -d     Produce debugging output\n"
         "  \t\t -r     Refine Peak Positions\n"
         "  \t\t -c     Count Cut-off (default is %d)\n"
         "  \t\t -t     Peak Threshold (default is %d)\n"
         "\n\tLocates signal peaks within SGA files.\n"
         "\n\tThe program reads a ChIP-seq data file in SGA format (<SGA File>),\n"
         "\tand detects signal peaks for ChIP-tag positions corresponding to\n"
         "\ta specific feature <feature name> if it is given.\n"
         "\tAdditional input parameters are the integration range (<window>),\n"
         "\tthe minimal distance amongst a group of high count values\n"
             "\t(<vicinity>), and the peak threshold (<threshold>).\n"
         "\tA value can be optionally specified as a cut-off for the feature counts.\n"
             "\tThe output is an SGA-formatted list containing signal peak locations\n\n",
         argv[0], Coff, Thres);
      return 1;
  }
  if (options.debug) {
    fprintf(stderr, " Arguments:\n");
    fprintf(stderr, " Selected Feature : %s\n", Feature);
    fprintf(stderr, " Integration range (Window) : %d\n\n", Win1);
    fprintf(stderr, " Minimal distance (Vicinity) : %d\n\n", Win2);
    fprintf(stderr, " Peak Threshold : %d\n\n", Thres);
  }
  /* Process Feature Specs */
  if (Feature == NULL) {
    ref_ft.ft = NULL;      /* Process all features */
    ref_ft.ft_str = '\0';
  } else {
    ref_ft.ft = malloc(FT_MAX * sizeof(char));
    char *s = Feature;
    int i = 0;
    while (*s != 0  && !isspace(*s)) {
      if (i >= FT_MAX) {
        fprintf(stderr, "Feature Description too long \"%s\" \n", Feature);
        return 1;
      }
      ref_ft.ft[i++] = *s++;
    }
    ref_ft.ft_str = '\0';
    while (isspace(*s++))
    ref_ft.ft_str = *s;
  }
  if (options.debug) {
    if (ref_ft.ft_str == '\0' && ref_ft.ft == NULL) {
      fprintf(stderr, "Feature Specs: ALL -> Process all features\n");
    } else if (ref_ft.ft_str == '\0') {
      fprintf(stderr, "Feature Specs: Feature name : %s\n", ref_ft.ft);
    } else {
      fprintf(stderr, "Feature Specs: Feature name/str : %s %c\n", ref_ft.ft, ref_ft.ft_str);
    }
  }
  if ( ref_ft.ft != NULL && (strcmp(ref_ft.ft, "+") == 0 || strcmp(ref_ft.ft, "-") == 0)) {
    strcpy(&ref_ft.ft_str, ref_ft.ft);
    strand_flag = 1;
    if (options.debug)
      fprintf(stderr, "Feature Specs: Process all features on str : %c\n", ref_ft.ft_str);
  }
  if (process_sga(argv[optind++]) != 0) {
    return 1;
  }
  return 0;
}
