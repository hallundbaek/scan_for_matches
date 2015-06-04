#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <time.h>
#include "scan_for_matches.h"

#define MAX_SEQ_LEN 250000000
#define MAX_PAT_LINE_LN 32000

/*  TEST program:

    scan_for_matches pattern < fasta_input > fasta_output
*/

/* ============================================================= */
/* ============================================================= */


char
compl (c)
     char c;
{
  switch (c) {
  case 'a':
    return 't';
  case 'A':
    return 'T';

  case 'c':
    return 'g';
  case 'C':
    return 'G';

  case 'g':
    return 'c';
  case 'G':
    return 'C';

  case 't':
  case 'u':
    return 'a';
  case 'T':
  case 'U':
    return 'A';

  case 'm':
    return 'k';
  case 'M':
    return 'K';

  case 'r':
    return 'y';
  case 'R':
    return 'Y';

  case 'w':
    return 'w';
  case 'W':
    return 'W';

  case 's':
    return 'S';
  case 'S':
    return 'S';

  case 'y':
    return 'r';
  case 'Y':
    return 'R';

  case 'k':
    return 'm';
  case 'K':
    return 'M';

  case 'b':
    return 'v';
  case 'B':
    return 'V';

  case 'd':
    return 'h';
  case 'D':
    return 'H';

  case 'h':
    return 'd';
  case 'H':
    return 'D';

  case 'v':
    return 'b';
  case 'V':
    return 'B';

  case 'n':
    return 'n';
  case 'N':
    return 'N';

  default:
    return c;
  }
}

match_t*
find_matches(fasta_file, line, ignore_file, protein, complements, show_overlaps, max_hits, stop_after, verbose)
    char *fasta_file, *line, *ignore_file;
    int protein, max_hits, stop_after, verbose, show_overlaps, complements;
{
  int hit_in_line;
  char tmp;
  char* ignore[20000];
  int ig_index = 0;
  int i, j, k, i1;
  char *hits[2000];
  char id[1000];
  char *p, *pc, *mp;
  match_t *first, *match;
  int ln;
  int got_gt;
  int is_first_match = 1;
  int is_non_empty_match = 0;
  char *data, *cdata;
  FILE *input;
  FILE *ig_fp;
  if (ig_fp = fopen(ignore_file, "r")) {
    while (fscanf (ig_fp, "%s", id) == 1) {
      if ((ignore[ig_index] = malloc (strlen (id) + 1)) == NULL) {
        fprintf (stderr, "memory allocation error\n");
        exit (1);
      }
      strcpy (ignore[ig_index++], id);
    }
    close (ig_fp);
    if (ig_index) {
      fprintf (stderr, "ignoring %d id(s)\n", ig_index);
    }
  }
  if (((data = malloc (MAX_SEQ_LEN + 1)) == NULL) ||
      ((cdata = malloc (MAX_SEQ_LEN + 1)) == NULL)) {
    fprintf (stderr, "failed to alloc memory\n");
  }
  first = malloc(sizeof(match_t));
  if ((protein && !parse_peptide_cmd (line))
      || (!protein && !parse_dna_cmd (line))) {
    fprintf (stderr, "failed to parse pattern: %s\n", line);
    exit (1);
  }
  if (!(input = fopen(fasta_file, "r"))) {
    input = stdin;
  }
  while ((max_hits > 0) && ((!got_gt && (fscanf (input, ">%s", id) == 1)) ||
         (got_gt && (fscanf (input, "%s", id) == 1)))) {
    while (getc (input) != '\n');
    for (p = data; ((i = getc (input)) != -1) && (i != '>');) {
      if ((i != ' ') && (i != '\n')) {
        *(p++) = i;
      }
    }
    if (i == '>')
      got_gt = 1;
    else
      got_gt = 0;

    *p = 0;

    for (i = 0; (i < ig_index) && (strcmp (ignore[i], id) != 0); i++);

    if (i == ig_index) {
      if (!protein) {
        comp_data (data, cdata);
      }
      else {
        strcpy (cdata, data);
      }

      ln = strlen (data);
      for (hit_in_line = 0, first_match (cdata, ln, hits, show_overlaps, &i);
           (max_hits > 0) && (i > 0); )
      {
        hit_in_line = 1;
        max_hits--;
        if (is_first_match) {
          match = first;
          is_first_match = 0;
        } else {
          match->next = (match_t*) malloc( sizeof( match_t ));
          match = match->next;
        }
        match->start = 1 + hits[0] - cdata;
        match->end   =  hits[i] - cdata;
        match->match = (char*) malloc(sizeof(char)*(match->end - match->start + i + 3));
        mp = match->match;
        for (i1 = 0; i1 < i; i1++) {
          j = hits[i1 + 1] - hits[i1];
          is_non_empty_match = 0;
          for (p = data + (hits[i1] - cdata); j; j--) {
            *(mp++) = *(p++);
            is_non_empty_match = 1;
          }
          if (is_non_empty_match || i1 + 1 == i) {
            *(mp++) = ' ';
          }
        }
        *(mp--) = '\0';
        if (verbose){
          printf (">%s:[%ld,%ld]\n", id, match->start, match->end);
          printf ( "%s\n", match->match );
        }
        if (!show_overlaps) {
          cont_match (hits, &i);
        } else {
          next_match (hits, &i);
        }
      }
      if (complements) {
        for (i = 0, j = ln - 1; i <= j; i++, j--) {
          tmp = compl (data[i]);
          data[i] = compl (data[j]);
          data[j] = tmp;
        }
        comp_data (data, cdata);
        for (first_match (cdata, ln, hits, show_overlaps, &i);
             (max_hits > 0) && (i > 0); cont_match (hits, &i))
        {
          hit_in_line = 1;
          max_hits--;
          if (is_first_match) {
            match = first;
            is_first_match = 0;
          } else {
            match->next = (match_t*) malloc( sizeof( match_t ));
            match = match->next;
          }
          match->start = 1 + hits[0] - cdata;
          match->end   =  hits[i] - cdata;
          match->match = (char*) malloc(sizeof(char)*(match->end - match->start + i + 3));
          mp = match->match;
          for (i1 = 0; i1 < i; i1++) {
            j = hits[i1 + 1] - hits[i1];
            is_non_empty_match = 0;
            for (p = data + (hits[i1] - cdata); j; j--) {
              *(mp++) = *(p++);
              is_non_empty_match = 1;
            }
            if (is_non_empty_match || i1 + 1 == i) {
              *(mp++) = ' ';
            }
          }
          *(--mp) = '\0';
          if (verbose){
            printf (">%s:[%ld,%ld]\n", id, match->start, match->end);
            printf ( "%s\n", match->match );
          }
        }
        if (!hit_in_line) {
          if (--stop_after == 0) {
            fprintf (stderr, "exceeded limit of lines failing to match\n");
            fclose(input);
            exit (1);
          }
        }
      }
    }
  }
  fclose(input);
  if (!is_first_match)
    match->next = NULL;
  return first;
}

int
main (argc, argv)
     int argc;
     char *argv[];
{
  extern char *optarg;
  extern int optind;

  char *data, *cdata;
  char line[MAX_PAT_LINE_LN];
  int got_gt = 0;
  FILE *ig_fp = NULL;
  FILE *fp;
  int i;
  char *hits[2000];
  char *pc;
  char id[1000];
  char fasta_file[1000];
  char ignore_file[1000];

  int stop_after = 100000000;
  int max_hits = 100000000;
  int command_line_pattern = 0;
  int show_overlaps = 0;
  int complements = 0;
  int protein = 0;
  int errflag = 0;
  int verbose = 1;
  int c;

  clock_t start_t, end_t;
  double total_time;
  /* Start the clock */
  start_t = clock();

    while ((c = getopt (argc, argv, "spcnmo:i:P:f:")) != -1)
    switch (c) {

    case 'o':
      show_overlaps = 1;
      break;

    case 'c':
      complements = 1;
      break;

    case 'p':
      protein = 1;
      break;

    case 'n':
      if (sscanf (optarg, "%d", &stop_after) != 1) {
        errflag = 1;
        fprintf (stderr, "invalid value on -n option (make it a positive integer)\n");
      }
      break;

    case 'm':
      if (sscanf (optarg, "%d", &max_hits) != 1) {
        errflag = 1;
        fprintf (stderr, "invalid value on -m option (make it a positive integer)\n");
      }
      break;

    case 'i':
      if ((ig_fp = fopen (optarg, "r")) == NULL) {
        errflag = 1;
        fprintf (stderr, "invalid file name for ids to ignore\n");
      }
      break;
    case 'P':
      if (sscanf(optarg, "%[^\t\n]", line) < 1) {
        errflag = 1;
        fprintf(stderr, "invalid pattern line\n");
      } else {
        command_line_pattern = 1;
      }
      break;
    case 'f':
      if ((sscanf(optarg, "%[^\t\n]", fasta_file) < 1 ) || !(fp = fopen(fasta_file, "r"))) {
        errflag = 1;
        fprintf(stderr, "invalid fasta file\n");
      } else {
        fclose(fp);
      }
      break;
    case 's':
      verbose = 0;
      break;
    }
  if (errflag || !command_line_pattern && ((optind >= argc) || ((fp = fopen (argv[optind], "r")) == NULL))) {
    fprintf (stderr, "errflag=%d optind=%d argc=%d\n", errflag, optind, argc);
    fprintf (stderr, "usage: scan_for_matches -c [for complementary strand] -p [protein][-n stop_after_n_misses] [-m max_hits] [-i file_of_ids_to_ignore] -o [overlapping hits] pattern_file < fasta_input > hits\n");
    exit (2);
  }

  /* Replaces comments and newlines with ' ' and
   * place pattern in the variable line.
   */
  if (!command_line_pattern) {
    pc = line;
    i = fgetc (fp);
    while ((i != EOF) && (pc < line + (MAX_PAT_LINE_LN - 1))) {
      if (i == '%') {
        while (((i = fgetc (fp)) != '\n') && (i != EOF));
        if (i != EOF) {
          i = ' ';
        }
      }
      else if (i == '\n') {
        i = ' ';
      }
      else {
        *(pc++) = i;
        i = fgetc (fp);
      }
    }
    *pc = 0;
    close (fp);
  }
  find_matches(fasta_file, line, ignore_file, protein, complements, show_overlaps, max_hits, stop_after, verbose);
  /*    printf("successfully completed\n"); */
  end_t = clock();
  total_time = (double)(end_t - start_t) / CLOCKS_PER_SEC;
//  printf("Total time taken by CPU:\n%f\n", total_time  );

  return (EXIT_SUCCESS);
}


