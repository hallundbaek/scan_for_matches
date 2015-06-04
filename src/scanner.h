#ifndef SCANNER_H
#define SCANNER_H

struct scan_state {
  punit_t *backtrack_pu;    // BR Back reference
  punit_t *current_pu;      // CR Current punit
  char *fasta_pos;        // SR Current place in data.
  char *fasta_start;     // start
  char *fasta_end;       // end
};

typedef struct scan_state scan_state_t;

typedef struct weight_vec {
  int aw, cw, gw, tw;
} wv_t; //not in parser.c

char *start_srch, *end_srch, *past_last; //not in parser.c

extern int comp_data(char *in, char *out);
extern int first_match(char *start, int len, char *hits[], int show_overlaps, int *matched_punits);
extern int cont_match(char *hits[], int *matched_punits);

#endif
