#ifndef PARSER_H
#define PARSER_H

#define PEPTIDE 1
#define DNA 2

#define BOOL  int
#define TRUE  1
#define FALSE 0

#define MAX_SOUGHT_CHARS 1000000
#define MAX_PUNITS 100
#define MAX_CODES 100    /* max length of a single punit being matched to */ 
#define MAX_WEIGHTS 2000  /* max weight storage */ 
#define MAX_TOTAL_CODES 600  /* max character info storage */
#define MAX_NAMES    50

#define EXACT_PUNIT    0  /* punit_t.type */
#define RANGE_PUNIT    1  /* punit_t.type */
#define COMPL_PUNIT    2  /* punit_t.type */
#define REPEAT_PUNIT   3  /* punit_t.type */
#define SIM_PUNIT      4  /* punit_t.type */
#define WEIGHT_PUNIT   5  /* punit_t.type */
#define OR_PUNIT       6  /* punit_t.type */
#define ANY_PUNIT      8  /* punit_t.type */
#define LLIM_PUNIT     9  /* length limit punit */
#define INV_REP_PUNIT 10  /* for real palindromic inverted repeats */
#define MATCH_START   11  /* ^ */
#define MATCH_END     12  /* $ */

#define A_BIT 0x01    /* char bitfield */
#define C_BIT 0x02    /* char bitfield */
#define G_BIT 0x04    /* char bitfield */
#define T_BIT 0x08    /* char bitfield */

#define MAX(A,B) (((A))>((B)) ? ((A)) : ((B)))

char punit_to_code[256]; 
char code_to_punit[256]; 

struct punit {
  int type;
  struct punit *nxt_punit, *prev_punit;
  struct punit *BR;
  int anchored;
  char *hit;
  int mlen;      /* how much matched */
  union {
    struct or_info_struct {
      struct punit *or1, *or2;
      char *SR;
      int alt;
    } or;
    struct exact_info_struct {
      int len;
      char *code;
    } exact;
    struct any_struct {
      long code_matrix;
    } any;
    struct char_info_struct {
      int len;
      char *code;
    } char_match;
    struct range_info_struct {
      int min, width;    /* width of x...y is y-x */
      int nxt;
    } range;
    struct compl_info_struct {
      int ins, del, mis, of, rule_set;
      int ins_left;
    } compl;
    struct repeat_info_struct {
      int ins, del, mis, of, nxtent;
      int ins_left;
    } repeat;
    struct sim_info_struct {
      int ins, del, mis, len, nxtent;
      int ins_left;
      char *code;
    } sim;
    struct llim_struct {
      int llim_vec[MAX_NAMES + 1];
      int bound;
    } llim;
    struct weight_info_struct {
      int len;
      int *vec;
      int cutoff;
      int tupsz;
      int maxwt;
    } wvec;
  } info;
};

typedef struct punit punit_t;

/* state is used to retain state for a search -- needed for parallel
   Prolog versions
*/
struct state {
  punit_t pu_s[MAX_PUNITS];
  char cv[MAX_TOTAL_CODES];
  int iv[21 * MAX_WEIGHTS / 2];
  punit_t *BR;
  punit_t *ad_pu_s;
  char *past_last;
};

struct punit *names[MAX_NAMES];

punit_t *BR1;

char *rule_sets[MAX_NAMES];


extern int parse_dna_cmd(char *line);
extern int parse_peptide_cmd(char *line);
extern int find_pivot_punit(int *max_before_length, int *max_after_length, int *index);

#endif
