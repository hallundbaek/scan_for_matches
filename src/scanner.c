/*  NOTES: June 30, 1992

    These result from considering what changes will need to be made
    to allow this to be effectively called from OR-parallel Prolog systems.


    1. The comments on what parse_dna_cmd and parse_peptide_cmd returned
       were wrong (they said 1, but it returns maximum # of matched chars).
       Why return max # matched?

    2. The routines should "auto-initialize".  This is probably not done,
       due to GenoGraphics and IGD requiring slighly differing 
       initializations.  These should be made consistent (code_to_punit goes
       into ggpunit.c

    3. Forms of parse_cmd, parse_dna_cmd, and parse_peptide_cmd should be
       written to return the PUs (and BR) in dynamically allocated memory.

    4. Forms of first_match, next_match, and pattern_match need to be
       written that operate on the dynamically allocated state.

For better or worse, I did these changes. 6/30/92 RAO
*/

/* ggpunit.c contains routines to parse a string and store it as a punit
   and to search encoded uncompressed DNA or uncoded protein sequence
*/
/* Style:
   Every function must have a prototype inside the #ifndef CC
   Where it is defined, each function must use the old style for
   typing parameters, i.e.:

char *char_pat(n,p)
int *n;
char *p;
{
   ....
}

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "parser.h"
#include "scanner.h"

#define MAX(A,B) (((A))>((B)) ? ((A)) : ((B)))

FILE *diag_file;

/* the following variable hides the relationship between parsing punits
   and executing searches.  These routines compose the interface:

   int parse_dna_cmd(line)      returns max # characters that could be matched
                                    if successful parse (sets ad_pu_s)
   char *line;                  returns 0 if unsuccessful parse

   int parse_peptide_cmd(line)  returns max # characters that could be matched
                                    if successful parse (sets ad_pu_s)
   char *line;                  returns 0 if unsuccessful parse

   int first_match(start,len,hits)  0 on miss; # pattern units on hit
   char *start;             start of search
   int len;                 length of search
   char *hits[];            on hit, set to vector of starts of punits +
                                just past the end of last punit, so if
                                the routine returns N, hits will have N+1
                                entries.  Each entry will the point where
                                a punit matched, except the last, which
                                is just past the last punit

   int next_match(hits)      continues search
                            (returns 0 or # pattern units to indicate success)

   int cont_match(hits)      continues search (non-overlapping)
                            (returns 0 or # pattern units to indicate success)

   char *hits[];              set on success
*/

char known_char[16] = { 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 };
char known_char_index[16] =
  { -1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1 };

int punit_sequence_type;


punit_t middle_punit;
punit_t start_punit = { .type = MATCH_START };
punit_t end_punit = { .type = MATCH_END };
int max_before_mlen;
int max_after_mlen;
int best_index;

extern punit_t *ad_pu_s;

/* Function prototypes for typechecking at compile time  */
#ifdef CC
#define PROTO(X) ()
#else
#define PROTO(X) X
#endif

void rev_compl_data
PROTO ((int rule_set, unsigned char *data, int len, unsigned char *result));
int loose_match
PROTO ((unsigned char *one_data, int one_len, unsigned char *two_data,
  int two_len, int max_ins, int max_del, int max_mis, int rule_set,
  int *match_range, int compl_flag, int *ins_left));
int max_mats PROTO ((struct punit * pu));
int cont_match PROTO ((char *hits[], int *matched_punits));
int pattern_match
PROTO ((punit_t * pu, char *start, char *end, char *hits[], int first,
  punit_t ** BR1));
int collect_hits PROTO ((struct punit * pu, char **revhits));
int punit_list_split PROTO ((punit_t *first_punit, int *max_before_len, int * max_after_len));
int pattern_optimized_match PROTO ((punit_t *first_punit, char **start, char *end, char *hits[], int first, punit_t **BR1));
/**
 * @brief first_match searches for a match from the given start param to the end of start + len
 *        and stores the matched pattern units, if any, in hits.
 *
 * @param start Pointer to start of search
 * @param len Length of the search
 * @param hits Contains a list of pointers to where the punits matched. hits will hold 
 *             N+1 entries if the routine returns N, where entry N+1 is just past the 
 *             last punit.
 *
 * @return Returns the number of pattern units matched. 0 on miss
 */
int
first_match (start, len, hits, show_overlaps, matched_punits)
     char *start;
     int len;
     char *hits[];
     int show_overlaps;
     int *matched_punits;
{
  if (ad_pu_s == NULL)
    return -1;
  start_srch = start;
  end_srch = start + (len - 1);
  BR1 = NULL;
  if (!show_overlaps) {
    punit_list_split(ad_pu_s, &max_before_mlen, &max_after_mlen);
  } else {
    best_index = 0;
  }
  if (best_index) {
    *matched_punits = pattern_optimized_match (ad_pu_s, &start_srch, end_srch, hits, 1 , &BR1);
  } else {
    *matched_punits = pattern_match (ad_pu_s, start_srch, end_srch, hits, 1, &BR1);
  }
  past_last = hits[*matched_punits];
  return *matched_punits;
}

/**
 * @brief Continues overlapping search
 *
 * @param hits Contains a list of pointers to where the punits matched. hits will hold
 *             N+1 entries if the routine returns N, where entry N+1 is just past the
 *             last punit.
 *
 * @return Returns the number of pattern units matched. 0 on miss
 */
int
next_match (hits, matched_punits)
     char *hits[];
     int *matched_punits;
{
  if (ad_pu_s == NULL)
    return -1;
  if (best_index) {
    *matched_punits = pattern_optimized_match (ad_pu_s, &start_srch, end_srch, hits, 0, &BR1);
  } else {
    *matched_punits = pattern_match (ad_pu_s, start_srch, end_srch, hits, 0, &BR1);
  }
  past_last = hits[*matched_punits];
  return (*matched_punits);
}

/**
 * @brief Continues non-overlapping search
 *
 * @param hits Contains a list of pointers to where the punits matched. hits will hold
 *             N+1 entries if the routine returns N, where entry N+1 is just past the
 *             last punit.
 *
 * @return Returns the number of pattern units matched. 0 on miss
 */
int
cont_match (hits, matched_punits)
     char *hits[];
     int *matched_punits;
{
  char *past_last1;
  past_last1 = past_last;
  if (best_index) {
    next_match (hits, matched_punits);
  } else {
    while ((next_match (hits, matched_punits) > 0) && (hits[0] < past_last1));
  }
  past_last = hits[*matched_punits];
  return (*matched_punits);
}

/**
 * @brief Used for OR-parallel Prolog systems (NEEDS MORE INFO)
 *
 * @param hits Contains a list of pointers to where the punits matched. hits will hold
 *             N+1 entries if the routine returns N, where entry N+1 is just past the
 *             last punit.
 * @param st Dynamically allocated state
 *
 * @return Returns the number of pattern units matched. 0 on miss
 */
int
cont_match_ns (hits, st)
     char *hits[];
     struct state *st;
{
  int i;
  char *past_last1;
  past_last1 = st->past_last;
  while (((i = next_match_ns (hits, st)) > 0) && (hits[0] < past_last1));
  st->past_last = hits[i];
  return (i);
}

/**
 * @brief Used for OR-parallel Prolog systems (NEEDS MORE INFO)
 *
 * @param start Pointer to start of search
 * @param len Length of the search
 * @param hits Contains a list of pointers to where the punits matched. hits will hold 
 *             N+1 entries if the routine returns N, where entry N+1 is just past the 
 *             last punit.
 * @param st Dynamically allocated state
 *
 * @return Returns the number of pattern units matched. 0 on miss
 */
int
first_match_ns (start, len, hits, st)
     char *start;
     int len;
     char *hits[];
     struct state *st;
{
  int i;

  if (st->ad_pu_s == NULL)
    return (-1);

  start_srch = start;
  end_srch = start + (len - 1);
  st->BR = NULL;
  i = pattern_match (st->ad_pu_s, start_srch, end_srch, hits, 1, &(st->BR));
  st->past_last = hits[i];
  if (!i)
    free_state (st);
  return i;
}

/**
 * @brief Used for OR-parallel Prolog systems (NEEDS MORE INFO)
 *
 * @param hits Contains a list of pointers to where the punits matched. hits will hold 
 *             N+1 entries if the routine returns N, where entry N+1 is just past the 
 *             last punit.
 * @param st Dynamically allocated state
 *
 * @return Returns the number of pattern units matched. 0 on miss
 */
int
next_match_ns (hits, st)
     char *hits[];
     struct state *st;
{
  int i;

  if (st->ad_pu_s == NULL)
    return (-1);
  i = pattern_match (st->ad_pu_s, start_srch, end_srch, hits, 0, &(st->BR));
  st->past_last = hits[i];
  if (!i)
    free_state (st);
  return i;
}

/**
 * @brief Finds the successor for a given punit if it has one. If pu is e.g. 
 *        in an OR-punit, it will backtrack out of the 'sub-tree' and find
 *        the next punit of the OR-punit.
 *
 * @param pu punit to find the successor of
 *
 * @return The next punit of pu. Returns NULL if no such punit exists
 */
struct punit *
next_punit (pu)
     struct punit *pu;
{
  struct punit *pu1, *pu2;

  if (pu->nxt_punit)
    return pu->nxt_punit;
  for (pu1 = pu, pu2 = pu1->prev_punit;
       (pu2 && ((pu2->nxt_punit == pu1) || (!pu2->nxt_punit)));
        pu1 = pu2, pu2 = pu1->prev_punit);
  if (pu2 && (pu2->nxt_punit != pu1))
    return pu2->nxt_punit;

  return NULL;
}

/* =================================================================== */
#define KnownChar(C)  (punit_sequence_type == PEPTIDE ? 1 : known_char[(C)])

#define Matches(C1,C2) (punit_sequence_type == PEPTIDE ? \
      ((C1 == C2) || (C2 == 'X')): \
                        (KnownChar((C1) & 15) && ((((C1) & 15) & ((C2) & 15)) == (C1 & 15))))
#define MatchRuleSet(RuleSet,C1,C2) (*(rule_sets[RuleSet] + C2 + known_char_index[(C1 & 15)]))
/*
int MatchRuleSet(RuleSet,C1,C2)
int RuleSet;
char C1,C2;
{
   int x,y;

   y = known_char_index[C1 & 15];

   if (y == -1)  return 0;

   x = *(rule_sets[RuleSet]+C2+y);
   return x;
}
*/

#define ExMatches(RuleSet,C1,C2) ((RuleSet == -1) ? \
                                  (Matches(C1,C2)) : \
          (MatchRuleSet(RuleSet,C1,C2)))


/**
 * @brief Returns the pointer to the first hit of the given punit
 *
 * @param pu punit whose hit poistion we want to return
 *
 * @return Pointer to the hit position of the punit.
 */
char *
position_first_match (pu)
     punit_t *pu;
{
  if (pu->type != OR_PUNIT) {
    return pu->hit;
  }
  else {
    if (pu->info.or.alt == 1) {
      return position_first_match (pu->info.or.or1);
    }
    else {
      return position_first_match (pu->info.or.or2);
    }
  }
}

/**
 * @brief Returns the length of the given punits' match
 *
 * @param pu punit whose hit length we want to return
 *
 * @return The match length of a given punit
 */
int
length_of_match (pu)
     punit_t *pu;
{
  punit_t *pu1;

  if (pu->type == OR_PUNIT) {
    if (pu->info.or.alt == 1) {
      pu1 = pu->info.or.or1;
    }
    else {
      pu1 = pu->info.or.or2;
    }
    for (pu->mlen = 0; pu1; pu1 = pu1->nxt_punit)
      pu->mlen += length_of_match (pu1);
  }
  return pu->mlen;
}

/**
 * @brief Will try to match a pattern to a sequence.
 *
 * @param pu The first pattern unit of the pattern.
 * @param start The pointer first base of the sequence.
 * @param end The pointer last base of the sequence.
 * @param hits Will hold the reported hits upon execution
 * @param first Whether or not it is the first match.
 * @param backtrack_punit The backtracking pattern unit.
 *
 * @return How many pattern units matched.
 */
int
pattern_match (pu, start, end, hits, first, backtrack_punit)
     punit_t **backtrack_punit;
     punit_t *pu;
     char *start, *end;
     char *hits[];
     int first;
{
  punit_t *tmp_backtrack_punit;
  int i, j;
  char *revhits[MAX_PUNITS];
  punit_t *current;
  char *seq_start, *seq_end;
  char *last;
  char *p1, *p2, *p3;
  int *pv1, *pv3;
  wv_t *pv2;
  int ln;
  int wval;
  struct punit *pu1;
  char scratch[4000];
  int first_backtrack = first;
  //flags
  int success, try, try_again, backtrack, backtrack_again, leave;
  try = 1;
  success = backtrack = backtrack_again = leave = try_again = 0;

  tmp_backtrack_punit = *backtrack_punit;
  seq_start = start; // First character of data.
  seq_end = end;   // Last character of data.

  if (!first){
    backtrack = 1;
    try = 0;
  }

  current = pu;
  int count = 0;
  while(!leave){
    //TRY segment
    if(try){
      try_again = 0;

      switch (current->type) {
        case MATCH_START:
          if (seq_start == start) {
            current->hit = seq_start;
            success = 1;
            break;
          }
          backtrack = 1;
          break;

        case MATCH_END:
          if (seq_start == end + 1) {
            current->hit = seq_start;
            success = 1;
            break;
          }
          backtrack = 1;
          break;

        case ANY_PUNIT:
          last = seq_end;
          if ((last > seq_start) && current->anchored)
            last = seq_start;
          while (seq_start <= last) {
            i = *seq_start;
            if (i >= 'A' && i <= 'Z' && ((1 << (i - 'A')) & current->info.any.code_matrix)) {
              break;
            }
            seq_start++;
          }
          if (seq_start > last) {
            backtrack = 1;
            break;
          }
          else {
            if ((current->hit = seq_start) < last) {
              current->BR = tmp_backtrack_punit;
              tmp_backtrack_punit = current;
            }
            seq_start++;
            success = 1;
            break;
          }
          break;

        case LLIM_PUNIT:
          for (ln = 0, i = current->info.llim.llim_vec[0]; i;) {
            pu1 = names[current->info.llim.llim_vec[i--]];
            ln += pu1->mlen;
          }
          if (ln < current->info.llim.bound) {
            current->hit = seq_start;
            success = 1;
            break;
          }
          backtrack = 1;
          break;

        case RANGE_PUNIT:
          if ((seq_start + (i = current->info.range.min - 1)) <= seq_end) {
            current->hit = seq_start;
            seq_start += i + 1;
            if (((seq_start <= seq_end) && (current->info.range.width)) ||
                        (!current->anchored && (seq_start <= seq_end))) {
              current->BR = tmp_backtrack_punit;
              tmp_backtrack_punit = current;
              current->info.range.nxt = current->info.range.min + 1;
            }
            success = 1;
            break;
          }
          else {
            backtrack = 1;
            break;
          }
          break;

        case EXACT_PUNIT:
          last = seq_end + 1 - current->info.exact.len;
          if ((last > seq_start) && current->anchored)
            last = seq_start;
          p1 = current->info.exact.code;
          ln = current->info.exact.len - 1;
          while (seq_start <= last) {
            if (Matches (*seq_start, *p1)) {
              p2 = seq_start + 1;
              p3 = p1 + 1;
              for (i = ln; i && Matches (*p2, *p3); i--, p3++, p2++);
              if (!i)
                break;
            }
            seq_start++;
          }
          if (seq_start > last) {
            backtrack = 1;
            break;
          }
          else {
            if ((current->hit = seq_start) < last) {
              current->BR = tmp_backtrack_punit;
              tmp_backtrack_punit = current;
            }
            seq_start += current->info.exact.len;
            success = 1;
            break;
          }
          break;

        case COMPL_PUNIT:
          pu1 = names[current->info.compl.of];
          p1 = pu1->hit;
          ln = pu1->mlen;
          if ((current->info.compl.rule_set != -1) || current->info.compl.ins ||
                           current->info.compl.del || current->info.compl.mis) {
            if (i = loose_match ((unsigned char *) p1, ln,
               (unsigned char *) seq_start, seq_end + 1 - seq_start,
               current->info.compl.ins,
               current->info.compl.del,
               current->info.compl.mis,
               current->info.compl.rule_set, &j, 1,
               &(current->info.compl.ins_left))) {
              i--;
              current->hit = seq_start;
              if ((current->hit = seq_start) < last) {
                current->BR = tmp_backtrack_punit;
                tmp_backtrack_punit = current;
              }
              seq_start += i;
              success = 1;
              break;
            }
          }
          else if (seq_end - seq_start >= ln - 1) {
            current->hit = seq_start;
            p1 = pu1->hit + (ln - 1);
            while (ln--) {
              //* FIX: 4/25/92: complement cannot match ambiguous char *
              if (!KnownChar ((*p1) & 15) || (((*(p1--) >> 4) & 15) != (*(seq_start++) & 15))){
                backtrack = 1;
                break;
              }
            }
            if(backtrack)
              break;
            success = 1;
            break;
          }
          backtrack = 1;
          break;

        case REPEAT_PUNIT:
          pu1 = names[current->info.repeat.of];
          p1 = pu1->hit;

          if (!(ln = pu1->mlen)) {
            current->hit = seq_start;
            success = 1;
            break;
          }

          if (i = loose_match ((unsigned char *) p1, ln,
             (unsigned char *) seq_start, seq_end + 1 - seq_start,
             current->info.repeat.ins,
             current->info.repeat.del,
             current->info.repeat.mis, -1, &j, 0,
             &(current->info.repeat.ins_left))) {
            i--;
            current->hit = seq_start;
            if ((current->hit = seq_start) < last) {
              current->BR = tmp_backtrack_punit;
              tmp_backtrack_punit = current;
            }
            seq_start += i;
            success = 1;
            break;
          }
          backtrack = 1;
          break;

        case INV_REP_PUNIT:
          pu1 = names[current->info.repeat.of];
          p1 = pu1->hit;

          if (!(ln = pu1->mlen)) {
            current->hit = seq_start;
            success = 1;
            break;
          }

          if (ln > 4000) {
            p3 = (char *) malloc (ln);
          }
          else {
            p3 = scratch;
          }
          for (i = 0, j = ln - 1; i < ln; i++, j--) {
            *(p3 + i) = *(p1 + j);
          }

          if (i = loose_match ((unsigned char *) p3, ln,
             (unsigned char *) seq_start, seq_end + 1 - seq_start,
             current->info.repeat.ins,
             current->info.repeat.del,
             current->info.repeat.mis, -1, &j, 0,
             &(current->info.repeat.ins_left))) {
            i--;
            current->hit = seq_start;
            if ((current->hit = seq_start) < last) {
              current->BR = tmp_backtrack_punit;
              tmp_backtrack_punit = current;
            }
            seq_start += i;
            success = 1;
            break;
          }
          backtrack = 1;
          break;

        case SIM_PUNIT:
          last = seq_end + 1 + current->info.sim.ins - current->info.sim.len;
          if ((last > seq_start) && current->anchored)
            last = seq_start;
          p1 = current->info.sim.code;
          ln = current->info.sim.len - 1;
          while (seq_start <= last) {
            if (i = loose_match ((unsigned char *) current->info.sim.code,
               current->info.sim.len,
               (unsigned char *) seq_start, seq_end + 1 - seq_start,
               current->info.sim.ins,
               current->info.sim.del, current->info.sim.mis, -1, &j, 0,
               &(current->info.sim.ins_left))) {
              i--;
              if ((current->hit = seq_start) < last) {
                current->BR = tmp_backtrack_punit;
                tmp_backtrack_punit = current;
              }
              seq_start += i;
              success=1;
              break;
            }
            seq_start++;
          }
          if(success)
            break;
          backtrack = 1;
          break;

        case WEIGHT_PUNIT:
          last = seq_end + 1 - current->info.wvec.len;
          if ((last > seq_start) && current->anchored)
            last = seq_start;
          pv1 = current->info.wvec.vec;
          ln = current->info.wvec.len;
          while (seq_start <= last) {
            if (current->info.wvec.tupsz == 4) {
              for (pv2 = (wv_t *) pv1, p1 = seq_start, i = ln, wval = 0; i; i--, pv2++) {
                switch (*(p1++) & 15) {
                case A_BIT:
                  wval += pv2->aw;
                  break;
                case C_BIT:
                  wval += pv2->cw;
                  break;
                case G_BIT:
                  wval += pv2->gw;
                  break;
                case T_BIT:
                  wval += pv2->tw;
                  break;
                case (A_BIT + C_BIT):
                  wval += (pv2->aw >> 1) + (pv2->cw >> 1);
                  break;
                case (A_BIT + G_BIT):
                  wval += (pv2->aw >> 1) + (pv2->gw >> 1);
                  break;
                case (A_BIT + T_BIT):
                  wval += (pv2->aw >> 1) + (pv2->tw >> 1);
                  break;
                case (C_BIT + G_BIT):
                  wval += (pv2->cw >> 1) + (pv2->gw >> 1);
                  break;
                case (C_BIT + T_BIT):
                  wval += (pv2->cw >> 1) + (pv2->tw >> 1);
                  break;
                case (G_BIT + T_BIT):
                  wval += (pv2->gw >> 1) + (pv2->tw >> 1);
                  break;
                case (C_BIT + G_BIT + T_BIT):
                  wval += (pv2->cw / 3) + (pv2->gw / 3) + (pv2->tw / 3);
                  break;
                case (A_BIT + G_BIT + T_BIT):
                  wval += (pv2->aw / 3) + (pv2->gw / 3) + (pv2->tw / 3);
                  break;
                case (A_BIT + C_BIT + T_BIT):
                  wval += (pv2->aw / 3) + (pv2->cw / 3) + (pv2->tw / 3);
                  break;
                case (A_BIT + C_BIT + G_BIT):
                  wval += (pv2->aw / 3) + (pv2->cw / 3) + (pv2->gw / 3);
                  break;
                case (A_BIT + C_BIT + G_BIT + T_BIT):
                  wval += (pv2->aw >> 2) +
                    (pv2->cw >> 2) + (pv2->gw >> 2) + (pv2->tw >> 2);
                  break;
                }
              }
            }
            else {  //    * handling 20-tuples  (or 21-tuples) *
              for (pv3 = pv1, p1 = seq_start, i = ln, wval = 0;
                   i; i--, pv3 += current->info.wvec.tupsz) {
                switch (*(p1++)) {
                case 'A':
                  wval += *(pv3 + 0);
                  break;
                case 'C':
                  wval += *(pv3 + 1);
                  break;
                case 'D':
                  wval += *(pv3 + 2);
                  break;
                case 'E':
                  wval += *(pv3 + 3);
                  break;
                case 'F':
                  wval += *(pv3 + 4);
                  break;
                case 'G':
                  wval += *(pv3 + 5);
                  break;
                case 'H':
                  wval += *(pv3 + 6);
                  break;
                case 'I':
                  wval += *(pv3 + 7);
                  break;
                case 'K':
                  wval += *(pv3 + 8);
                  break;
                case 'L':
                  wval += *(pv3 + 9);
                  break;
                case 'M':
                  wval += *(pv3 + 10);
                  break;
                case 'N':
                  wval += *(pv3 + 11);
                  break;
                case 'P':
                  wval += *(pv3 + 12);
                  break;
                case 'Q':
                  wval += *(pv3 + 13);
                  break;
                case 'R':
                  wval += *(pv3 + 14);
                  break;
                case 'S':
                  wval += *(pv3 + 15);
                  break;
                case 'T':
                  wval += *(pv3 + 16);
                  break;
                case 'V':
                  wval += *(pv3 + 17);
                  break;
                case 'W':
                  wval += *(pv3 + 18);
                  break;
                case 'Y':
                  wval += *(pv3 + 19);
                  break;
                default:
                  if (current->info.wvec.tupsz > 20)
                    wval += *(pv3 + 20);
                  break;
                }
              }
            }
            if ((wval > current->info.wvec.cutoff) && (wval < current->info.wvec.maxwt)) {
              break;
            }
            seq_start++;
          }
          if (seq_start > last) {
            backtrack = 1;
            break;
          } else {
            if ((current->hit = seq_start) < last) {
              current->BR = tmp_backtrack_punit;
              tmp_backtrack_punit = current;
            }
            seq_start += current->info.wvec.len;
            success = 1;
            break;
          }
          break;

        case OR_PUNIT:
          current->BR = tmp_backtrack_punit;
          tmp_backtrack_punit = current;
          current->info.or.SR = seq_start;
          current->hit = seq_start;
          current->info.or.alt = 1;
          current = current->info.or.or1;
          try_again = 1;
          break;
      }
      if(try_again){
        try = 1;
      } else {
        try = 0;
      }
    }

    //BACKTRACK segment
    if(backtrack){
      backtrack_again = 0;

      if (!tmp_backtrack_punit)
        return 0;
      else {
        current = tmp_backtrack_punit;
        tmp_backtrack_punit = current->BR;
        seq_start = current->hit;
        switch (current->type) {
          case RANGE_PUNIT:
            if ((current->info.range.nxt <= current->info.range.min + current->info.range.width) &&
                (seq_start + current->info.range.nxt - 1 <= seq_end)) {
              seq_start = seq_start + current->info.range.nxt++;
              tmp_backtrack_punit = current;
              success = 1;
              break;
            }
            else if (((++current->hit + current->info.range.min - 1) <= seq_end) && !current->anchored) {
              current->info.range.nxt = current->info.range.min + 1;
              seq_start = current->hit + current->info.range.min;
              tmp_backtrack_punit = current;
              success = 1;
              break;
            }
            else {
              backtrack_again = 1;
              break;
            }
            break;

          case INV_REP_PUNIT:
            if (current->info.repeat.ins_left && first_backtrack && 0) {
              current->info.repeat.ins_left--;
              current->mlen--;
              seq_start = current->hit + current->mlen;
              success = 1;
              break;
            } else {
              backtrack_again = 1;
              break;
            }
          case REPEAT_PUNIT:
            if (current->info.repeat.ins_left && first_backtrack && 0) {
              current->info.repeat.ins_left--;
              current->mlen--;
              seq_start = current->hit + current->mlen;
              success = 1;
              break;
            } else {
              backtrack_again = 1;
              break;
            }
          case COMPL_PUNIT:
            if (current->info.compl.ins_left && first_backtrack && 0) {
              current->info.compl.ins_left--;
              current->mlen--;
              seq_start = current->hit + current->mlen;
              success = 1;
              break;
            } else {
              backtrack_again = 1;
              break;
            }
          case SIM_PUNIT:
            if (current->info.sim.ins_left && first_backtrack) {
              current->info.sim.ins_left--;
              current->mlen--;
              seq_start = current->hit + current->mlen;
              success = 1;
              break;
            }
          case LLIM_PUNIT:
          case ANY_PUNIT:
          case EXACT_PUNIT:
          case WEIGHT_PUNIT:
            seq_start++;
            try=1;
            break;

          case OR_PUNIT:
            if (current->info.or.alt == 1) {
              current->info.or.alt = 2;
              current = current->info.or.or2;
              try=1;
              break;
            }
            else if ((!current->anchored) && (current->info.or.alt == 2)) {
              seq_start++;
              try=1;
              break;
            }
            else {
              backtrack_again = 1;
              break;
            }
            break;
        }
      }
      if(backtrack_again){
        backtrack = 1;
      } else {
        backtrack = 0;
      }
    }

    //SUCESS segment
    if(success){
      current->mlen = seq_start - current->hit;
      pu1=next_punit(current);
      if (pu1) {
        current = pu1;
        try = 1;
      } else {
        leave = 1;
      }
      success = 0;
    }
  }

  //LEAVE segment
  i = collect_hits (current, revhits);
  for (j = 0; i--; j++)
    hits[j] = revhits[i];
  hits[j] = seq_start;
  *backtrack_punit = tmp_backtrack_punit;
  return j;
}

/**
 * @brief Writes the hits of the matches to revhits in reverse order.
 *
 * @param pu The last pattern unit of the pattern.
 * @param revhits Where the reported hits will be placed.
 *
 * @return How many hits there are.
 */
int
collect_hits (pu, revhits)
     struct punit *pu;
     char **revhits;
{
  struct punit *last;
  char **p = revhits;

  last = NULL;
  while (pu) {
    if (pu->type != OR_PUNIT) {
      *(p++) = pu->hit;
      last = pu;
      pu = pu->prev_punit;
    }
    else {
      if (pu->nxt_punit == last) {
        if (pu->info.or.alt == 1)
          pu = pu->info.or.or1;
        else
          pu = pu->info.or.or2;
        while (pu->nxt_punit)
          pu = pu->nxt_punit;
        last = NULL;
      }
      else {
        last = pu;
        pu = pu->prev_punit;
      }
    }
  }
  return p - revhits;
}

#define MIN(A,B) ((A) < (B) ? (A) : (B))

/**
 * @brief Creates the reverse complement of data and saves the pointer to
 *        the start in result. When reverse complements with specified rule
 *        set, the character is set to the index into the reule set for the 
 *        previously matched char.
 *
 * @param rule_set Whether or not a rule set is specified, is -1 when not.
 * @param data The data to reverse complement
 * @param len Length of the resulting data
 * @param result When finished points to the start of the reverse complemented
 *        data
 */
void
rev_compl_data (rule_set, data, len, result)
     int rule_set;
     unsigned char *data, *result;
     int len;
{
  unsigned int i;

  while (len > 0) {
    i = *(data++);
    if (rule_set == -1)
      result[--len] = (unsigned char) ((i >> 4) & 15);
    else
      result[--len] = (unsigned char) known_char_index[(i & 15)] << 2;
  }
}

#define Stack(N) {stack[nxtent].p1=pat_data; stack[nxtent].p2=seq_data; \
                  stack[nxtent].n1=pat_len;  stack[nxtent].n2=seq_len;   \
                  stack[nxtent].mis=max_mis; stack[nxtent].ins=max_ins; \
                  stack[nxtent].del=max_del; stack[nxtent++].next_choice=N;}

/**
 * @brief Performs a loose match (where a maximum number of mis, del or ins is
 *        specified, e.g. ACGCT[mis,del,ins]) on one_data against the region
 *        pointed at by two_data.
 *
 * loose_match works by first reverse complementing the data if specified.
 * Then it's checking if we are dealing with the special case where max_ins=max_del=0.
 * If we do, we try to match char by char and if, at some point, it does not match, the number
 * of allowed mismatches is decremented, both data pointers a incremented and we continue we
 * run out of allowed mismatches or a match is found.
 *
 * If a maximum number of insertions and deletions is specified, loose_match keeps track of a stack
 * of states, where a new state is added each time a choice is made.
 * When insertions and/or deletions are allowed, we start by checking if we have an exact match
 * of the chars at pat_data and seq_data. If we do, both pointers are incremented and we continue
 * seaching for a new match.
 *
 * If not, we start by using our up the allowed number of mismatches and saving a state depending
 * on whether ins or del are allowed. This state is then saved in a stack, which will be used
 * to backtrack to, if no possible match was found with the choices made.
 *
 * With insertions allowed, these will be used up next. A state will be saved if deletions are
 * available, the number of allowed insertion is decremented and the pointer to the pattern char
 * is incremented.
 *
 * If deletions are allowed, these will be used up last. The number of allowed deletions is
 * decremented and the pointer to the data is incremented.
 *
 * If a point is reached where no possible match could be found, the most recent state is popped
 * from the stack and the choice set by next_choice is then tried instead.
 *
 * @param pat_data Pointer to the pattern to be matched against the data region.
 * @param pat_len Length of the pattern
 * @param seq_data The data region to be matched by pat_data
 * @param seq_len Length of the data region
 * @param max_ins Maximum number of insertions allowed
 * @param max_del Maximum number of deletions allowed
 * @param max_mis Maximum number of mismatches allowed
 * @param rule_set A specified rule set, if any.
 * @param match_range How much more could be matched besides the minimum length
 *                    returned.
 * @param compl_flag Flag that specifies if the data should be reverse complemented,
 *                   or not.
 *
 * @return The length of the least amount of two_data that could be matched.
 */
int
loose_match (pat_data, pat_len, seq_data, seq_len,
       max_ins, max_del, max_mis, rule_set, match_range, compl_flag, ins_left)
     unsigned char *pat_data, *seq_data;
     int pat_len, seq_len;
     int max_ins, max_del, max_mis;
     int *match_range;
     int compl_flag, rule_set;
     int *ins_left;
{
  unsigned char result[MAX_CODES];
  unsigned char *start_seq_data = seq_data;
  int i, nxtent, match;
  struct stackent {
    unsigned char *p1, *p2;
    int n1, n2;
    int mis, ins, del;
    int next_choice;
  } stack[100];
  match = 0;
  if (compl_flag) {
    /* FIX: 4/25/92 RAO; complements of ambiguous regions fail */
    for (i = 0; i < pat_len; i++)
      if (!KnownChar (*(pat_data + i) & 15))
        return 0;
    rev_compl_data (rule_set, pat_data, pat_len, result);
    pat_data = result;
  }
  /* special-case for ins=del=0 */
  if ((max_ins == 0) && (max_del == 0)) {
    if (pat_len > seq_len) {
      return 0;
    }
    for (i = pat_len; i >= 1; i--) {
      if (!KnownChar ((*seq_data) & 15) || (!ExMatches (rule_set, *seq_data, *pat_data) && (--max_mis < 0)))
        return 0;
      else {
        seq_data++;
        pat_data++;
      }
    }
    *ins_left = 0;
    return (seq_data - start_seq_data) + 1;
  }

  nxtent = 0;
  while (seq_len || nxtent) {
    if (seq_len && pat_len && KnownChar ((*seq_data) & 15) &&
                  ExMatches (rule_set, *seq_data, *pat_data)) {
      seq_data++;
      pat_data++;
      seq_len--;
      if (!(--pat_len)){
        match = 1;
        break;
      }
    }
    else if (max_mis && (pat_len >= 1) && (seq_len >= 1)) {
      if (max_ins) {
        Stack (1);
      }
      else if (max_del) {
        Stack (2);
      }
      max_mis--;
      pat_data++;
      seq_data++;
      pat_len--;
      seq_len--;
      if (!pat_len) {
        match = 1;
        break;
      }
    }
    else if ((max_ins) && (pat_len >= 1)) {
      if ((max_del) && (seq_len >= 1)) {
        Stack (2);
      }
      max_ins--;
      pat_data++;
      pat_len--;
      if (!pat_len){
        match = 1;
        break;
      }
    }
    else if ((max_del) && (seq_len >= 1)) {
      max_del--;
      seq_data++;
      seq_len--;
      if (!pat_len){
        match = 1;
        break;
      }
    }
    else if (nxtent--) {
      pat_data = stack[nxtent].p1;
      seq_data = stack[nxtent].p2;
      pat_len = stack[nxtent].n1;
      seq_len = stack[nxtent].n2;
      max_mis = stack[nxtent].mis;
      max_ins = stack[nxtent].ins;
      max_del = stack[nxtent].del;
      if (stack[nxtent].next_choice == 1)
      {
        if (max_del)
          stack[nxtent++].next_choice = 2;
        max_ins--;
        pat_data++;
        pat_len--;
        if (!pat_len){
          match = 1;
          break;
        }
      }
      else {
        max_del--;
        seq_data++;
        seq_len--;
        if (!pat_len){
          match = 1;
          break;
        }
      }
    }
    else
      return 0;
  }
  if (max_ins >= pat_len && !match) {
    max_ins -= pat_len;
    match = 1;
  }
  if (match) {
    *ins_left = max_ins;
    return (seq_data - start_seq_data) + 1;
  }
  return 0;
}

#ifdef REVERSE_PUNIT

#define FIXBUFF *buff += strlen(*buff)

/* WARNING: VERY UNRELIABLE */
/**
 * @brief 
 *
 * @param buff
 */
void
reverse_punit (buff)
     char *buff;
{
  punit_t *pu;

  strcpy (buff, "Reversed: ");
  buff += strlen (buff);
  for (pu = ad_pu_s; pu; pu = pu->nxt_punit)
    reverse_one_punit (pu, &buff);
}

#define NEXTBUFF(X) **buff=((X)), *buff += 1, **buff='\0'

/**
 * @brief 
 *
 * @param pu
 * @param buff
 */
void
reverse_one_punit (pu, buff)
     punit_t *pu;
     char **buff;
{
  punit_t *pu1;
  int i;

  for (i = 0; i < MAX_NAMES; i++)
    if (pu == names[i]) {
      sprintf (*buff, " p%d=", i);
      FIXBUFF;
      break;
    }
  if (i == MAX_NAMES)
    NEXTBUFF (' ');
  switch (pu->type) {
  case RANGE_PUNIT:
    (void) sprintf (*buff, "%d...%d",
        pu->info.range.min,
        pu->info.range.min + pu->info.range.width);
    break;

  case EXACT_PUNIT:
    for (i = 0; i < pu->info.exact.len; i++)
      NEXTBUFF (display_seq_char[0xF & pu->info.exact.code[i]]);
    break;

  case LLIM_PUNIT:
    exit (2);
    break;

  case COMPL_PUNIT:
    (void) sprintf (*buff, "~p%d[%d,%d,%d]",
        pu->info.compl.of,
        pu->info.compl.mis,
        pu->info.compl.ins, pu->info.compl.del);
    break;

  case REPEAT_PUNIT:
    (void) sprintf (*buff, "p%d[%d,%d,%d]",
        pu->info.repeat.of,
        pu->info.repeat.mis,
        pu->info.repeat.ins, pu->info.repeat.del);
    break;

  case SIM_PUNIT:
    for (i = 0; i < pu->info.sim.len; i++)
      NEXTBUFF (display_seq_char[0xF & pu->info.sim.code[i]]);

    (void) sprintf (*buff, "[%d,%d,%d]",
        pu->info.sim.mis, pu->info.sim.ins, pu->info.sim.del);
    break;

  case WEIGHT_PUNIT:
    if (pu->info.wvec.tupsz != 4)
      exit (3);
    NEXTBUFF ('{');

    for (i = 0; i < pu->info.wvec.len; i++) {
      wv_t *wts = (wv_t *) pu->info.wvec.vec + i;
      if (i)
  NEXTBUFF (',');
      sprintf (*buff, "(%d,%d,%d,%d)", wts->aw, wts->cw, wts->gw, wts->tw);
      FIXBUFF;
    }

    NEXTBUFF ('}');
    (void) sprintf (*buff, ">%d", pu->info.wvec.cutoff);
    break;


  case OR_PUNIT:
    NEXTBUFF ('(');
    for (pu1 = pu->info.or.or1; pu1; pu1 = pu1->nxt_punit)
      reverse_one_punit (pu1, buff);
    NEXTBUFF ('|');
    for (pu1 = pu->info.or.or2; pu1; pu1 = pu1->nxt_punit)
      reverse_one_punit (pu1, buff);
    NEXTBUFF (')');
    break;
  }
  FIXBUFF;
}
#endif


/**
 * @brief Will set the middle_punit to the best pivot pattern unit and will
 *        remove its previous and next punit.
 *
 * @param first_punit First pattern unit of the pattern.
 * @param max_before_len Will be set to the maximum matching length of the
 *                       punits before the pivot pattern unit.
 * @param max_after_len Will be set to the maximum matching length of Then
 *                      punits after the pivot punit.
 *
 * @return 0 if first punit is best, or if there are no pivot is found.
 *         otherwise returns 1.
 */
int
punit_list_split (first_punit, max_before_len, max_after_len)
    punit_t  *first_punit;
    int *max_before_len, *max_after_len;
{
  punit_t *iterator;
  int i;
  if (!find_pivot_punit(max_before_len, max_after_len, &best_index)){
    return 0;
  }
  iterator = first_punit;
  for(i = 1; i < best_index; i++){
    iterator = iterator->nxt_punit;
  }
  middle_punit = *(iterator->nxt_punit);
  middle_punit.nxt_punit = NULL;
  middle_punit.prev_punit = NULL;
  return 1;
}


/**
 * @brief Will perform a pattern optmised match. middle_punit needs to be set
 *        with punit_list_split before this is run.
 *
 * @param first_punit The first pattern unit of the pattern.
 * @param start_c Pointer to the first character of the sequence.
 * @param end_c Pointer to the last character of the sequence.
 * @param hits Will hold the hits made.
 * @param first Whether it is the first run or not.
 * @param BR1 The backtracking punit, not actually used.
 *
 * @return How many bases matched, if no match was found returns 0.
 */
int
pattern_optimized_match (first_punit, start_c, end_c, hits, first, BR1)
     punit_t **BR1;
     punit_t *first_punit;
     char **start_c, *end_c;
     char *hits[];
     int first;
{
  int i,j,k;
  *BR1 = NULL;
  if (!first){
    *start_c = middle_punit.hit + middle_punit.mlen;
  }
  middle_punit.anchored = 0;
  char *start_save = *start_c;
  char *start_tmp;
  char *mid_hit;
  while (*start_c <= end_c) {
    if (pattern_match (&middle_punit, *start_c, end_c, &mid_hit, 1, BR1)) {
      start_tmp = *start_c;
      *BR1 = NULL;
      if (mid_hit - max_before_mlen > start_save) {
        *start_c = mid_hit - max_before_mlen;
      } else if (mid_hit - max_before_mlen < past_last) {
        *start_c = past_last;
      } else {
        *start_c = start_save;
      }
      if (i = pattern_match(ad_pu_s, *start_c, mid_hit + max_after_mlen + 1, hits, 1, BR1)) {
        past_last = hits[i + 1];
        return i;
      } else {
        *start_c = mid_hit + 1;
      }
    } else {
      return 0;
    }
  }
  return 0;
}
