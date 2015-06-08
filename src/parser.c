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

/* the following variable hides the relationship between parsing punits
   and executing searches.  These routines compose the interface:

   int parse_dna_cmd(line)      returns max # characters that could be matched
                                    if successful parse (sets ad_pu_s)
   char *line;                  returns 0 if unsuccessful parse

   int parse_peptide_cmd(line)  returns max # characters that could be matched
                                    if successful parse (sets ad_pu_s)
   char *line;                  returns 0 if unsuccessful parse
*/

int initialized = 0;

int punit_sequence_type;

/* state is used to retain state for a search -- needed for parallel
   Prolog versions
*/
struct state *
get_state () {
  return ((struct state *) malloc (sizeof (struct state)));
}

void
free_state (st)
     struct state *st;
{
  free (st);
}

struct punit *ad_pu_s = NULL;  /* CONNECTS PARSER TO SCAN ROUTINES */

/* ============================================================= */

/* prototypes  */
#ifdef CC
#define PROTO(X) ()
#else
#define PROTO(X) X
#endif

//void free_state PROTO ((struct state * st));
void reverse_punit PROTO ((char *buff));
void reverse_one_punit PROTO ((punit_t * pu, char **buff));
char *num PROTO ((int *n, char *p));
int max_mat PROTO ((struct punit * pu));
void set_anchors PROTO ((struct punit * pu));
void set_anchors_on PROTO ((struct punit * pu));
punit_t *parser PROTO ((punit_t ** pu_s, char **cv, int **iv, char *line));
char *punit_cons_list PROTO ((punit_t ** pu_s, char **cv, int **iv, char *p));
char *punit_cons PROTO ((punit_t ** pu_s, char **cv, int **iv, char *p));
char *punit_list PROTO ((punit_t ** pu_s, char **cv, int **iv, char *p));
char *punit_parse PROTO ((punit_t ** pu_s, char **cv, int **iv, char *p));
char *or_pat PROTO ((punit_t ** pu_s, char **cv, int **iv, char *p));
char *char_pat PROTO ((char **enc_pat, char *p));
char *range_pat PROTO ((int *min, int *max, char *p));
char *elipses PROTO ((char *p));
char *wt_pat
PROTO ((int *maxwt, int *cutoff, int *tupsz, int **t, int *n, char *p));
char *wt_template PROTO ((int *tupsz, int **t, int *n, char *p));
char *n_tuple PROTO ((int *n, int *t, char *p));
char *repeat_pat PROTO ((int *n, int *mis, int *ins, int *del, char *p));
char *inv_rep_pat PROTO ((int *n, int *mis, int *ins, int *del, char *p));
char *compl_pat
PROTO ((int *rs, int *n, int *mis, int *ins, int *del, char *p));
char *rule_id PROTO ((int *rs, char *p));
char *compl_id PROTO ((int *rs, int *n, char *p));
char *sim_pat PROTO ((char **enc_pat, int *mis, int *ins, int *del, char *p));
char *misinsdel PROTO ((int *mis, int *ins, int *del, char *p));
char *dna_pat PROTO ((char **enc_pat, char *p));
char *name_assgn PROTO ((int *n, char *p));
char *name_id PROTO ((int *n, char *p));
char *white_space PROTO ((char *p));
int parse_cmd PROTO ((char *line));
char *wt_pat
PROTO ((int *maxwt, int *cutoff, int *tupsz, int **t, int *n, char *p));
char *punit_cons_list PROTO ((punit_t ** pu_s, char **cv, int **iv, char *p));
int max_mats PROTO ((struct punit * pu));

/* ============================================================= */

/**
* @brief Given a character pointer, this will return the pointer to the next non-whitespace character.
*
* @param p Pointer to a character.
*
* @return Pointer to the first non-whitespace character after the given pointer.
*/
char *
white_space (p)
     char *p;
{
  while (isspace (*p))
    p++;
  return (p);
}

/**
* @brief Converts a string number to an integer.
*
* @param n This will, after execution, point to the converted integer.
* @param p Pointer to the first character of the string integer.
*
* @return a pointer thats points after the string, returns NULL if invalid number.
*/
char *
num (n, p)
     int *n;
     char *p;
{
  int i;
  int sign;

  if (*p == '-') {
    sign = -1;
    p++;
  }
  else
    sign = 1;

  if ((*p < '0') || (*p > '9'))
    return NULL;
  else {
    i = *(p++) - '0';
    while ((*p >= '0') && (*p <= '9'))
      i = (i * 10) + (*(p++) - '0');
    *n = i * sign;
    return p;

  }
}

/**
* @brief Checks if p is a valid name id e.g. p1, p2, ... , pn, where 0 <= n <= MAX_NAMES.
*
* @param n Will hold the pattern number, e.g. if p points to "p5", n will point to the integer 5.
* @param p The pointer to the first character of a pattern name.
*
* @return The pointer to the character after the pattern name. NULL if unsuccessful.
*/
char *
name_id (n, p)
     int *n;
     char *p;
{
  if ((*(p++) == 'p') && (p = num (n, p)) && (*n >= 0) && (*n <= MAX_NAMES))
    return p;
  else
    return NULL;
}

/**
* @brief Assign a name id if formatted correctly e.g. p1=, whitespace is allowed before and after =.
*
* @param n Will hold the pattern number, e.g. if p points to "p5", n will point to the integer 5.
* @param p The pointer to the first character of a pattern name.
*
* @return The pointer to the character after the =. NULL if unsuccessful.
*/
char *
name_assgn (n, p)
     int *n;
     char *p;
{
  if ((p = name_id (n, p)) && (p = white_space (p)) && (*(p++) == '='))
    return p;
  else
    return NULL;
}

/**
* @brief Encrypts the characters in a literal pattern, to their respective bit values,
*        e.g. ACG is stored in enc_pat as [A_BIT, C_BIT, G_BIT]. The value of these
*        bits are described in scanner.h
*
* @param enc_pat Is the pointer to the list of all encrypted patterns, the new pattern will be placed at this pointer.
* @param p Is the pointer to the first character of a literal pattern.
*
* @return the pointer to the immediate character after the literal pattern. Returns NULL if the pattern cannot be parsed.
*/
char *
dna_pat (enc_pat, p)
     char **enc_pat;
     char *p;
{
  char *p1;

  if (punit_sequence_type != DNA)
    return (NULL);
  for (p1 = *enc_pat;
       (*p != '\0') && (*p != ' ') && (*p != '\t') && (*p != '[')
       && (*p != '\n') && (*p != ')'); p++) {
    *p1 = punit_to_code[*p];
    if (*p1)
      p1++;
    else
      return NULL;
  }
  if (p1 > *enc_pat) {
    *p1 = '\0';
    *enc_pat = p1 + 1;
    return p;
  }
  else {
    return NULL;
  }
}

/**
* @brief Checks if given miss, insertions and deletions are formatted correctly.
*        E.g. [0,2,4].
*
* @param mis Will upon successful execution hold the given number of misses.
* @param ins Will upon successful execution hold the given number of insertions.
* @param del Will upon successful execution hold the given number of deletions.
* @param p The pointer to the first character of the specification.
*
* @return The pointer to the immediate character after the specification,
*         returns NULL if the specification is wrongly formatted.
*/
char *
misinsdel (mis, ins, del, p)
     int *mis, *ins, *del;
     char *p;
{
  if ((*(p++) == '[') && (p = num (mis, p)) && (*(p++) == ',') &&
      (p = num (ins, p)) && (*(p++) == ',') &&
      (p = num (del, p)) && (*(p++) == ']'))
    return p;
  else
    return NULL;
}

/**
* @brief Parses a DNA or PEPTIDE pattern of type AGC[1,2,3] or AGCUCAG.
*
* @param enc_pat Where the encrypted character bits of the pattern will be held.
* @param mis Will hold the given number of misses, 0 if none are given.
* @param ins Will hold the given number of insertions, 0 if none are given.
* @param del Will hold the given number of deletions, 0 if none are given.
* @param p The pointer to the first character of the pattern.
*
* @return The pointer to the immediate character after the pattern.
*         Returns NULL if parsing fails.
*/
char *
sim_pat (enc_pat, mis, ins, del, p)
     char **enc_pat, *p;
     int *mis, *ins, *del;
{
  char *p1;

  if (((punit_sequence_type == DNA) && (p1 = dna_pat (enc_pat, p))) ||
      ((punit_sequence_type == PEPTIDE) && (p1 = char_pat (enc_pat, p)))) {
    if (p = misinsdel (mis, ins, del, p1))
      return p;
    else {
      *mis = *ins = *del = 0;
      return p1;
    }
  }
  else
    return NULL;
}

/**
* @brief Parses complement pattern name ids, e.g. ~p1 or r1~p1.
*
* @param rs Will hold the rule set id, if none is given, will be -1.
* @param n Will hold the pattern id.
* @param p The pointer to the first character of the pattern.
*
* @return The pointer to the immediate character after the pattern.
*         Returns NULL if parsing fails.
*/
char *
compl_id (rs, n, p)
     int *rs;
     int *n;
     char *p;
{
  if (punit_sequence_type != DNA)
    return NULL;

  if (*p == '~') {
    *rs = -1;
    return name_id (n, p + 1);
  }
  else if ((p = rule_id (rs, p)) && (*p == '~'))
    return name_id (n, p + 1);
  else
    return NULL;
}

/**
* @brief Parses reverse complement patterns of the following formats:
*        ~p1=[1,2,3], r1~p1[1,2,3], ~p1 and r1~p1.
*
* @param rs Will hold the id of the parsed ruleset, -1 if none is given.
* @param n Will hold the id of the pattern name.
* @param mis Will hold the given misses, 0 if none are specified.
* @param ins Will hold the given insertions, 0 if none are specified.
* @param del Will hold the given deletions, 0 if none are specified.
* @param p The pointer to the first character of the complement pattern.
*
* @return The pointer to the immediate character after the pattern.
*         Returns NULL if parsing fails.
*/
char *
compl_pat (rs, n, mis, ins, del, p)
     int *rs, *n, *mis, *ins, *del;
     char *p;
{
  char *p1;

  if (punit_sequence_type != DNA)
    return (NULL);

  if (p = compl_id (rs, n, p)) {
    if (p1 = misinsdel (mis, ins, del, p))
      return p1;
    else {
      *mis = *ins = *del = 0;
      return p;
    }
  }
  else
    return NULL;
}

/**
* @brief Parses patterns of the format p1[0,1,2] or p1.
*
* @param n Will hold the id of the pattern name.
* @param mis Will hold the given misses, 0 if none are specified.
* @param ins Will hold the given insertions, 0 if none are specified.
* @param del Will hold the given deletions, 0 if none are specified.
* @param p The pointer to the first character of the pattern.
*
* @return The pointer to the immediate character after the pattern.
*         Returns NULL if parsing fails.
*/
char *
repeat_pat (n, mis, ins, del, p)
     int *n, *mis, *ins, *del;
     char *p;
{
  char *p1;

  if (p = name_id (n, p))
    if (p1 = misinsdel (mis, ins, del, p))
      return p1;
    else {
      *mis = *ins = *del = 0;
      return p;
    }
  else
    return NULL;
}

/**
* @brief Parses patterns of the format <p1[0,1,2] or <p1.
*
* @param n Will hold the id of the pattern name.
* @param mis Will hold the given misses, 0 if none are specified.
* @param ins Will hold the given insertions, 0 if none are specified.
* @param del Will hold the given deletions, 0 if none are specified.
* @param p The pointer to the first character of the pattern.
*
* @return The pointer to the immediate character after the pattern.
*         Returns NULL if parsing fails.
*/
char *
inv_rep_pat (n, mis, ins, del, p)
     int *n, *mis, *ins, *del;
     char *p;
{
  char *p1;

  if (*(p++) != '<')
    return NULL;

  if (p = name_id (n, p))
    if ((punit_sequence_type == DNA) && (p1 = misinsdel (mis, ins, del, p)))
      return p1;
    else {
      *mis = *ins = *del = 0;
      return p;
    }
  else
    return NULL;
}

/**
* @brief Parses tuples of the format (0,6,8,9).
*
* @param n Will hold the length of the tuple.
* @param t Will point to the list of integers in the tuple.
* @param p The pointer to the first character of the tuple.
*
* @return The pointer to the immediate character after the tuple.
*         Returns NULL if parsing fails.
*/
char *
n_tuple (n, t, p)
     int *n, *t;
     char *p;
{
  if ((*(p++) == '(') && (p = white_space (p)) && (p = num (t, p))) {
    *n = 1;
    while (TRUE) {
      if (*p == ',') {
        p++;
        p = white_space (p);
        if (p = num (t + (*n), p)) {
          (*n)++;
          p = white_space (p);
        }
        else
          return NULL;
      }
      else if (*(p++) == ')')
        return p;
      else
        return NULL;
    }
  }
  else
    return NULL;
}

/**
* @brief Parses a weight matrix of the following format:
*        {(16,0,84,0),(57,10,29,4),(0,80,0,20),(95,0,0,5),
*         (0,100,0,0),(18,60,20,2),(0,0,100,0),(0,50,50,0)}
*
* @param tupsz Will hold the size of the individual tuples in the matrix,
*              these being 4, 20 or 21.
* @param t Will hold the weight matrix values, where the first value of
*          each column will be at t[n*tupsz].
* @param n Will hold the amount of tuples in the weight matrix.
* @param p The pointer to the first character of the weight matrix.
*
* @return The pointer to the immediate character after the weight matrix.
*         Returns NULL if parsing fails.
*/
char *
wt_template (tupsz, t, n, p)
     int *tupsz, **t, *n;
     char *p;
{
  int sz;

  *n = 0;
  if ((*(p++) == '{') && (p = white_space (p))
  && (p = n_tuple (tupsz, *t, p))) {
    *t += *tupsz;
    (*n)++;
    while ((p = white_space (p)) && (*p == ',')) {
      p++;
      if ((p = white_space (p)) &&
         (p = n_tuple (&sz, *t, p)) && (sz == *tupsz)) {
        (*n)++;
        *t += sz;
      }
      else
        return NULL;
    }
    if ((*(p++) == '}')
    && ((*tupsz == 4) || (*tupsz == 20) || (*tupsz == 21)))
      return p;
    else
      return NULL;
  }
  else
    return NULL;
}

/**
* @brief Parses weight matrices of the forms:
*            x > {...} > y
*        OR
*                {...} > y
*
* @param maxwt Will hold the maximum total weights for a match.
* @param cutoff Will Hold the minimum total weights for a match.
* @param tupsz Will hold the size of the individual tuples in the matrix,
*              these being 4, 20 or 21.
* @param t Will hold the weight matrix values, where the first value of
*          each column will be at t[n*tupsz].
* @param n Will hold the amount of tuples in the weight matrix.
* @param p The pointer to the first character of the weight matrix.
*
* @return The pointer to the immediate character after the weight matrix.
*         Returns NULL if parsing fails.
*/
char *
wt_pat (maxwt, cutoff, tupsz, t, n, p)
     int *maxwt, *cutoff, *tupsz, **t, *n;
     char *p;
{
  char *p1;

  if ((p1 = num (maxwt, p)) && (p1 = white_space (p1)) && (*(p1++) == '>') &&
      (p1 = white_space (p1)) && (p1 = wt_template (tupsz, t, n, p1))
      && (p1 = white_space (p1)) && (*(p1++) == '>')
      && (p1 = white_space (p1)) && (p1 = num (cutoff, p1))) {
    return p1;
  }

  *maxwt = 1000000;

  if ((p = wt_template (tupsz, t, n, p)) && (p = white_space (p)) &&
      (*(p++) == '>') && (p = white_space (p)) && (p = num (cutoff, p)))
    return p;
  else
    return NULL;
}

/**
* @brief Parses the exact pattern: "..."
*
* @param p The pointer to the first ".".
*
* @return The pointer to the immediate character after the three dots.
*         Returns NULL if it is not successfully matched.
*/
char *
elipses (p)
     char *p;
{
  if ((*(p++) == '.') && (*(p++) == '.') && (*(p++) == '.'))
    return p;
  else
    return NULL;
}


/**
* @brief Matches the word in str to the pointer p.
*
* @param str Is the word that needs to be matched at p.
* @param p Is the pointer to the first letter of a given word in the pattern.
*
* @return The pointer to the character immediately after the word.
*         Returns NULL if the word is not successfully matched.
*/
char *
word (str, p)
     char *str;
     char *p;
{
  while (*str) {
    if (*(str++) != *(p++))
      return NULL;
  }
  return p;
}

/**
* @brief Matches patterns of the form length(p1+p2+p6) < 14.
*
* @param v Will hold the pattern ids.
* @param max Will hold the maximum length of the pattern.
* @param p The pointer to the first character of the pattern.
*
* @return The pointer to the character immediately after the length pattern.
*         Returns NULL if the pattern is cannot be parsed.
*/
char *
llim_pat (v, max, p)
     int *v, *max;
     char *p;
{
  char *p1;

  *v = 1;
  if ((p = word ("length(", p)) && (p = name_id (v + 1, p))) {
    while ((p1 = white_space (p)) && (*(p1++) == '+')
     && (p1 = white_space (p1))
     && (p1 = name_id ((v + 1) + ((*v)++), p1)))
      p = p1;

    if ((p = white_space (p)) && (*(p++) == ')') && (p = white_space (p))
  && (*(p++) == '<') && (p = white_space (p)) && (p = num (max, p)))
      return p;
    else
      return NULL;
  }
  return NULL;
}

/**
* @brief Parses a range pattern of the form 2...5 or 2 ... 5.
*
* @param min Will hold the minimum of the range.
* @param max Will hold the maximum of the range.
* @param p The pointer to the character holding the first number of the range.
*
* @return The pointer to the character immediately after the range.
*         Returns NULL if the range pattern is of from format.
*/
char *
range_pat (min, max, p)
     int *min, *max;
     char *p;
{
  if ((p = num (min, p)) && (p = white_space (p)) && (p = elipses (p)) &&
      (p = white_space (p)) && (p = num (max, p)))
    return p;
  else
    return NULL;
}

/**
* @brief Matches exact patterns when the sequence type is a peptide.
*
* @param enc_pat Will hold a copy of the pattern given.
* @param p The pointer to the first character of the pattern.
*
* @return The pointer to the character immediately after the exact pattern.
          Returns NULL if the patterns contains wrong characters.
*/
char *
char_pat (enc_pat, p)
     char **enc_pat, *p;
{
  char *p1;

  if (punit_sequence_type == PEPTIDE && isalpha (*p)) {
    for (p1 = *enc_pat; isalpha (*p); p++) {
      *(p1++) = (char) toupper ((int) *p);
    }

    *p1 = '\0';
    *enc_pat = p1 + 1;
    return p;
  }
  else
    return NULL;
}

/**
* @brief Matches patterns of the form any(HQD) and notany(HK).
*
* @param acm Will hold an alphabetical character mask, where the bits
             corresponding to their position in the alphabet is 1 if the
             character is allowed and 0 otherwise.
* @param p The pointer to the first character of the pattern.
*
* @return The pointer to the character after the pattern.
*         Returns NULL if the pattern is not matched.
*/
char *
any_pat (acm, p)
     char *p;
     long *acm;
{
  char *p1;
  int i;
  long cm;

  if (punit_sequence_type == PEPTIDE && strncmp ("any(", p, 4) == 0) {
    p += 4;
    for (cm = 0; isalpha (*p); p++) {
      cm |= 1 << (toupper ((int) *p) - 'A');
    }
    if (*(p++) == ')') {
      *acm = cm;
      return p;
    }
    else
      return NULL;
  }
  else if (punit_sequence_type == PEPTIDE && strncmp ("notany(", p, 7) == 0) {
    p += 7;
    for (cm = 0x3ffffff; isalpha (*p); p++) {
      cm &= ~((long) (1 << (toupper ((int) *p) - 'A')));
    }
    if (*(p++) == ')') {
      *acm = cm;
      return p;
    }
    else
      return NULL;
  }
  else
    return NULL;
}

/**
* @brief Parses OR patterns of the following format ( pattern | pattern ).
*
* @param pu_s Will be the pattern unit of the OR pattern.
* @param cv Will hold the exact or fuzzy matches patterns.
* @param iv Will hold the weight matrixes generated in sub patterns.
* @param p The pointer to the first character of the OR pattern.
*
* @return The pointer to the character immediately after the OR pattern.
*         Returns NULL if the OR pattern is not successfully parsed.
*/
char *
or_pat (pu_s, cv, iv, p)
     punit_t **pu_s;
     char **cv, *p;
     int **iv;
{
  punit_t *or_pu, *p1, *p2;
  char *punit_cons_list ();

  or_pu = (*pu_s)++;
  if ((*(p++) == '(') && (p = white_space (p)) &&
      (p1 = *pu_s) && (p = punit_cons_list (pu_s, cv, iv, p)) &&
      (p = white_space (p)) && (*(p++) == '|') && (p = white_space (p)) &&
      (p2 = *pu_s) && (p = punit_cons_list (pu_s, cv, iv, p)) &&
      (p = white_space (p)) && (*(p++) == ')')) {
    or_pu->type = OR_PUNIT;
    or_pu->info.or.or1 = p1;
    or_pu->info.or.or2 = p2;
    p1->prev_punit = p2->prev_punit = or_pu;
    or_pu->nxt_punit = or_pu->prev_punit = NULL;
    return (p);
  }
  else {
    *pu_s = or_pu;
     return NULL;
  }
}

/**
* @brief Will convert a DNA character to an integer.
*
* @param p The pointer to the character.
* @param c Will be 0 if a or A, 1 if c or C, 2 if g or G,
           and 3 if t, T, u or U.
*
* @return The pointer to the character after the DNA character.
          Returns NULL if the character cannot be parsed.
*/
char *
dna_char (p, c)
     char *p, *c;
{
  switch (*p) {
  case 'a':
  case 'A':
    *c = 0;
    return p + 1;

  case 'c':
  case 'C':
    *c = 1;
    return p + 1;

  case 'g':
  case 'G':
    *c = 2;
    return p + 1;

  case 't':
  case 'T':
  case 'u':
  case 'U':
    *c = 3;
    return p + 1;

  default:
    return NULL;
  }
}

/**
* @brief Creates a bitwise relationship between two consecutive characters in
*        the DNA string in the rs list.
*
* This works by setting the rs[dna_char(char1) * 4 + dna_char(char2)] = 1.
* For example will TG set rs[3*4 + 2] = rs[14] = 1, there will then be exactly
* 16 different positions, corresponding to all possible sets of DNA characters.
*
* @param p The pointer to two consecutive DNA characters in a ruleset.
* @param rs The list that will hold the ruleset.
*
* @return The pointer to the character after the two consecutive DNA
*         characters. Returns NULL if they cannot be successfully parsed.
*/
char *
bond (p, rs)
     char *p;
     char rs[16];
{
  char c1, c2;

  if ((p = dna_char (p, &c1)) && (p = dna_char (p, &c2))) {
    rs[(c1 << 2) + c2] = 1;
    return p;
  }
  else
    return NULL;
}

/**
* @brief Parses bonds of the form AA,AC,AG,AT,CA,CA,...,TT.
*
* @param p The pointer to the first character in the bonds list.
* @param rs The ruleset to be populated with rules.
*
* @return The pointer to the character after the bonds.
*         Returns NULL if it is not successfully parsed.
*/
char *
parse_bonds (p, rs)
     char *p;
     char rs[16];
{
  if (p = bond (p, rs)) {
    while ((*p == ',') && (p = bond (p + 1, rs)));
    return p;
  }
  else
    return NULL;
}

/**
* @brief Parses rulesets of the form: {AA,AC,AG,AT,CA,CC,CG,CT,GA,...,TT}
*
* This works by setting the rs[dna_char(char1) * 4 + dna_char(char2)] = 1.
* For example will TG set rs[3*4 + 2] = rs[14] = 1, there will then be exactly
* 16 different positions, corresponding to all possible sets of DNA characters.
* These positions will correspond to if the bond is allowed, in the previus example
* T would be allowed to bond with G.
*
* @param p Pointer to the first character of the ruleset.
* @param n The number of the rule set, if r1, then this would be 1.
*
* @return The pointer to the character after the last }.
*         Returns NULL if it is not successfully parsed.
*/
char *
bond_set (p, n)
     char *p;
     int n;
{
  char rs[16];
  int i, j;

  for (i = 0; i < 16; i++)
    rs[i] = 0;

  if ((*(p++) == '{') && (p = parse_bonds (p, rs)) && (*(p++) == '}')) {
    if (!rule_sets[n]) {
      if ((rule_sets[n] = malloc (16)) == NULL) {
        fprintf (stderr, "memory allocation failure\n");
        exit (1);
      }
    }
    for (i = 0; i < 16; i++)
      *(rule_sets[n] + i) = rs[i];
    return p;
  }
  else
    return NULL;
}

/**
* @brief Parses rule ids of the form r1.
*
* @param n When parsed this will contain the rule set id.
* @param p Pointer to the first character of the ruleset.
*
* @return Pointer to the first character after the rule id.
*         Returns NULL if it is not successfully parsed.
*/
char *
rule_id (n, p)
     int *n;
     char *p;
{
  if ((*(p++) == 'r') && (p = num (n, p)) && (*n >= 0) && (*n <= MAX_NAMES))
    return p;
  else
    return NULL;
}

/**
* @brief Parses rule sets of the form r1={AG,CT,GG}
*
* This works by setting the rs[dna_char(char1) * 4 + dna_char(char2)] = 1.
* For example will TG set rs[3*4 + 2] = rs[14] = 1, there will then be exactly
* 16 different positions, corresponding to all possible sets of DNA characters.
* These positions will correspond to if the bond is allowed, in the previus example
* T would be allowed to bond with G.
*
* @param p Pointer to the first character of the ruleset id.
*
* @return Pointer to the character after the last character of the ruleset.
*         Returns NULL if it is not successfully parsed.
*/
char *
parse_rule_set (p)
     char *p;
{
  int n;

  if (punit_sequence_type != DNA)
    return (NULL);

  if ((p = rule_id (&n, p)) && (p = white_space (p)) && (*(p++) == '=')
      && (p = bond_set (p, n)))
    return p;
  else
    return NULL;
}

/**
* @brief Parses a pattern unit at position p.
*
* @param pu_s Will hold parsed pattern units.
* @param cv Will hold characters of exact and fuzzy patterns.
* @param iv Will hold a potential weight table.
* @param p The pointer to the first character of the pattern unit.
*
* @return The pointer to the immediate character after the pattern unit.
*         Returns NULL if it cannot be parsed.
*/
char *
punit_parse (pu_s, cv, iv, p)
     punit_t **pu_s;
     char **cv, *p;
     int **iv;
{
  punit_t *pu1;
  int i, j, k, n, rs;
  char *p1, *p2;
  int *i1;
  int tupsz;

  while (p1 = parse_rule_set (p))
    p = white_space (p1);

  pu1 = (*pu_s)++;
  if (p1 = name_assgn (&i, p)) {
    if (names[i])
      return (NULL);    /* ensure only one occurrence of pN= */
    names[i] = (struct punit *) pu1;
    p = p1;
  }
  p = white_space (p);
  if (*p == '^') {
    pu1->type = MATCH_START;
    p1 = p + 1;
  }
  else if (*p == '$') {
    pu1->type = MATCH_END;
    p1 = p + 1;
  }
  else if (p1 = range_pat (&i, &j, p)) {
    pu1->type = RANGE_PUNIT;
    pu1->info.range.min = i;
    pu1->info.range.width = j - i;
  }
  else if (p1 = any_pat (&(pu1->info.any.code_matrix), p)) {
    pu1->type = ANY_PUNIT;
  }
  else if (p1 =
     llim_pat (pu1->info.llim.llim_vec, &(pu1->info.llim.bound), p)) {
    pu1->type = LLIM_PUNIT;
  }
  else if (p1 = inv_rep_pat (&n, &i, &j, &k, p)) {
    pu1->type = INV_REP_PUNIT;
    pu1->info.repeat.of = n;
    pu1->info.repeat.mis = i;
    pu1->info.repeat.ins = j;
    pu1->info.repeat.del = k;
  }
  else if (p1 = repeat_pat (&n, &i, &j, &k, p)) {
    pu1->type = REPEAT_PUNIT;
    pu1->info.repeat.of = n;
    pu1->info.repeat.mis = i;
    pu1->info.repeat.ins = j;
    pu1->info.repeat.del = k;
  }
  else if ((p2 = *cv) && (p1 = sim_pat (cv, &i, &j, &k, p))) {
    if ((i == 0) && (j == 0) && (k == 0)) {
      pu1->type = EXACT_PUNIT;
      pu1->info.exact.code = p2;
      pu1->info.exact.len = strlen (p2);
    }
    else {
      pu1->type = SIM_PUNIT;
      pu1->info.sim.code = p2;
      pu1->info.sim.len = strlen (p2);
      pu1->info.sim.mis = i;
      pu1->info.sim.ins = j;
      pu1->info.sim.del = k;
    }
  }
  else if (p1 = compl_pat (&rs, &n, &i, &j, &k, p)) {
    pu1->type = COMPL_PUNIT;
    pu1->info.compl.rule_set = rs;
    pu1->info.compl.of = n;
    pu1->info.compl.mis = i;
    pu1->info.compl.ins = j;
    pu1->info.compl.del = k;
  }
  else if ((i1 = *iv) && (p1 = wt_pat (&k, &i, &tupsz, iv, &j, p))) {
    pu1->type = WEIGHT_PUNIT;
    pu1->info.wvec.vec = i1;
    pu1->info.wvec.len = j;
    pu1->info.wvec.cutoff = i;
    pu1->info.wvec.maxwt = k;
    pu1->info.wvec.tupsz = tupsz;
  }
  else {
    *pu_s = pu1;
    return NULL;
  }
  return (p1);
}

/**
* @brief Parses a list of pattern units seperated by whitespace.
*
* @param pu_s Will hold parsed pattern units.
* @param cv Will hold characters of exact and fuzzy patterns.
* @param iv Will hold a potential weight table.
* @param p The pointer to the first character of the first pattern unit.
*
* @return The pointer to the immediate character after the pattern unit.
*         Returns NULL if it cannot be parsed.
*/
char *
punit_list (pu_s, cv, iv, p)
     punit_t **pu_s;
     char **cv;
     int **iv;
     char *p;
{
  punit_t *last, *next;
  char *p1;

  last = NULL;

  p = white_space (p);
  next = *pu_s;
  while (p1 = punit_parse (pu_s, cv, iv, p)) {
    next->prev_punit = last;
    if (last)
      last->nxt_punit = next;
    last = next;
    next = *pu_s;
    p = white_space (p1);
  }
  if (!last)
    return NULL;
  else {
    last->nxt_punit = NULL;
    return p;
  }
}

/**
* @brief Parses OR patterns or parse lists, depending on what pattern is given.
*
* @param pu_s Will hold parsed pattern units.
* @param cv Will hold characters of exact and fuzzy patterns.
* @param iv Will hold a potential weight table.
* @param p The pointer to the first character of the first pattern unit.
*
* @return The pointer to the character immediately after the parsed pattern.
*/
char *
punit_cons (pu_s, cv, iv, p)
     punit_t **pu_s;
     char **cv, *p;
     int **iv;
{
  char *punit_list ();
  char *or_pat ();
  char *p1;

  while (p1 = parse_rule_set (p))
    p = white_space (p1);

  if (p1 = or_pat (pu_s, cv, iv, p))
    return p1;
  else
    return punit_list (pu_s, cv, iv, p);
}

/**
* @brief Parses a list of patterns units where disjunctions are allowed.
*
* @param pu_s Will hold parsed pattern units.
* @param cv Will hold characters of exact and fuzzy patterns.
* @param iv Will hold a potential weight table.
* @param p The pointer to the first character of the first pattern unit.
*
* @return The pointer to the character immediately after the list.
*         Returns NULL if it is not successfully parsed.
*/
char *
punit_cons_list (pu_s, cv, iv, p)
     punit_t **pu_s;
     char **cv;
     int **iv;
     char *p;
{
  punit_t *last, *next;
  char *p1;

  last = NULL;

  p = white_space (p);
  next = *pu_s;
  while (p1 = punit_cons (pu_s, cv, iv, p)) {
    next->prev_punit = last;
    if (last)
      last->nxt_punit = next;

    while (next->nxt_punit)
      next = next->nxt_punit;

    last = next;
    next = *pu_s;
    p = white_space (p1);
  }
  if (!last)
    return NULL;
  else {
    last->nxt_punit = NULL;
    return p;
  }
}

/**
* @brief Parses a string of character into parse units.
*
* @param pu_s Will hold parsed pattern units.
* @param cv Will hold characters of exact and fuzzy patterns.
* @param iv Will hold a potential weight table.
* @param p The pointer to the first character of the first pattern unit.
*
* @return The first parse unit parsed.
*         Returns NULL if the string is not successfully parsed.
*/
punit_t *
parser (pu_s, cv, iv, line)
     punit_t **pu_s;
     char **cv;
     int **iv;
     char *line;
{
  punit_t *first;

  first = *pu_s;
  first->prev_punit = NULL;

  line = punit_cons_list (pu_s, cv, iv, white_space (line));
  if (line && (line = white_space (line)) && *line == '\0')  /* read all of it */
    return first;
  else
    return NULL;
}

/**
* @brief Sets the achored value of pu to 1 and runs set_anchored_on recursively
*        on its two OR 'children if the type is OR_PUNIT or simply on the next
*        if not.
* @param pu The parse unit to run the function on.
*/
void
set_anchors_on (pu)
     struct punit *pu;
{
  pu->anchored = 1;
  if (pu->type == OR_PUNIT) {
    set_anchors_on (pu->info.or.or1);
    set_anchors_on (pu->info.or.or2);
  }
  if (pu = pu->nxt_punit) {
    set_anchors_on (pu);
  }
}

/**
* @brief Sets the achored value of pu to 0 and runs set_anchored_on recursively
*        on its two OR 'children if the type is OR_PUNIT or simply on the next
*        if not.
* @param pu The parse unit to run the function on.
*/
void
set_anchors (pu)
     struct punit *pu;
{
  pu->anchored = 0;
  if (pu->type == OR_PUNIT) {
    set_anchors (pu->info.or.or1);
    set_anchors (pu->info.or.or2);
  }

  if (pu = pu->nxt_punit) {
    set_anchors_on (pu);
  }
}


/**
* @brief Set sequence type to DNA, and parses the line.
*
* @param line The line of pattern units to be parsed.
*
* @return Maximum matching size of the pattern units parsed.
*         Returns 0 if the pattern cannot be parsed.
*/
int
parse_dna_cmd (line)
     char *line;
{
  punit_sequence_type = DNA;
  return (parse_cmd (line));
}

/**
* @brief Set sequence type to PEPTIDE, and parses the line.
*
* @param line The line of pattern units to be parsed.
*
* @return Maximum matching size of the pattern units parsed.
*         Returns 0 if the pattern cannot be parsed.
*/
int
parse_peptide_cmd (line)
     char *line;
{
  punit_sequence_type = PEPTIDE;
  return (parse_cmd (line));
}

/**
* @brief Initializes the arrays holding the parse units, weight table and
*        character code tables, Then parses the given line.
*
* @param line The line of pattern units to be parsed.
*
* @return Maximum matching size of the pattern units parsed.
*         Returns 0 if the pattern cannot be parsed.
*/
int
parse_cmd (line)
     char *line;
{
  static punit_t pu_s[MAX_PUNITS];
  static char cv[MAX_TOTAL_CODES];
  static int iv[21 * MAX_WEIGHTS / 2];

  return (parse_cmd_ns (line, pu_s, cv, iv));
}

/**
* @brief Gets the current state, sets the sequence type to DNA,
*        and parses the given line.
*
* @param line The line of pattern units to be parsed.
* @param st The state to use.
*
* @return Maximum matching size of the pattern units parsed.
*         Returns 0 if the pattern cannot be parsed.
*/
int
parse_dna_cmd_ns (line, st)
     char *line;
     struct state **st;
{
  struct state *st1;
  int rc;

  st1 = get_state ();
  *st = st1;

  punit_sequence_type = DNA;
  rc = parse_cmd_ns (line, st1->pu_s, st1->cv, st1->iv);
  st1->ad_pu_s = ad_pu_s;
  return rc;
}

/**
* @brief Gets the current state, sets the sequence type to PEPTIDE,
*        and parses the given line.
*
* @param line The line of pattern units to be parsed.
* @param st The state to use.
*
* @return Maximum matching size of the pattern units parsed.
*         Returns 0 if the pattern cannot be parsed.
*/
int
parse_peptide_cmd_ns (line, st)
     char *line;
     struct state **st;
{
  struct state *st1;
  int rc;

  st1 = get_state ();
  *st = st1;

  punit_sequence_type = PEPTIDE;
  rc = parse_cmd_ns (line, st1->pu_s, st1->cv, st1->iv);
  st1->ad_pu_s = ad_pu_s;
  return rc;
}

/**
* @brief Parses a string of character into parse units.
*
* @param line The line of pattern units to be parsed.
* @param pu_s Will hold parsed pattern units.
* @param cv Will hold characters of exact and fuzzy patterns.
* @param iv Will hold a potential weight table.
*
* @return Maximum matching size of the pattern units parsed.
*         Returns 0 if the pattern cannot be parsed.
*/
int
parse_cmd_ns (line, pu_s, cv, iv)
     char *line;
     punit_t pu_s[MAX_PUNITS];
     char cv[MAX_TOTAL_CODES];
     int iv[21 * MAX_WEIGHTS / 2];
{
  char *ad_cv;
  int *ad_iv;

  int i;

  if (!initialized)
    build_conversion_tables ();

  for (i = 0; i < MAX_NAMES; i++)
    names[i] = NULL;

  ad_iv = iv;
  ad_cv = cv;
  ad_pu_s = (struct punit *) pu_s;

  if (((int) strlen (line)) >= ((int) MAX_SOUGHT_CHARS)) {
    ad_pu_s = NULL;
    return 0;
  }
  else {
    if (parser (&ad_pu_s, &ad_cv, &ad_iv, line) == NULL) {
/*          fprintf(stderr,"Sorry, could not parse %s\n",line);   */
      ad_pu_s = NULL;
      return 0;
    }
    else {
      set_anchors (pu_s);
      ad_pu_s = (struct punit *) pu_s;
      return max_mats (pu_s);
    }
  }
}

/**
* @brief Given a parse unit, returns the maximum matching lenght of this.
*
* @param pu Pattern unit to be calculated maximum matching length of.
*
* @return The maximum matching length of the parse unit.
*/
int
max_mat (pu)
     struct punit *pu;
{
  struct punit *p1;
  long x;
  int i;

  switch (pu->type) {
  case EXACT_PUNIT:
    return pu->info.exact.len;
    break;

  case LLIM_PUNIT:
    return 0;
    break;

  case RANGE_PUNIT:
    return pu->info.range.min + pu->info.range.width;
    break;

  case ANY_PUNIT:
    return 1;
    break;

  case COMPL_PUNIT:
    p1 = names[pu->info.compl.of];
    return (max_mat (p1) + pu->info.compl.ins);
    break;

  case REPEAT_PUNIT:
  case INV_REP_PUNIT:
    p1 = names[pu->info.repeat.of];
    return (max_mat (p1) + pu->info.repeat.ins);
    break;

  case SIM_PUNIT:
    return (pu->info.sim.len + pu->info.sim.ins);
    break;

  case WEIGHT_PUNIT:
    return pu->info.wvec.len;
    break;

  case OR_PUNIT:
    return MAX (max_mat (pu->info.or.or1), max_mat (pu->info.or.or2));
    break;
  }
  return (0);
}

/**
* @brief Given a parse unit list, returns the maximum matching lenght of these.
*
* @param pu Pattern unit list to be calculated maximum matching length of.
*
* @return The maximum matching length of the parse unit list.
*/
int
max_mats (pu)
     struct punit *pu;
{
  int sum = 0;

  while (pu) {
    sum += max_mat (pu);
    pu = pu->nxt_punit;
  }
  return sum;
}

/**
* @brief Builds a conversion table, this table will be used to convert DNA
*        Characters to bit versions.
*
* @return Always 0.
*/
int
build_conversion_tables () {
  int the_char;

  for (the_char = 0; the_char < 256; the_char++) {
    switch (tolower (the_char)) {
    case 'a':
      punit_to_code[the_char] = A_BIT;
      break;
    case 'c':
      punit_to_code[the_char] = C_BIT;
      break;
    case 'g':
      punit_to_code[the_char] = G_BIT;
      break;
    case 't':
      punit_to_code[the_char] = T_BIT;
      break;
    case 'u':
      punit_to_code[the_char] = T_BIT;
      break;
    case 'm':
      punit_to_code[the_char] = (A_BIT | C_BIT);
      break;
    case 'r':
      punit_to_code[the_char] = (A_BIT | G_BIT);
      break;
    case 'w':
      punit_to_code[the_char] = (A_BIT | T_BIT);
      break;
    case 's':
      punit_to_code[the_char] = (C_BIT | G_BIT);
      break;
    case 'y':
      punit_to_code[the_char] = (C_BIT | T_BIT);
      break;
    case 'k':
      punit_to_code[the_char] = (G_BIT | T_BIT);
      break;
    case 'b':
      punit_to_code[the_char] = (C_BIT | G_BIT | T_BIT);
      break;
    case 'd':
      punit_to_code[the_char] = (A_BIT | G_BIT | T_BIT);
      break;
    case 'h':
      punit_to_code[the_char] = (A_BIT | C_BIT | T_BIT);
      break;
    case 'v':
      punit_to_code[the_char] = (A_BIT | C_BIT | G_BIT);
      break;
    case 'n':
      punit_to_code[the_char] = (A_BIT | C_BIT | G_BIT | T_BIT);
      break;
    default:
      punit_to_code[the_char] = 0;
      break;
    }
    if (punit_to_code[the_char] & A_BIT)
      punit_to_code[the_char] |= T_BIT << 4;
    if (punit_to_code[the_char] & C_BIT)
      punit_to_code[the_char] |= G_BIT << 4;
    if (punit_to_code[the_char] & G_BIT)
      punit_to_code[the_char] |= C_BIT << 4;
    if (punit_to_code[the_char] & T_BIT)
      punit_to_code[the_char] |= A_BIT << 4;
  }

  for (the_char = 0; the_char < 256; the_char++) {
    switch (the_char & 15) {
    case A_BIT:
      code_to_punit[the_char] = 'A';
      break;

    case C_BIT:
      code_to_punit[the_char] = 'C';
      break;

    case G_BIT:
      code_to_punit[the_char] = 'G';
      break;

    case T_BIT:
      code_to_punit[the_char] = 'T';
      break;

    case A_BIT + C_BIT:
      code_to_punit[the_char] = 'M';
      break;

    case A_BIT + G_BIT:
      code_to_punit[the_char] = 'R';
      break;

    case A_BIT + T_BIT:
      code_to_punit[the_char] = 'W';
      break;

    case C_BIT + G_BIT:
      code_to_punit[the_char] = 'S';
      break;

    case C_BIT + T_BIT:
      code_to_punit[the_char] = 'Y';
      break;

    case G_BIT + T_BIT:
      code_to_punit[the_char] = 'K';
      break;

    case C_BIT + G_BIT + T_BIT:
      code_to_punit[the_char] = 'B';
      break;

    case A_BIT + G_BIT + T_BIT:
      code_to_punit[the_char] = 'D';
      break;

    case A_BIT + C_BIT + T_BIT:
      code_to_punit[the_char] = 'H';
      break;

    case A_BIT + C_BIT + G_BIT:
      code_to_punit[the_char] = 'V';
      break;

    case A_BIT + C_BIT + G_BIT + T_BIT:
      code_to_punit[the_char] = 'N';
      break;
    }
  }
  initialized = 1;
  return (0);
}

/**
* @brief Takes a string of dna characters and converts them to their
*        bit versions
*
* @param in String of DNA characters to be converted.
* @param out Upon completion holds the string of DNA bit versions.
*
* @return Always 0.
*/
int
comp_data (in, out)
     char *in, *out;
{
  while (*in)
    *(out++) = punit_to_code[*(in++)];
  *out = '\0';
  return (0);
}

/**
 * @brief Finds the best pivot pattern unit for pattern order optimisation.
 *
 * @param max_before_length This will upon successful execution hold the length
 *                          of max bases before the pivot pattern unit.
 * @param max_after_length This will upon successful execution hold the length
 *                         of max bases  the pivot pattern unit.
 * @param index This will upon successful execution hold the index where the
 *              best pivot pattern unit is.
 *
 * @return The index of the best pivot pattern unit.
 */
int
find_pivot_punit(max_before_length, max_after_length, index)
int *max_before_length, *max_after_length , *index;
{
  punit_t *current;
  *index = 0;
  *max_before_length = 0;
  int tmp_index = *index;
  int tmp_max_length = 0;
  int exact_length = 0;
  int good_choice = 0;
  int smallest_k = 100;
  int del_ins, k;
  for (current = ad_pu_s; current; current = current->nxt_punit) {
    switch (current->type) {
      case EXACT_PUNIT:
        if (current->info.exact.len > exact_length || !good_choice) {
          *index = tmp_index;
          *max_before_length = tmp_max_length;
          exact_length = current->info.exact.len;
          good_choice = 1;
        }
        break;
      case SIM_PUNIT:
        k = current->info.sim.ins + current->info.sim.del;
        //Only mismatches, insertions or deletions.
        if ((current->info.sim.mis & !current->info.sim.del & !current->info.sim.ins) ||
            (!current->info.sim.mis & current->info.sim.del & !current->info.sim.ins) ||
            (!current->info.sim.mis & !current->info.sim.del & current->info.sim.ins) &&
            ((current->info.sim.len - k) > exact_length )){
          *index = tmp_index;
          *max_before_length = tmp_max_length;
          exact_length = current->info.sim.len - k;
          good_choice = 1;
        }
        break;
    }
    tmp_max_length += max_mat(current);
    tmp_index += 1;
  }
  *max_after_length = tmp_max_length - *max_before_length;
  return *index;
}
