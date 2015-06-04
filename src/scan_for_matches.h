#ifndef SCAN_FOR_MATCHES_H
#define SCAN_FOR_MATCHES_H

struct match {
  char* match;
  long int start;
  long int end;
  struct match* next;
};

typedef struct match match_t;

match_t*
find_matches(char *fasta_file, char *line, char *ignore_file, int protein, 
             int complements, int show_overlaps, int max_hits, int stop_after, int verbose);


//FUNCTION TO EXPORT GOES IN HERE

#endif
