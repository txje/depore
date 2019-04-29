#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../incl/klib/kvec.h" // C dynamic vector

#ifndef ALN_H
#define ALN_H

typedef kvec_t(char) charvec;

typedef struct aln_result {
  int score;
  int qstart;
  int qend;
  int tstart;
  int tend;
  int failed; // boolean flag
} result;

#define MATCH 'M'
#define INS 'I'
#define DEL 'D'
#define MISMATCH 'X'

#define LOW -2000000000 // almost the lowest 32-bit integer

result align_full_matrix(char* query, char* target, int qlen, int tlen, charvec *path, int semilocal, int match, int mismatch, int ins, int del);

#endif
