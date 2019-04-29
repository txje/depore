#include "aln.h"

// query along y-axis, target along x
result align_full_matrix(char* query, char* target, int qlen, int tlen, charvec *path, int semilocal, int match, int mismatch, int ins, int del) {

  if(tlen == 0 || qlen == 0) {
    result res;
    res.failed = 1;
    // other fields are unset and unreliable
    return res;
  }

  // must be dynamically allocated to prevent stack overflows (on *some* systems)
  int* score_matrix[qlen];
  char* direction_matrix[qlen];
  int i;
  for(i = 0; i < qlen; i++) {
    score_matrix[i] = (int*)malloc(tlen * sizeof(int));
    direction_matrix[i] = (char*)malloc(tlen * sizeof(char));
  }

  int x, y;
  int max_x = 0, max_y = 0; // used to compute semilocal path (can pick best score that doesn't hit the end)

  for(y = 0; y < qlen; y++) {
    for(x = 0; x < tlen; x++) {
      int match_score, ins_score, del_score;

      // match
      if(query[y] == target[x]) {
        match_score = match;
      } else {
        match_score = mismatch;
      }
      if(y > 0 && x > 0) {
        match_score = score_matrix[y-1][x-1] + match_score;
      }

      // ins
      if(y == 0) {
        ins_score = ins;
      } else {
        ins_score = score_matrix[y-1][x] + ins;
      }

      // del
      if(x == 0) {
        del_score = del;
      } else {
        del_score = score_matrix[y][x-1] + del;
      }

      // compare
      if(match_score >= ins_score && match_score >= del_score) {
        score_matrix[y][x] = match_score;
        if(query[y] == target[x]) {
          direction_matrix[y][x] = MATCH;
        } else {
          direction_matrix[y][x] = MISMATCH;
        }
      } else if(ins_score >= del_score) {
        score_matrix[y][x] = ins_score;
        direction_matrix[y][x] = INS;
      } else {
        score_matrix[y][x] = del_score;
        direction_matrix[y][x] = DEL;
      }

      if(semilocal && score_matrix[y][x] > score_matrix[max_y][max_x]) {
        max_y = y;
        max_x = x;
      }
    }
  }

  // compute maximum score position
  if(!semilocal) { // global
    max_x = 0;
    max_y = qlen - 1;
    for(x = 1; x < tlen; x++) { // check last row
      if(score_matrix[qlen - 1][x] > score_matrix[max_y][max_x]) {
        max_y = qlen - 1;
        max_x = x;
      }
    }
    for(y = 0; y < qlen; y++) { // check last column
      if(score_matrix[y][tlen - 1] > score_matrix[max_y][max_x]) {
        max_y = y;
        max_x = tlen - 1;
      }
    }
  }

  x = max_x;
  y = max_y;
  int lastx, lasty;
  while(x >= 0 && y >= 0) {
    if(path != NULL)
      kv_push(char, *path, direction_matrix[y][x]);
    lastx = x;
    lasty = y;
    if(direction_matrix[y][x] == MATCH || direction_matrix[y][x] == MISMATCH) {
      x--;
      y--;
    } else if(direction_matrix[y][x] == INS) {
      y--;
    } else if(direction_matrix[y][x] == DEL) {
      x--;
    }
  }
  
  //invert the path, since we added it in reverse order
  if(path != NULL) {
    char tmp;
    for(i = 0; i < kv_size(*path)/2; i++) {
      tmp = kv_A(*path, i);
      kv_A(*path, i) = kv_A(*path, kv_size(*path)-1-i);
      kv_A(*path, kv_size(*path)-1-i) = tmp;
    }
  }

  result res;
  res.score = score_matrix[max_y][max_x];
  res.qstart = lasty;
  res.qend = max_y;
  res.tstart = lastx;
  res.tend = max_x;
  // end positions are INCLUSIVE

  for(i = 0; i < qlen; i++) {
    free(score_matrix[i]);
    free(direction_matrix[i]);
  }

  return res;
}
