#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>
#include <string.h>

#include "../incl/klib/kseq.h"
#include "ssw.h"
#include "adapters.h"
#include "ssw_write.c"
#include "aln.h"

#define SEQPREFIX 50
#define ADAPTERLEN 100
#define SEQPOSTFIX 50

// have to reorder params to make this work with kseq
int fileread(FILE* f, char* buffer, int size) {
  return fread(buffer, 1, size, f);
}

// init kseq struct
KSEQ_INIT(FILE*, fileread);

void version() {
  printf("poretrim version 0.1\n");
}

void usage() {
  printf("Usage: poretrim [options]\n");
  printf("Options:\n");
  printf("  -r: Read FASTA/Q[.gz]\n");
  printf("  -o: Output prefix\n");
  printf("  -v, --verbose: verbose\n");
  printf("  -h, --help: show this\n");
  printf("  --version: show version information\n");
}

static struct option long_options[] = {
// if these are the same as a single-character option, put that character in the 4th field instead of 0
  { "verbose", no_argument, 0, 'v' },
  { "help",    no_argument, 0, 'h' },
  { "version", no_argument, 0, 0 },
  { 0, 0, 0, 0}
};

/* This table is used to transform nucleotide letters into numbers. */
static const int8_t nt_table[128] = {
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


char* rc(char* s, int l, char compl[256]) {
  char* r = malloc((l+1) * sizeof(char));
  int i;
  for(i = 0; i < l; i++)
    r[i] = compl[s[l-1-i]];
  r[i] = '\0';
  return r;
}


// porechop default: match = 3, mismatch = -6, gap open = -5, gap extend = -2
s_align* align(char* ref_seq, char* read_seq, int ref_len, int read_len) {
  int32_t l, m, k, match = 3, mismatch = -5, gap_open = -1, gap_extension = -1;  // default parameters for genome sequence alignment
  s_profile* profile;
  int8_t* num = (int8_t*)malloc(read_len);  // the read sequence represented in numbers
  int8_t* ref_num = (int8_t*)malloc(ref_len);
  s_align* result;

  // initialize scoring matrix for genome sequences
  //  A  C  G  T  N (or other ambiguous code)
  //  3 -6 -6 -6  0 A
  // -6  3 -6 -6  0 C
  // -6 -6  3 -6  0 G
  // -6 -6 -6  3  0 T
  //  0  0  0  0  0 N (or other ambiguous code)
  int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
  for (l = k = 0; l < 4; ++l) {
    for (m = 0; m < 4; ++m)
      mat[k++] = l == m ? match : mismatch;
    mat[k++] = 0; // ambiguous base: no penalty
  }
  for (m = 0; m < 5; ++m) mat[k++] = 0;

  for (m = 0; m < read_len; ++m) num[m] = nt_table[(int)read_seq[m]];
  profile = ssw_init(num, read_len, mat, 5, 2);
  for (m = 0; m < ref_len; ++m) ref_num[m] = nt_table[(int)ref_seq[m]];

  // Only the 8 bit of the flag is set. ssw_align will always return the best alignment beginning position and cigar.
  result = ssw_align(profile, ref_num, ref_len, gap_open, gap_extension, 1, 0, 0, read_len);

  init_destroy(profile);
  free(mat);
  free(ref_num);
  free(num);
  return result;
}


int main(int argc, char *argv[]) {
  int verbose = 0;
  char* read_fasta;
  char* output_prefix;

  int opt, long_idx;
  opterr = 0;
  while ((opt = getopt_long(argc, argv, "r:o:vh", long_options, &long_idx)) != -1) {
    switch (opt) {
      case 'r':
        read_fasta = optarg;
        break;
      case 'o':
        output_prefix = optarg;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'h':
        usage();
        return 0;
        break;
      case '?':
        if (optopt == 'r' || optopt == 'o')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      case 0:
        // as long as all the long arguments have characters too, I don't think this section will be used
        if (long_idx == 0) verbose = 1; // --verbose
        else if (long_idx == 1) {usage(); return 0;} // --help
        else if (long_idx == 2) {version(); return 0;} // --version
        break;
      default:
        usage();
        return 1;
    }
  }

  if(read_fasta == NULL) {
    fprintf(stderr, "-r read FASTA is required\n");
    usage();
    return 1;
  }
  if(output_prefix == NULL) {
    fprintf(stderr, "-o output prefix is required\n");
    usage();
    return 1;
  }

  // load FASTA file

  FILE* fp;
  kseq_t* seq;
  int i, l;
  char *adapter_seq, *rc_adapter_seq;
  //s_align* a;

  // set up complement array
  char compl[256];
  for(i = 0; i < 256; i++) {
    compl[i] = 'N';
  }
  compl[65] = 'T';
  compl[67] = 'G';
  compl[71] = 'C';
  compl[84] = 'A';
  compl[97] = 't';
  compl[99] = 'g';
  compl[103] = 'c';
  compl[116] = 'a';

  fp = fopen(read_fasta, "r");
  seq = kseq_init(fp);
  printf("Reading fasta file: %s\n", read_fasta);

  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    if(verbose) {
      fprintf(stderr, "Reading %s (%i bp).\n", seq->name.s, l);
    }

    // align read to known adapters, barcodes
    //
    /*
    adapter_seq = malloc(sizeof(char) * (RBK004_pre_l + 1));
    for(i = 0; i < RBK004_pre_l; i++) adapter_seq[i] = RBK004_pre[i];
    adapter_seq[i] = '\0';
    a = align(seq->seq.s, adapter_seq, l, RBK004_pre_l);
    fprintf(stderr, "RBK004_pre (%d:%d) aligned to %s (%d:%d) of %d with score %u\n", a->read_begin1, a->read_end1, seq->name.s, a->ref_begin1, a->ref_end1, l, a->score1);
    fprintf(stderr, "  cigar: ");
    for(i = 0; i < a->cigarLen; i++)
      fprintf(stderr, "%u%c ", cigar_int_to_len(a->cigar[i]), cigar_int_to_op(a->cigar[i]));
    fprintf(stderr, "\n");
    ssw_write(a, seq->seq.s, adapter_seq, nt_table);
    align_destroy(a);
    free(adapter_seq);

    adapter_seq = malloc(sizeof(char) * (RBK004_post_l + 1));
    for(i = 0; i < RBK004_post_l; i++) adapter_seq[i] = RBK004_post[i];
    adapter_seq[i] = '\0';
    a = align(seq->seq.s, adapter_seq, l, RBK004_post_l);
    fprintf(stderr, "RBK004_post (%d:%d) aligned to %s (%d:%d) of %d with score %u\n", a->read_begin1, a->read_end1, seq->name.s, a->ref_begin1, a->ref_end1, l, a->score1);
    fprintf(stderr, "  cigar: ");
    for(i = 0; i < a->cigarLen; i++)
      fprintf(stderr, "%u%c ", cigar_int_to_len(a->cigar[i]), cigar_int_to_op(a->cigar[i]));
    fprintf(stderr, "\n");
    ssw_write(a, seq->seq.s, adapter_seq, nt_table);
    align_destroy(a);
    free(adapter_seq);
    */

    int bc;
    int scores[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
    int positions[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
    int max = 0;
    int delta = 0;
    for(bc = 1; bc <= 12; bc++) {

      adapter_seq = malloc(sizeof(char) * (RBK004_pre_l + RBK004_bc_l + RBK004_post_l + 1));
      for(i = 0; i < RBK004_pre_l; i++) adapter_seq[i] = RBK004_pre[i];
      for(i = 0; i < RBK004_bc_l; i++) adapter_seq[i+RBK004_pre_l] = RBK004_bcs[bc-1][i];
      for(i = 0; i < RBK004_post_l; i++) adapter_seq[i+RBK004_pre_l+RBK004_bc_l] = RBK004_post[i];
      adapter_seq[RBK004_pre_l + RBK004_bc_l + RBK004_post_l] = '\0';
      //rc_adapter_seq = rc(adapter_seq, RBK004_pre_l + RBK004_bc_l + RBK004_post_l, compl);

      /*
      a = align(seq->seq.s, adapter_seq, l, RBK004_pre_l+RBK004_bc_l+RBK004_post_l);
      fprintf(stderr, "RBK004_bcs[bc-1](w/pre,post) (%d:%d) aligned to %s (%d:%d) of %d with score %u\n", a->read_begin1, a->read_end1, seq->name.s, a->ref_begin1, a->ref_end1, l, a->score1);
      fprintf(stderr, "  cigar: ");
      for(i = 0; i < a->cigarLen; i++)
        fprintf(stderr, "%u%c ", cigar_int_to_len(a->cigar[i]), cigar_int_to_op(a->cigar[i]));
      fprintf(stderr, "\n");
      ssw_write(a, seq->seq.s, adapter_seq, nt_table);
      align_destroy(a);
      */

      charvec path;
      kv_init(path);
      //result r = align_full_matrix(adapter_seq, seq->seq.s, RBK004_pre_l+RBK004_bc_l+RBK004_post_l, l, &path, 0, 1, -1, -1, -1);
      result r = align_full_matrix(adapter_seq, seq->seq.s, RBK004_pre_l+RBK004_bc_l+RBK004_post_l, l, &path, 0, 1, -2, -1, -1); // semilocal
      if(r.failed) {
        fprintf(stderr, "alignment failed\n");
        continue;
      }
      scores[bc-1] = r.score;
      positions[bc-1] = r.tstart;
      if(r.score > max) {
        delta = r.score - max;
        max = r.score;
      } else if(max - r.score < delta) {
        delta = max - r.score;
      }
      //fprintf(stderr, "RBK004 barcode %d (%d:%d) aligned to %s (%d:%d) of %d with score %u\n", bc, r.qstart, r.qend, seq->name.s, r.tstart, r.tend, l, r.score);
      int q = r.qstart;
      int t = r.tstart;
      char* qstr = malloc((kv_size(path) + 1) * sizeof(char));
      char* tstr = malloc((kv_size(path) + 1) * sizeof(char));
      char* astr = malloc((kv_size(path) + 1) * sizeof(char));
      qstr[kv_size(path)] = '\0';
      tstr[kv_size(path)] = '\0';
      astr[kv_size(path)] = '\0';
      //fprintf(stderr, "  cigar (%d): ", kv_size(path));
      for(i = 0; i < kv_size(path); i++) {
        //fprintf(stderr, "%c", (char)kv_A(path, i));
        if(kv_A(path, i) == 'M') {
          qstr[i] = adapter_seq[q++];
          tstr[i] = seq->seq.s[t++];
          astr[i] = '|';
        } else if(kv_A(path, i) == 'X') {
          qstr[i] = adapter_seq[q++];
          tstr[i] = seq->seq.s[t++];
          astr[i] = '*';
        } else if(kv_A(path, i) == 'I') {
          qstr[i] = adapter_seq[q++];
          tstr[i] = '-';
          astr[i] = '*';
        } else if(kv_A(path, i) == 'D') {
          qstr[i] = '-';
          tstr[i] = seq->seq.s[t++];
          astr[i] = '*';
        }
      }
      //fprintf(stderr, "\n");
      /*
      fprintf(stderr, "q: %s\n", qstr);
      fprintf(stderr, "   %s\n", astr);
      fprintf(stderr, "t: %s\n", tstr);
      */

      free(adapter_seq);
    }
    fprintf(stderr, "%s\n", seq->name.s);
    fprintf(stderr, "score (position): ");
    for(i = 0; i < 12; i++) {
      if(i > 0) {
        fprintf(stderr, ", ");
      }
      fprintf(stderr, "%d (%d)", scores[i], positions[i]);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "max: %d, delta: %d\n", max, delta);
    
    fprintf(stderr, "%s%s", pad, prefix);

    /*!	@typedef	structure of the alignment result
      typedef struct {
        uint16_t score1;
        uint16_t score2;
        int32_t ref_begin1;
        int32_t ref_end1;
        int32_t	read_begin1;
        int32_t read_end1;
        int32_t ref_end2;
        uint32_t* cigar;
        int32_t cigarLen;
      } s_align;
    */

    // automatically detect experiment type (expected adapters, barcodes)
    // demultiplex (liberally), trim adapters (and barcodes), split reads with middle adapters
  }

  kseq_destroy(seq);
  fclose(fp);
}
