#include "memory.h"
#include "huffman.h"

#define USE_HUFFMAN 1

#define OPT_DEAMORTIZE 1

typedef struct {
  i64 freq;
  int rank;
} rank_code;

typedef struct {
  i64 n; // number of bytes
  int bs; // block size
  int sbs; // small-block size
//  int cmin, cmax; // min/max of size (in bytes) of compressed blocks
  darray *da;

  u16 opt;

  int logn; // number of bytes to represent n
  int k; // number of characters in super-alphabet
  int sigma; // size of super-alphabet
#if USE_HUFFMAN
  i64 *freq[3];
//  uchar *freq[3];
  Huffman2 *huf[3];
#else
  rank_code *code[3];
#endif

  i64 phase;
  i64 step; // current rewriting position
  i64 bu; // number of blocks to rewrite for each update

  i64 counter; // for deamortization
  i64 nbytes; // number of bytes written
//  double divergence;

} CRAM;

CRAM *cram_initialize(char *filename, int bs, int sbs, int bu, u16 opt);
i64 cram_usedmemory(CRAM *cram);
void cram_read(CRAM *cram, i64 s, i64 t, uchar *buf);
void cram_write(CRAM *cram, i64 s, i64 t, uchar *buf);

void cram_save(CRAM *cram, uchar *filename, uchar *memfilename);
CRAM *cram_load(uchar *filename, uchar *memfilename);
void cram_free(CRAM *cram);
