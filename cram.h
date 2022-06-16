#include "memory.h"
#include "huffman.h"
#include "buddy.h"

typedef struct {
  i64 n; // number of bytes
  int bs; // block size
  int sbs; // small-block size
  darray *da;
  

  int paddingLen;
  
  int k; // number of characters in super-alphabet
  int sigma; // size of super-alphabet
  i64 *freq[3];
  
  Huffman *huf[3];

  int currentTreeIndex;

  i64 step; // current rewriting position
  i64 bu; // number of blocks to rewrite for each update
  i64 nb;    

  //additional modules
  CodeAllocator *codeAllocator;
  AdditionalCodeTable *additionalCodeTable[2];
  i64 changeMinFreq[2];
  //u64 *encodedTree;
  int *encodedTree;
  i64 *pivotFreq;
  int canAddCode;
  u64 *currentCode[2];
  int *currentLen[2];
  int changeMultiple;
  int writeBlockUsingSaved;
  u64 newCode[512];
  int newCodeLen[512];
  int newCodeCh[512];
  int newCodeIdx;
} CRAM;

CRAM *cram_initialize(char *filename, int bs, int bu,int mode,int changeMultiple);
i64 cram_usedmemory(CRAM *cram);
void cram_read(CRAM *cram, i64 s, i64 t, uchar *buf);
void cram_write(CRAM *cram,i64 b,uchar *buf);
int cram_read_block(CRAM *cram, i64 b, uchar *buf);
void cram_save(CRAM *cram, uchar *filename, uchar *memfilename);
CRAM *cram_load(uchar *filename, uchar *memfilename);
void cram_free(CRAM *cram);
